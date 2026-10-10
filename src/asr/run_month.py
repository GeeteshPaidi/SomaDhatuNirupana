"""
Per-month pipeline for the ASR rework (TRACKER 1.8 / 2.x).

  1. Read every L1 FITS header + raw counts for one month.
  2. Nightside background template: mean of counts / (EXPOSURE * n_SCD/16)
     over all valid nightside exposures (SOLARANG >= 90, SAT_ALT > 0).
  3. Energy calibration (gain, offset, sigma0): joint Gaussian fit of the
     Mg/Al/Si and Ca K lines in the stacked, background-subtracted dayside
     spectrum (all-16-SCD exposures), with sigma(E)^2 = sigma0^2 + F w E.
  4. Cash fit of every valid dayside spectrum (no denoising, no selection
     other than geometry validity) with src/asr/fitter.py.

Outputs (outputs/asr/months/<YYYY_MM>/): calib.json, background.npy,
fits.parquet (one row per dayside spectrum: header metadata, fitted
parameters, uncertainties, GoF, quality flag).
"""
import argparse
import glob
import json
import os
import sys
import time
from multiprocessing import Pool

import numpy as np
import pandas as pd
import xraylib
from astropy.io import fits
from scipy import optimize

sys.path.insert(0, os.path.dirname(__file__))
import fitter as ft  # noqa: E402
import solar  # noqa: E402
import xrf_model as xm  # noqa: E402

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
N_CH = 2048


def read_one(path):
    try:
        with fits.open(path, memmap=False) as h:
            hd, d = h[1].header, h[1].data
            counts = np.zeros(N_CH, dtype=np.float32)
            counts[np.clip(d["CHANNEL"].astype(int), 0, N_CH - 1)] = d["COUNTS"]
            scd = [s for s in str(hd.get("SCD_USED", "")).split(",") if s.strip()]
            return {
                "source_file": os.path.basename(path),
                "counts": counts,
                "exposure": float(hd["EXPOSURE"]), "solarang": float(hd["SOLARANG"]),
                "sat_alt": float(hd["SAT_ALT"]), "n_scd": len(scd), "temp": float(hd["TEMP"]),
                "mid_utc": str(hd["MID_UTC"]), "sw_versn": float(hd["SW_VERSN"]),
                "lat": [float(hd[f"V{k}_LAT"]) for k in range(4)],
                "lon": [float(hd[f"V{k}_LON"]) for k in range(4)],
                "bore_lat": float(hd["BORE_LAT"]), "bore_lon": float(hd["BORE_LON"]),
            }
    except Exception:  # noqa: BLE001
        return None


R_MOON = 1737.4


def shadow_angle(sat_alt):
    """SOLARANG beyond which the spacecraft itself is in the lunar shadow."""
    return 90.0 + np.degrees(np.arccos(R_MOON / (R_MOON + max(sat_alt, 1.0))))


def background_templates(rows, hi_band=(8.0, 20.0), gain=0.0135, offset=0.0, margin=1.0):
    """Two nightside templates (counts / s / 16-SCD-equivalent):
      P  particle background: spacecraft in lunar shadow (SOLARANG > shadow + margin)
      T  solar-driven instrumental component (Al Ka + continuum seen only when the
         spacecraft is sunlit): mean sunlit-nightside spectrum minus P scaled to the
         same 8-20 keV rate; normalised to its own mean level (gamma = 1).
    """
    E = gain * np.arange(N_CH) + offset
    hb = (E >= hi_band[0]) & (E < hi_band[1])
    ok = [r for r in rows if r["solarang"] >= 90 and r["sat_alt"] > 0 and r["exposure"] > 0 and r["n_scd"] > 0]
    sh = [r for r in ok if r["solarang"] > shadow_angle(r["sat_alt"]) + margin]
    sl = [r for r in ok if r["solarang"] < shadow_angle(r["sat_alt"]) - margin]

    def mean_rate(rs):
        tot, t = np.zeros(N_CH), 0.0
        for r in rs:
            tot += r["counts"]
            t += r["exposure"] * r["n_scd"] / 16.0
        return tot / max(t, 1e-9), t

    P, tP = mean_rate(sh)
    S, tS = mean_rate(sl)
    T = S - P * (S[hb].sum() / max(P[hb].sum(), 1e-12)) if len(sl) else np.zeros(N_CH)
    info = {"n_night_shadow": len(sh), "t_night_shadow_s": tP, "n_night_sunlit": len(sl), "t_night_sunlit_s": tS,
            "sunlit_over_shadow_hiband": float(S[hb].sum() / max(P[hb].sum(), 1e-12)) if len(sl) else None}
    return P, T, info


CAL_LINES = [("Mg", xraylib.KA_LINE), ("Al", xraylib.KA_LINE), ("Al", xraylib.KB_LINE),
             ("Si", xraylib.KA_LINE), ("Si", xraylib.KB_LINE), ("Ca", xraylib.KA_LINE), ("Ca", xraylib.KB_LINE)]


def calibrate(net, gain0=0.0135, offset0=0.0):
    """Joint Gaussian-line fit to the stacked net spectrum -> gain, offset, sigma0."""
    E_l = np.array([xraylib.LineEnergy(xraylib.SymbolToAtomicNumber(e), l) for e, l in CAL_LINES])
    kb = {"Al": 1, "Si": 1, "Ca": 1}
    ka_idx = {e: i for i, (e, l) in enumerate(CAL_LINES) if l == xraylib.KA_LINE}
    rr = np.array([xraylib.RadRate(xraylib.SymbolToAtomicNumber(e), l) for e, l in CAL_LINES])
    win1 = np.arange(int(1.05 / gain0), int(2.10 / gain0))
    win2 = np.arange(int(3.25 / gain0), int(4.40 / gain0))
    ch = np.r_[win1, win2]
    y = net[ch]
    w = 1.0 / np.maximum(np.abs(y), 1.0)

    def model(p):
        gain, off, s0 = p[:3]
        amps = dict(zip(["Mg", "Al", "Si", "Ca"], p[3:7]))
        b = p[7:]
        E = gain * ch + off
        out = np.zeros_like(E)
        for i, (e, l) in enumerate(CAL_LINES):
            a = amps[e] * (rr[i] / rr[ka_idx[e]])
            s = np.sqrt(s0 ** 2 + xm.FANO * xm.W_PAIR * E_l[i])
            out += a * gain / (s * np.sqrt(2 * np.pi)) * np.exp(-0.5 * ((E - E_l[i]) / s) ** 2)
        lo = ch < win2[0]
        out[lo] += b[0] + b[1] * (E[lo] - 1.5)
        out[~lo] += b[2] + b[3] * (E[~lo] - 3.8)
        return out

    peak = lambda e: net[int(e / gain0)] * 0.06 * np.sqrt(2 * np.pi) / gain0  # noqa: E731
    p0 = np.r_[gain0, offset0, 0.045, peak(1.25), peak(1.49), peak(1.74), peak(3.69), 0, 0, 0, 0]
    res = optimize.least_squares(lambda p: (model(p) - y) * np.sqrt(w), p0,
                                 bounds=(np.r_[0.0130, -0.1, 0.01, np.zeros(4), -np.inf * np.ones(4)],
                                         np.r_[0.0142, 0.1, 0.15, np.full(4, np.inf), np.inf * np.ones(4)]))
    gain, off, s0 = res.x[:3]
    return {"gain": float(gain), "offset": float(off), "sigma0": float(s0),
            "fwhm_1.49keV_eV": float(2.3548e3 * np.sqrt(s0 ** 2 + xm.FANO * xm.W_PAIR * 1.4867)),
            "fwhm_5.9keV_eV": float(2.3548e3 * np.sqrt(s0 ** 2 + xm.FANO * xm.W_PAIR * 5.895)),
            "success": bool(res.success)}


_F = None


def make_fitter(calib, bkg, inst, kT=None):
    """kT=None -> CHIANTI table with per-spectrum GOES temperature;
    kT>0 -> legacy continuum-only spectrum at fixed kT (keV)."""
    det = xm.Detector(gain=calib["gain"], offset=calib["offset"], sigma0=calib["sigma0"])
    if kT is None:
        fm = xm.ForwardModel(secondary=False, e_max=15, de=0.005, solar=solar.SolarTable(), T_mk=8.0)
    else:
        fm = xm.ForwardModel(kT=kT, secondary=False, e_max=15, de=0.005)
    return ft.SpectrumFitter(fm, det, bkg, inst_rate=inst)


def _init(kT, calib, bkg, inst):
    global _F
    _F = make_fitter(calib, bkg, inst, kT)


def _fit(r):
    o = _F.fit(r["counts"], r["exposure"], r["n_scd"], r["solarang"], T_mk=r.get("T_fit"))
    o.pop("corr", None)
    o["source_file"] = r["source_file"]
    return o


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--month", default="2024/12")
    ap.add_argument("--kT", type=float, default=None,
                    help="legacy fixed continuum kT (keV); default: GOES T + CHIANTI")
    ap.add_argument("--T_scale", type=float, default=1.0,
                    help="sensitivity: multiply the GOES temperature by this factor")
    ap.add_argument("--T_fixed", type=float, default=None, help="sensitivity: fixed CHIANTI T (MK)")
    ap.add_argument("--workers", type=int, default=12)
    ap.add_argument("--limit", type=int, default=None, help="max dayside spectra to fit (random subset)")
    ap.add_argument("--out_root", default=os.path.join(ROOT, "outputs", "asr", "months"))
    args = ap.parse_args()
    tag = args.month.replace("/", "_")
    suffix = ""
    if args.kT is not None:
        suffix = f"_kT{args.kT}"
    elif args.T_fixed is not None:
        suffix = f"_Tfix{args.T_fixed}"
    elif args.T_scale != 1.0:
        suffix = f"_Tx{args.T_scale}"
    out_dir = os.path.join(args.out_root, tag + suffix)
    os.makedirs(out_dir, exist_ok=True)

    t0 = time.perf_counter()
    paths = sorted(glob.glob(os.path.join(ROOT, "data", "raw", *args.month.split("/"), "**", "*.fits"),
                        recursive=True))
    with Pool(args.workers) as p:
        rows = [r for r in p.imap(read_one, paths, chunksize=256) if r is not None]
    n_read = len(rows)
    print(f"{tag}: {len(paths)} files, {n_read} readable ({time.perf_counter()-t0:.0f}s)", flush=True)

    bkg, inst, binfo = background_templates(rows)
    np.save(os.path.join(out_dir, "background.npy"), bkg)
    np.save(os.path.join(out_dir, "instrumental.npy"), inst)
    valid = [r for r in rows if r["sat_alt"] > 0 and r["exposure"] > 0 and r["n_scd"] > 0]
    day = [r for r in valid if r["solarang"] < 90]
    stack_rows = [r for r in day if r["n_scd"] == 16]
    stack = np.sum([r["counts"] for r in stack_rows], axis=0)
    t_stack = sum(r["exposure"] for r in stack_rows)
    net = stack - bkg * t_stack
    calib = calibrate(net)
    calib.update(month=tag, n_files=len(paths), n_readable=n_read, n_valid_geom=len(valid),
                 n_dayside=len(day), **binfo,
                 n_stack=len(stack_rows), kT=args.kT)
    print(json.dumps(calib, indent=1), flush=True)
    np.savez(os.path.join(out_dir, "stack.npz"), stack=stack, t_stack=t_stack, bkg=bkg)

    if args.limit and len(day) > args.limit:
        rng = np.random.default_rng(0)
        day = [day[i] for i in sorted(rng.choice(len(day), args.limit, replace=False))]
    T, EM, XB = solar.goes_at([r["mid_utc"] for r in day])
    for r, t, em, xb in zip(day, T, EM, XB):
        r["goes_T_mk"], r["goes_EM"], r["goes_xrsb"] = float(t), float(em), float(xb)
        if args.kT is None:
            r["T_fit"] = args.T_fixed if args.T_fixed is not None else t * args.T_scale
    calib["frac_with_goes"] = float(np.mean(np.isfinite(T))) if len(T) else 0.0
    t1 = time.perf_counter()
    with Pool(args.workers, initializer=_init, initargs=(args.kT, calib, bkg, inst)) as p:
        res = list(p.imap(_fit, day, chunksize=16))
    fit_s = time.perf_counter() - t1
    meta = pd.DataFrame([{k: v for k, v in r.items() if k not in ("counts", "T_fit")} for r in day])
    df = meta.merge(pd.DataFrame(res), on="source_file")
    df["month"] = tag
    df.to_parquet(os.path.join(out_dir, "fits.parquet"), index=False)
    calib.update(n_fitted=len(df), fit_wall_s=fit_s, ms_per_fit_wall=fit_s / max(len(df), 1) * 1000)
    with open(os.path.join(out_dir, "calib.json"), "w") as f:
        json.dump(calib, f, indent=1)
    print(f"fitted {len(df)} in {fit_s:.0f}s; quality==0: {(df['quality']==0).mean():.3f}", flush=True)


if __name__ == "__main__":
    main()
