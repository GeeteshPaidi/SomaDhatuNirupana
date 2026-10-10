"""
Task 1.7 GATE: synthetic recovery test.

Generates CLASS-like spectra with the forward model (known calibration,
real nightside background template, Poisson noise), over a grid of signal
levels, random lunar compositions and incidence angles, then fits them.
Optionally generates with one solar temperature and fits with another to
measure the "no XSM" systematic.

Signal level is set by the Si Ka line counts in the fit window. The
scattered-continuum amplitude is set so that the continuum at 3 keV is
SC_FRAC_3KEV x the Si Ka peak height, matching the 2024-12 dayside stack
(continuum ~0.09 counts/ch/spectrum vs Si peak ~23 counts/ch/spectrum).

Outputs: outputs/asr/synth/recovery_<tag>.parquet + summary CSV.
"""
import argparse
import os
import sys
import time
from multiprocessing import Pool

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
import fitter as ft  # noqa: E402
import xrf_model as xm  # noqa: E402

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
CAL = dict(gain=0.01363, offset=0.017, sigma0=0.045)
SC_FRAC_3KEV = 0.004
GAMMA_TRUE = 0.5
FIXED = None  # None -> fitter default (Ti, Fe fixed); {} -> all free
# Lunar composition ranges (weight ratios to Si), from PLAN.md sanity ranges / Apollo soils.
RANGES = {"Mg": (0.10, 0.35), "Al": (0.20, 0.80), "Ca": (0.30, 0.60), "Ti": (0.005, 0.30), "Fe": (0.10, 0.80)}

_G = {}


def _fm(kT, secondary, use_solar):
    if use_solar:  # kT is a temperature in MK for the CHIANTI table
        import solar
        return xm.ForwardModel(secondary=secondary, e_max=15, de=0.005, solar=solar.SolarTable(), T_mk=kT)
    return xm.ForwardModel(kT=kT, secondary=secondary, e_max=15, de=0.005)


def _init(kT_gen, kT_fit, bkg_path, inst_path=None, use_solar=False):
    if bkg_path.endswith(".npz"):  # legacy single month-mean nightside template, no instrumental term
        d = np.load(bkg_path)
        bkg, inst = d["night"] / d["tn"], None
    else:
        bkg = np.load(bkg_path)
        inst = np.load(inst_path) if inst_path else None
    det = xm.Detector(**CAL)
    _G["gen"] = ft.SpectrumFitter(_fm(kT_gen, True, use_solar), det, bkg, inst_rate=inst)
    _G["fit"] = ft.SpectrumFitter(_fm(kT_fit, False, use_solar), det, bkg, inst_rate=inst,
                                  fixed_ratios=FIXED if FIXED is not None else None)


def truth_theta(G, ratios, s1, si_counts):
    """Return theta giving `si_counts` Si Ka counts and the stack-like scatter level."""
    r = np.array([ratios[e] for e in xm.FIT_ELEMENTS])
    th = np.r_[1.0, r, 0.0, 0.0, 0.0]
    I, comp = G.fm.primary_fast(ratios, s1)
    i_si = [i for i, ln in enumerate(G.fm.lines) if ln.el == "Si" and ln.label == "Ka"][0]
    th[0] = si_counts / (I[i_si] / G.i_norm)
    lines = G.source(th, s1)
    si_peak = lines.max()
    sc = G._scatter(comp, s1)
    k3 = np.argmin(np.abs(G.E - 3.0))
    th[6] = SC_FRAC_3KEV * si_peak / (th[0] * sc[k3])
    th[7] = 1.0
    th[8] = GAMMA_TRUE if G.use_inst else 0.0
    return th


def generate_counts(G, th, s1, exposure=8.0, n_scd=16, rng=None, secondary=True):
    """Expected counts using the generator model (with secondary fluorescence)."""
    r = dict(zip(xm.FIT_ELEMENTS, th[1:6]))
    if secondary:
        psi1 = np.degrees(np.arcsin(s1))
        I = G.fm.line_intensities(r, psi1)
    else:
        I, _ = G.fm.primary_fast(r, s1)
    comp = xm.composition_from_ratios(r)
    src = th[0] * G.eff * ((I / G.i_norm) @ G.R + th[6] * G._scatter(comp, s1))
    bk = G.background_counts(exposure, n_scd)
    m = src + th[7] * bk[0] + th[8] * bk[1]
    counts = np.zeros(G.det.n_channels)
    counts[G.sel] = rng.poisson(m)  # hi-band channels included; fitter sums them
    return counts


def one(args):
    seed, si_counts = args
    rng = np.random.default_rng(seed)
    ratios = {e: rng.uniform(*RANGES[e]) for e in xm.FIT_ELEMENTS}
    solarang = rng.uniform(10, 80)
    s1 = np.sin(np.radians(90 - solarang))
    G, F = _G["gen"], _G["fit"]
    th = truth_theta(G, ratios, s1, si_counts)
    counts = generate_counts(G, th, s1, rng=rng)
    t0 = time.perf_counter()
    o = F.fit(counts, 8.0, 16, solarang)
    o.pop("corr", None)
    o["fit_ms"] = (time.perf_counter() - t0) * 1000
    o.update(seed=seed, si_counts=si_counts, solarang=solarang)
    for e in xm.FIT_ELEMENTS:
        o[f"true_{e}"] = ratios[e]
    return o


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=200, help="spectra per signal level")
    ap.add_argument("--levels", default="100,300,1000,3000,10000")
    ap.add_argument("--kT_gen", type=float, default=1.07)
    ap.add_argument("--kT_fit", type=float, default=1.07)
    ap.add_argument("--workers", type=int, default=12)
    ap.add_argument("--bkg", default=os.path.join(ROOT, "outputs", "asr", "stack_2024_12.npz"))
    ap.add_argument("--inst", default=None, help="instrumental template .npy (enables gamma)")
    ap.add_argument("--tag", default="")
    ap.add_argument("--solar", action="store_true", help="CHIANTI table; kT_gen/kT_fit are T in MK")
    ap.add_argument("--out_dir", default=os.path.join(ROOT, "outputs", "asr", "synth"))
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    levels = [int(x) for x in args.levels.split(",")]
    jobs = [(1000 * li + i, lv) for li, lv in enumerate(levels) for i in range(args.n)]
    with Pool(args.workers, initializer=_init, initargs=(args.kT_gen, args.kT_fit, args.bkg, args.inst, args.solar)) as p:
        res = pd.DataFrame(p.map(one, jobs, chunksize=4))
    tag = f"gen{args.kT_gen}_fit{args.kT_fit}{args.tag}"
    res.to_parquet(os.path.join(args.out_dir, f"recovery_{tag}.parquet"), index=False)

    rows = []
    for (lv), g in res.groupby("si_counts"):
        row = {"si_counts": lv, "n": len(g), "conv_frac": g["success"].mean(),
               "median_cond": g["cond"].median(), "median_gof_z": g["gof_z"].median(),
               "frac_quality0": (g["quality"] == 0).mean(), "median_fit_ms": g["fit_ms"].median()}
        for e in xm.FIT_ELEMENTS:
            rel = g[f"r_{e}"] / g[f"true_{e}"] - 1
            pull = (g[f"r_{e}"] - g[f"true_{e}"]) / g[f"sd_r_{e}"]
            row[f"{e}_bias_med"] = rel.median()
            row[f"{e}_scatter_rmad"] = 1.4826 * np.median(np.abs(rel - rel.median()))
            row[f"{e}_pull_std"] = 1.4826 * np.median(np.abs(pull - pull.median()))
            row[f"{e}_med_relunc"] = (g[f"sd_r_{e}"] / g[f"true_{e}"]).median()
        rows.append(row)
    summ = pd.DataFrame(rows)
    summ.to_csv(os.path.join(args.out_dir, f"summary_{tag}.csv"), index=False)
    pd.set_option("display.width", 250)
    print(summ.round(3).T.to_string())


if __name__ == "__main__":
    main()
