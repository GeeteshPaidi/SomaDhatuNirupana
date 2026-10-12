"""
Cross-check of the GOES-derived solar spectrum against Chandrayaan-2 XSM for one
month (pre-declared test; analysis choices fixed before running).

1. XSM Level-2 1-min spectra (exposure >= 30 s, open filter) are fitted in
   1.5-10 keV (chi-square with statistical + systematic errors, response
   CH2xsmrspwitharea_open20191214v01 from CALDB 20210628) with the same CHIANTI
   table used for CLASS: isothermal (T, EM) and two-temperature (T1, EM1, T2, EM2).
2. Mapping-set CLASS spectra (quality 0, gamma = 0 model) whose mid-time lies in
   an XSM minute are re-fitted with the identical fitter, calibration and
   background, replacing the GOES temperature by the XSM-derived spectrum.
3. Per-spectrum r_XSM / r_GOES for Mg/Si and Al/Si on pairs with both fits at
   quality 0, and T_XSM / T_GOES.

Inputs: data/xsm/ (PRADAN XSM daily zips bundled in one zip, and the XSM
archive bundle with the CALDB), outputs/asr/months_v3/<month>/.
Outputs: outputs/asr/xsm/.
"""
import argparse
import glob
import io
import json
import os
import sys
import time
import zipfile
from multiprocessing import Pool

import numpy as np
import pandas as pd
from astropy.io import fits
from scipy import optimize

sys.path.insert(0, os.path.dirname(__file__))
import fitter as ft  # noqa: E402
import run_month as rm  # noqa: E402
import solar  # noqa: E402

ROOT = rm.ROOT
OUT = os.path.join(ROOT, "outputs", "asr", "xsm")
XSM_DIR = os.path.join(ROOT, "data", "xsm")
RSP_NAME = "caldb/CH2xsmrspwitharea_open20191214v01.rsp"
MET_REF = pd.Timestamp("2017-01-01")  # MJDREF 57754, UTC
FIT_BAND = (1.5, 10.0)
MIN_EXP = 30.0


# -- XSM data ----------------------------------------------------------------
def _caldb():
    for f in glob.glob(os.path.join(XSM_DIR, "*.zip")):
        z = zipfile.ZipFile(f)
        for n in z.namelist():
            if n.endswith(".zip") and "caldb" in n:
                return zipfile.ZipFile(io.BytesIO(z.read(n)))
    raise FileNotFoundError("XSM CALDB not found in data/xsm")


def load_response():
    h = fits.open(io.BytesIO(_caldb().read(RSP_NAME)))
    m = h["MATRIX"].data
    R = np.array([np.asarray(row, float) for row in m["MATRIX"]])  # (n_E, 512), cm2 incl. area
    eb = h["EBOUNDS"].data
    return m["ENERG_LO"].astype(float), m["ENERG_HI"].astype(float), R, \
        0.5 * (eb["E_MIN"] + eb["E_MAX"]).astype(float)


def load_spectra(month):
    """All XSM L2 1-min spectra of a month from the daily zips (read in memory)."""
    y, mo = month.split("/")
    rows = []
    for f in glob.glob(os.path.join(XSM_DIR, "*.zip")):
        z = zipfile.ZipFile(f)
        for n in z.namelist():
            if not (n.startswith(f"ch2_xsm_{y}{mo}") and n.endswith(".zip")):
                continue
            inner = zipfile.ZipFile(io.BytesIO(z.read(n)))
            pha = [x for x in inner.namelist() if x.endswith("_level2.pha")]
            if not pha:
                continue
            h = fits.open(io.BytesIO(inner.read(pha[0])))
            d = h["SPECTRUM"].data
            for i in range(len(d)):
                rows.append({"t0": float(d["TSTART"][i]), "t1": float(d["TSTOP"][i]),
                             "exposure": float(d["EXPOSURE"][i]), "filt": int(d["FILT_STATUS"][i]),
                             "counts": np.asarray(d["COUNTS"][i], float),
                             "stat_err": np.asarray(d["STAT_ERR"][i], float),
                             "sys_err": np.asarray(d["SYS_ERR"][i], float)})
    df = pd.DataFrame(rows)
    df["start"] = MET_REF + pd.to_timedelta(df.t0, "s")
    df["stop"] = MET_REF + pd.to_timedelta(df.t1, "s")
    return df.sort_values("start").reset_index(drop=True)


# -- XSM spectral model --------------------------------------------------------
class XsmModel:
    """Count rate per XSM channel for CHIANTI isothermal components."""

    def __init__(self, table, e_lo, e_hi, R, ch_e):
        self.tab = table
        # table bins (5 eV) -> response bins (10 eV) inside the table range
        te = table.e_edges
        ok = (e_lo >= te[0] - 1e-6) & (e_hi <= te[-1] + 1e-6)
        self.R = R[ok]
        lo, hi = e_lo[ok], e_hi[ok]
        i0 = np.round((lo - te[0]) / 0.005).astype(int)
        self.i0, self.dE = i0, hi - lo
        self.sel = (ch_e >= FIT_BAND[0]) & (ch_e <= FIT_BAND[1])
        # folded rate per unit EM (1e48) at each table temperature
        self.F = np.array([self._fold(table.S[k]) for k in range(len(table.t_mk))])
        self.logF = np.log(np.maximum(self.F, 1e-300))

    def _rebin(self, s):
        return 0.5 * (s[self.i0] + s[self.i0 + 1])

    def _fold(self, s):
        return (self._rebin(s) * self.dE) @ self.R

    def rate1(self, T_mk):
        lt = np.clip(np.log(T_mk), self.tab.lt[0], self.tab.lt[-1])
        i = int(np.clip(np.searchsorted(self.tab.lt, lt) - 1, 0, len(self.tab.lt) - 2))
        w = (lt - self.tab.lt[i]) / (self.tab.lt[i + 1] - self.tab.lt[i])
        return np.exp((1 - w) * self.logF[i] + w * self.logF[i + 1])

    def fit(self, sp):
        c = sp["counts"][self.sel]
        err = np.sqrt(sp["stat_err"][self.sel] ** 2 + (sp["sys_err"][self.sel] * c) ** 2)
        err = np.maximum(err, 1.0)
        t = sp["exposure"]

        def r1(p):
            return (t * 10 ** p[1] * self.rate1(10 ** p[0]) [self.sel] - c) / err

        best = None
        for lt0 in (0.5, 0.8, 1.1):
            g = self.rate1(10 ** lt0)[self.sel] * t
            p0 = [lt0, np.log10(max(c.sum(), 1.0) / max(g.sum(), 1e-30))]
            r = optimize.least_squares(r1, p0, bounds=([np.log10(2.0), -6], [np.log10(40.0), 6]))
            if best is None or r.cost < best.cost:
                best = r
        out = {"T1_mk": 10 ** best.x[0], "EM1_1e48": 10 ** best.x[1],
               "chi2_1t": 2 * best.cost, "dof_1t": int(self.sel.sum() - 2)}

        def r2(p):
            m = t * (10 ** p[1] * self.rate1(10 ** p[0]) + 10 ** p[3] * self.rate1(10 ** p[2]))
            return (m[self.sel] - c) / err

        best2 = None
        for dT, f in ((0.6, 0.5), (0.8, 0.2), (1.4, 0.05)):
            p0 = [np.log10(max(out["T1_mk"] * dT, 2.01)), np.log10(out["EM1_1e48"] * (1 - f) / dT),
                  np.log10(min(out["T1_mk"] / dT * 1.3, 39.0)), np.log10(out["EM1_1e48"] * f)]
            try:
                r = optimize.least_squares(r2, p0, bounds=([np.log10(2.0), -6, np.log10(2.0), -6],
                                                           [np.log10(40.0), 6, np.log10(40.0), 6]))
            except ValueError:
                continue
            if best2 is None or r.cost < best2.cost:
                best2 = r
        p = best2.x
        (Ta, EMa), (Tb, EMb) = sorted([(10 ** p[0], 10 ** p[1]), (10 ** p[2], 10 ** p[3])])
        out.update(Tc_mk=Ta, EMc_1e48=EMa, Th_mk=Tb, EMh_1e48=EMb, chi2_2t=2 * best2.cost,
                   dof_2t=int(self.sel.sum() - 4))
        return out


# -- two-temperature incident spectrum for the CLASS fitter --------------------
class TwoTempTable(solar.SolarTable):
    """SolarTable whose 'temperature' is a tuple (T1, EM1, T2, EM2)."""

    def spectrum(self, key):
        if np.ndim(key) == 0:
            return super().spectrum(key)
        T1, EM1, T2, EM2 = key
        return EM1 * super().spectrum(T1) + EM2 * super().spectrum(T2)


_F = None


def _init(calib, bkg, two_t):
    global _F
    _F = rm.make_fitter(calib, bkg, None, gamma_zero=True)
    if two_t:
        tab = TwoTempTable()
        _F.fm.solar = tab
        _F.fm.T_mk = None


def _fit(r):
    o = _F.fit(r["counts"], r["exposure"], r["n_scd"], r["solarang"], T_mk=r["T_key"])
    o.pop("corr", None)
    o["source_file"] = r["source_file"]
    return {k: (str(v) if k == "T_mk" else v) for k, v in o.items()}


def _xsm_fit(sp):
    return _XM.fit(sp)


_XM = None


def _xinit(model):
    global _XM
    _XM = model


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--month", default="2024/12")
    ap.add_argument("--months_root", default=os.path.join(ROOT, "outputs", "asr", "months_v3"))
    ap.add_argument("--workers", type=int, default=12)
    ap.add_argument("--limit", type=int, default=None, help="test on a subset of CLASS spectra")
    a = ap.parse_args()
    tag = a.month.replace("/", "_")
    os.makedirs(OUT, exist_ok=True)
    t0 = time.perf_counter()

    # 1. XSM spectra and fits
    xs = load_spectra(a.month)
    xs = xs[(xs.exposure >= MIN_EXP) & (xs.filt == 0)].reset_index(drop=True)
    e_lo, e_hi, R, ch_e = load_response()
    model = XsmModel(solar.SolarTable(), e_lo, e_hi, R, ch_e)

    # CLASS mapping-set spectra inside XSM minutes
    mdir = os.path.join(a.months_root, tag)
    cf = pd.read_parquet(os.path.join(mdir, "fits.parquet"))
    cf = cf[cf.quality == 0].copy()
    tt = pd.to_datetime(cf.mid_utc).values
    k = np.searchsorted(xs.start.values, tt, side="right") - 1
    inside = (k >= 0) & (tt < xs.stop.values[np.clip(k, 0, None)])
    cf = cf[inside].copy()
    cf["xsm_idx"] = k[inside]
    print(f"XSM minutes: {len(xs)}; CLASS quality-0 spectra inside XSM minutes: {len(cf)}", flush=True)
    if a.limit:
        cf = cf.sample(n=min(a.limit, len(cf)), random_state=0)
    need = sorted(cf.xsm_idx.unique())

    with Pool(a.workers, initializer=_xinit, initargs=(model,)) as p:
        xr = p.map(_xsm_fit, [xs.loc[i, ["counts", "stat_err", "sys_err", "exposure"]].to_dict() for i in need],
                   chunksize=8)
    xf = pd.DataFrame(xr, index=need)
    xf["start"] = xs.loc[need, "start"].values
    xf.to_parquet(os.path.join(OUT, f"xsm_fits_{tag}.parquet"))
    print(f"XSM fits: {len(xf)} minutes ({time.perf_counter()-t0:.0f}s); "
          f"median T1 {xf.T1_mk.median():.2f} MK, median chi2/dof 1T {np.median(xf.chi2_1t/xf.dof_1t):.2f}, "
          f"2T {np.median(xf.chi2_2t/xf.dof_2t):.2f}", flush=True)

    # 2. CLASS re-fits
    calib = json.load(open(os.path.join(mdir, "calib.json")))
    bkg = np.load(os.path.join(mdir, "background.npy"))
    day = os.path.join(ROOT, "data", "raw", *a.month.split("/"))
    rows = []
    for _, c in cf.iterrows():
        d = c.source_file.split("_")[3][6:8]
        r = rm.read_one(os.path.join(day, d, c.source_file))
        if r is None:
            continue
        x = xf.loc[c.xsm_idx]
        r.update(goes_T=c.goes_T_mk, T1=x.T1_mk, two=(x.Tc_mk, x.EMc_1e48, x.Th_mk, x.EMh_1e48))
        rows.append(r)
    print(f"read {len(rows)} CLASS spectra", flush=True)
    res = {}
    for name, key, two in [("goes", "goes_T", False), ("xsm1t", "T1", False), ("xsm2t", "two", True)]:
        for r in rows:
            r["T_key"] = r[key]
        t1 = time.perf_counter()
        with Pool(a.workers, initializer=_init, initargs=(calib, bkg, two)) as p:
            out = list(p.imap(_fit, rows, chunksize=16))
        res[name] = pd.DataFrame(out).set_index("source_file")
        print(f"  {name}: {len(out)} fits in {time.perf_counter()-t1:.0f}s; "
              f"quality 0: {(res[name].quality == 0).mean():.3f}", flush=True)

    # 3. comparison
    m = cf.set_index("source_file")[["r_Mg", "r_Al", "goes_T_mk", "xsm_idx"]]
    for n, d in res.items():
        m = m.join(d[["r_Mg", "r_Al", "quality"]].add_suffix(f"_{n}"))
    m["T1_xsm"] = xf.T1_mk.reindex(m.xsm_idx).values
    m.to_parquet(os.path.join(OUT, f"class_refits_{tag}.parquet"))
    summ = {"month": tag, "n_xsm_minutes_month": int(len(xs)), "n_class_inside": int(len(cf)),
            "n_xsm_minutes_used": int(len(xf)),
            "xsm_chi2dof_1t_median": float(np.median(xf.chi2_1t / xf.dof_1t)),
            "xsm_chi2dof_2t_median": float(np.median(xf.chi2_2t / xf.dof_2t)),
            "T_xsm1t_over_goes": list(np.nanpercentile(m.T1_xsm / m.goes_T_mk, [16, 50, 84]).round(4)),
            "T_xsm1t_median_mk": float(np.median(m.T1_xsm)), "T_goes_median_mk": float(np.median(m.goes_T_mk))}
    base = m.quality_goes == 0
    summ["goes_refit_reproduces"] = {e: float(np.nanmax(np.abs(m[f"r_{e}_goes"][base] / m[f"r_{e}"][base] - 1)))
                                     for e in ("Mg", "Al")}
    for n in ("xsm1t", "xsm2t"):
        ok = base & (m[f"quality_{n}"] == 0)
        summ[f"{n}_n_pairs"] = int(ok.sum())
        summ[f"{n}_frac_quality0"] = float((m[f"quality_{n}"] == 0).mean())
        for e in ("Mg", "Al"):
            q = (m[f"r_{e}_{n}"][ok] / m[f"r_{e}_goes"][ok]).values
            summ[f"{n}_{e}_ratio_p16_p50_p84"] = list(np.percentile(q, [16, 50, 84]).round(4))
    summ["wall_s"] = time.perf_counter() - t0
    with open(os.path.join(OUT, f"summary_{tag}{'_test' if a.limit else ''}.json"), "w") as f:
        json.dump(summ, f, indent=1)
    print(json.dumps(summ, indent=1))


if __name__ == "__main__":
    main()
