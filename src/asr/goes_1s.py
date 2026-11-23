"""
R1-3 sensitivity: exposure-matched GOES temperatures from the 1 s XRS data.

For each dayside CLASS exposure of a month (taken from an existing
fits.parquet), the GOES-18 1 s XRS-A and XRS-B fluxes are averaged over the
exposure window [mid - t/2, mid + t/2] (good-quality samples only) and
converted to T and EM with the same White et al. (2005) / CHIANTI method used
for the 1 min data. The result is written as a TEM table indexed by the
exposure mid-time, which run_month.py --goes_tem can use in place of
goes_tem.parquet.

Usage: python src/asr/goes_1s.py --month 2024_12 --fits outputs/asr/months_v4/2024_12/fits.parquet
"""
import argparse
import glob
import os

import numpy as np
import pandas as pd

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
OUT = os.path.join(ROOT, "outputs", "asr", "solar")


def load_1s(month):
    import sunpy.timeseries as ts
    y, m = month.split("_")
    files = sorted(glob.glob(os.path.join(ROOT, "data", "goes", "flx1s", f"sci_xrsf-l2-flx1s_g18_d{y}{m}??_*.nc")))
    if not files:
        raise SystemExit(f"no 1 s GOES files for {month}")
    frames = []
    for f in files:
        df = ts.TimeSeries(f).to_dataframe()
        good = (df["xrsa_quality"] == 0) & (df["xrsb_quality"] == 0) & (df["xrsa"] > 0) & (df["xrsb"] > 0)
        frames.append(df.loc[good, ["xrsa", "xrsb"]])
    g = pd.concat(frames).sort_index()
    g.index = g.index.tz_localize(None) if g.index.tz is not None else g.index
    return g, ts.TimeSeries(files[0])


def rescale_to_avg1m(g, month):
    """Scale each 1 s sample so that its minute mean equals the avg1m flux of
    that minute (per channel). The flx1s and avg1m products differ in XRS-A
    calibration (flx1s ~14% higher in the median for 2024-12, up to ~2.6x in
    quiet periods; XRS-B within 1%), so without this the test would mix a
    product-calibration offset into the time-resolution effect. Minutes without
    a good avg1m value are dropped."""
    import sunpy.timeseries as ts
    y = month.split("_")[0]
    f = glob.glob(os.path.join(ROOT, "data", "goes", f"sci_xrsf-l2-avg1m_g18_y{y}_*.nc"))[0]
    a = ts.TimeSeries(f).to_dataframe()
    a.index = a.index.tz_localize(None) if a.index.tz is not None else a.index
    a = a[(a.xrsa_quality == 0) & (a.xrsb_quality == 0) & (a.xrsa > 0) & (a.xrsb > 0)]
    minute = g.index.floor("min")
    mean1s = g.groupby(minute).transform("mean")
    ref = a[["xrsa", "xrsb"]].reindex(minute).set_axis(g.index)
    out = g * (ref / mean1s)
    return out.dropna()


def window_means(g, mid, expo):
    """Mean XRS-A/B over each exposure window via cumulative sums."""
    t = g.index.values.astype("datetime64[ns]").astype(np.int64)
    ca = np.r_[0.0, np.cumsum(g["xrsa"].values)]
    cb = np.r_[0.0, np.cumsum(g["xrsb"].values)]
    half = (np.asarray(expo) * 0.5e9).astype(np.int64)
    m = mid.values.astype("datetime64[ns]").astype(np.int64)
    i0 = np.searchsorted(t, m - half, side="left")
    i1 = np.searchsorted(t, m + half, side="right")
    n = i1 - i0
    with np.errstate(invalid="ignore", divide="ignore"):
        a = (ca[i1] - ca[i0]) / n
        b = (cb[i1] - cb[i0]) / n
    a[n == 0] = np.nan
    b[n == 0] = np.nan
    return a, b, n


def main():
    from sunkit_instruments import goes_xrs
    from sunpy.timeseries.sources.goes import XRSTimeSeries
    ap = argparse.ArgumentParser()
    ap.add_argument("--month", required=True)
    ap.add_argument("--fits", required=True, help="fits.parquet giving the exposures (mid_utc, exposure)")
    ap.add_argument("--raw", action="store_true", help="do not rescale the 1 s fluxes to avg1m")
    a = ap.parse_args()

    ex = pd.read_parquet(a.fits, columns=["mid_utc", "exposure", "goes_T_mk"])
    mid = pd.to_datetime(ex["mid_utc"]).dt.tz_localize(None)
    g, proto = load_1s(a.month)
    if not a.raw:
        g = rescale_to_avg1m(g, a.month)
    xa, xb, n = window_means(g, mid, ex["exposure"].values)
    ok = np.isfinite(xa) & np.isfinite(xb) & (xa > 0) & (xb > 0)
    df = pd.DataFrame({"xrsa": xa[ok], "xrsb": xb[ok]}, index=pd.DatetimeIndex(mid[ok].values, name="time"))
    df = df[~df.index.duplicated()].sort_index()
    tem = goes_xrs.calculate_temperature_em(XRSTimeSeries(df, proto.meta, proto.units)).to_dataframe()
    out = df.join(tem)
    out["n_samples"] = pd.Series(n[ok], index=pd.DatetimeIndex(mid[ok].values)).groupby(level=0).first()
    path = os.path.join(OUT, f"goes_tem_1s{'_raw' if a.raw else ''}_{a.month}.parquet")
    out.to_parquet(path)

    # summary: exposure-matched 1 s T vs nearest-minute 1 min T
    t1m = ex["goes_T_mk"].values[ok]
    t1s = out["temperature"].reindex(pd.DatetimeIndex(mid[ok].values)).values
    both = np.isfinite(t1m) & np.isfinite(t1s)
    r = t1s[both] / t1m[both]
    print(f"{a.month}: {ok.sum()}/{len(ex)} exposures with 1 s GOES; median samples {np.median(n[ok]):.0f}")
    print(f"  T_1s/T_1min: median {np.median(r):.4f}, 5-95% [{np.percentile(r,5):.3f}, {np.percentile(r,95):.3f}], "
          f"|dT|>10%: {np.mean(np.abs(r-1)>0.1):.3%}")
    print("wrote", path)


if __name__ == "__main__":
    main()
