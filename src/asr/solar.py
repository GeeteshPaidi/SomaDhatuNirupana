"""
Time-resolved solar spectrum for the ASR rework (TRACKER Q3 decision: GOES).

1. GOES XRS 1-min science data (NOAA NCEI, GOES-16 for 2020-2022, GOES-18
   for 2024-2026) -> isothermal temperature T and emission measure EM per
   minute with sunkit_instruments.goes_xrs.calculate_temperature_em
   (White et al. 2005 method, CHIANTI response, coronal abundances).
2. Isothermal CHIANTI photon spectrum (lines + continuum, coronal
   abundances) from sunkit_spex.legacy.thermal.thermal_emission, tabulated
   on a temperature grid; the solar spectrum for an exposure is the table
   interpolated (in log T) at the GOES temperature of its mid-time.

Outputs: outputs/asr/solar/goes_tem.parquet, outputs/asr/solar/spectrum_table.npz
"""
import argparse
import glob
import os

import numpy as np
import pandas as pd

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
OUT = os.path.join(ROOT, "outputs", "asr", "solar")
TABLE = os.path.join(OUT, "spectrum_table.npz")
TEM = os.path.join(OUT, "goes_tem.parquet")

# Excitation / scatter grid (keV bin edges). Lower bound respects the
# sunkit_spex table limit (1.0003 keV).
E_EDGES = np.round(np.arange(1.01, 20.0 + 1e-9, 0.005), 6)
T_GRID_MK = np.geomspace(2.0, 40.0, 61)

# CLASS months -> GOES satellite/year files
GOES_FILES = {2020: "g16", 2021: "g16", 2022: "g16", 2024: "g18", 2025: "g18", 2026: "g18"}


def build_table():
    import astropy.units as u
    from sunkit_spex.legacy import thermal
    e = E_EDGES * u.keV
    rows = []
    for T in T_GRID_MK:
        f = thermal.thermal_emission(e, T * u.MK, 1e48 / u.cm ** 3)  # ph / (s cm2 keV) at 1 AU
        rows.append(np.asarray(f.value, float))
    S = np.array(rows)
    os.makedirs(OUT, exist_ok=True)
    np.savez(TABLE, e_edges=E_EDGES, t_mk=T_GRID_MK, flux_em48=S)
    return S


class SolarTable:
    """S(E; T) for EM = 1e48 cm^-3, photons / (s cm^2 keV) at 1 AU."""

    def __init__(self, path=TABLE):
        d = np.load(path)
        self.e_edges, self.t_mk, self.S = d["e_edges"], d["t_mk"], d["flux_em48"]
        self.e_mid = 0.5 * (self.e_edges[1:] + self.e_edges[:-1])
        self.lt = np.log(self.t_mk)

    def spectrum(self, T_mk):
        """Bin-averaged photon flux on e_mid for temperature T (log-T interpolation,
        clipped to the table range)."""
        lt = np.clip(np.log(T_mk), self.lt[0], self.lt[-1])
        i = int(np.clip(np.searchsorted(self.lt, lt) - 1, 0, len(self.lt) - 2))
        w = (lt - self.lt[i]) / (self.lt[i + 1] - self.lt[i])
        return (1 - w) * self.S[i] + w * self.S[i + 1]

    def on_grid(self, T_mk, E):
        """Spectrum resampled (bin-average preserved) onto grid centres E with
        the same 5 eV spacing; E must lie within the table."""
        s = self.spectrum(T_mk)
        idx = np.clip(np.round((E - self.e_mid[0]) / 0.005).astype(int), 0, len(s) - 1)
        return s[idx]


def build_goes_tem():
    import sunpy.timeseries as ts
    from sunkit_instruments import goes_xrs
    frames = []
    for year, sat in GOES_FILES.items():
        f = glob.glob(os.path.join(ROOT, "data", "goes", f"sci_xrsf-l2-avg1m_{sat}_y{year}_*.nc"))
        if not f:
            print(f"missing GOES file for {year}")
            continue
        g = ts.TimeSeries(f[0])
        df = g.to_dataframe()
        good = (df["xrsa_quality"] == 0) & (df["xrsb_quality"] == 0) & (df["xrsa"] > 0) & (df["xrsb"] > 0)
        from sunpy.timeseries.sources.goes import XRSTimeSeries
        g2 = XRSTimeSeries(df[good][["xrsa", "xrsb"]], g.meta, g.units)
        tem = goes_xrs.calculate_temperature_em(g2).to_dataframe()
        out = df[good][["xrsa", "xrsb"]].join(tem)
        out["sat"] = sat
        frames.append(out)
        print(year, sat, len(out), "minutes; T median", float(np.nanmedian(out["temperature"])))
    tem = pd.concat(frames)
    tem.index.name = "time"
    tem.to_parquet(TEM)
    return tem


def goes_at(times, tem=None, window_min=1):
    """Nearest-minute GOES T (MK) and EM (cm^-3) for an array of UTC times;
    NaN if no good GOES minute within `window_min`."""
    if tem is None:
        tem = pd.read_parquet(TEM)
    t = pd.to_datetime(pd.Series(times)).dt.tz_localize(None)
    idx = tem.index.get_indexer(t, method="nearest", tolerance=pd.Timedelta(minutes=window_min))
    T = np.where(idx >= 0, tem["temperature"].values[idx], np.nan)
    EM = np.where(idx >= 0, tem["emission_measure"].values[idx], np.nan)
    B = np.where(idx >= 0, tem["xrsb"].values[idx], np.nan)
    return T, EM, B


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--table", action="store_true")
    ap.add_argument("--goes", action="store_true")
    a = ap.parse_args()
    if a.table:
        S = build_table()
        print("table", S.shape)
    if a.goes:
        build_goes_tem()
