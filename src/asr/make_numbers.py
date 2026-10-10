"""
Task 5.9: dump every headline number used in the manuscript to
outputs/asr/headline_numbers.json and to paper_asr/numbers.tex (LaTeX
macros), plus the generated tables in paper_asr/tables/. The manuscript
uses only these macros, so text and outputs cannot drift apart.
"""
import ast
import glob
import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
A = os.path.join(ROOT, "outputs", "asr")
P = os.path.join(ROOT, "paper_asr")
N = {}


def put(name, value, fmt=None):
    if value is None or (isinstance(value, float) and not np.isfinite(value)):
        N[name] = "??"
    elif fmt:
        N[name] = format(value, fmt)
    else:
        N[name] = str(value)


def comma(n):
    return f"{int(n):,}".replace(",", "{,}")


def ci(s):
    v = ast.literal_eval(s) if isinstance(s, str) else s
    return f"[{v[0]:.2f}, {v[1]:.2f}]".replace("-", "$-$")


def header_numbers():
    h = pd.read_csv(os.path.join(A, "header_survey.csv"))
    put("NHeaderSample", comma(len(h)))
    put("PctPartialExp", 100 * np.mean(h.EXPOSURE < 8.0), ".1f")
    put("PctFullExp", 100 * np.mean(h.EXPOSURE >= 8.0), ".1f")
    put("PctAllSCD", 100 * np.mean(h.n_scd == 16), ".0f")
    put("PctNotAllSCD", 100 * np.mean(h.n_scd < 16), ".0f")


def stage_numbers():
    s = json.load(open(os.path.join(A, "maps", "summary_2deg.json")))
    st = s["stages"]
    put("NRawFiles", comma(st["raw_L1_files"]))
    put("NDayside", comma(st["dayside_valid_geometry"]))
    put("NGoesOK", comma(st["goes_T_reliable_B1"]))
    put("NConverged", comma(st["fit_converged"]))
    put("NMapped", comma(st["mapping_set"]))
    put("PctMapped", 100 * st["mapping_set"] / st["dayside_valid_geometry"], ".1f")
    months = [d for d in glob.glob(os.path.join(A, "months", "*")) if len(os.path.basename(d)) == 7]
    put("NMonths", len(months))
    fp = s["footprint_size"]
    put("FootprintA", fp["side01_km_median"], ".0f")
    put("FootprintB", fp["side12_km_median"], ".0f")
    c = s["coverage"]
    put("CovAllArea", c["coverage_all_area_pct"], ".1f")
    put("CovUsedArea", c["coverage_used_area_pct"], ".1f")
    put("CovUsedCells", c["coverage_used_cells_pct"], ".2f")
    put("CovFiveArea", c["Mg_cells_with_ge5_footprints_area_pct"], ".1f")
    pr = s["mapping_set_ratio_percentiles"]
    for e in ["Mg", "Al"]:
        put(f"Med{e}", pr[e][2], ".2f")
        put(f"Pfive{e}", pr[e][0], ".2f")
        put(f"Pninetyfive{e}", pr[e][4], ".2f")
    # per-month table
    rows = []
    ms = pd.read_parquet(os.path.join(A, "maps", "mapping_set.parquet"), columns=["month"])
    nmap = ms["month"].value_counts()
    for d in sorted(months):
        cj = json.load(open(os.path.join(d, "calib.json")))
        tag = os.path.basename(d)
        rows.append({"month": tag.replace("_", "-"), "files": cj["n_files"], "dayside": cj["n_dayside"],
                     "mapped": int(nmap.get(tag, 0)), "gain": 1000 * cj["gain"], "offset": 1000 * cj["offset"],
                     "fwhm": cj["fwhm_1.49keV_eV"], "bkg": "own" if cj.get("bkg_template_source", tag) == tag else "Dec 2024"})
    t = pd.DataFrame(rows)
    os.makedirs(os.path.join(P, "tables"), exist_ok=True)
    with open(os.path.join(P, "tables", "months.tex"), "w") as f:
        f.write("\\begin{tabular}{lrrrrrrl}\n\\toprule\nMonth & Files & Dayside & Mapped & $g$ (eV/ch) & $o$ (eV) & FWHM$_{1.49}$ (eV) & Umbra template\\\\\n\\midrule\n")
        for r in rows:
            f.write(f"{r['month']} & {comma(r['files'])} & {comma(r['dayside'])} & {comma(r['mapped'])} & "
                    f"{r['gain']:.3f} & {r['offset']:.1f} & {r['fwhm']:.0f} & {r['bkg']}\\\\\n")
        f.write("\\bottomrule\n\\end{tabular}\n")


def calib_numbers():
    cj = json.load(open(os.path.join(A, "months", "2024_12", "calib.json")))
    put("NStackDec", comma(cj["n_stack"]))
    put("GainDec", 1000 * cj["gain"], ".2f")
    put("OffsetDec", 1000 * cj["offset"], ".1f")
    put("SigmaZeroDec", 1000 * cj["sigma0"], ".0f")
    put("FwhmAlDec", cj["fwhm_1.49keV_eV"], ".0f")


def background_numbers():
    d = pd.read_parquet(os.path.join(A, "bkg_check3_2024_12.parquet"))
    n = d[d.sa >= 90]
    put("AlSunlit", n.al[n.sc_sunlit].median(), ".1f")
    put("AlUmbra", n.al[~n.sc_sunlit].median(), ".2f")
    put("NSunlitDec", comma(int(n.sc_sunlit.sum())))
    put("NUmbraDec", comma(int((~n.sc_sunlit).sum())))
    s = n[n.sc_sunlit]
    put("RhoAlSoft", s.al.corr(s.lo, method="spearman"), ".2f")
    put("RhoAlHard", s.al.corr(s.hi, method="spearman"), ".2f")
    g = pd.read_csv(os.path.join(A, "inst_vs_goes_2024_12.csv"))
    put("RhoAlGoes", g.al.corr(g.phi, method="spearman"), ".2f")


def solar_numbers():
    import solar
    import xrf_model as xm
    st = solar.SolarTable()
    e = st.e_mid
    import astropy.units as u  # noqa: F401
    from sunkit_spex.legacy import thermal
    edges = st.e_edges * u.keV
    fr = {"Mg": [], "Si": []}
    for T in [3, 4, 6, 10, 15]:
        f = np.asarray(thermal.thermal_emission(edges, T * u.MK, 1e48 / u.cm ** 3).value)
        c = np.asarray(thermal.continuum_emission(edges, T * u.MK, 1e48 / u.cm ** 3).value)
        for el, ed in [("Mg", 1.305), ("Si", 1.839)]:
            m = e > ed
            fr[el].append(100 * (1 - c[m].sum() / f[m].sum()))
    put("LineFracMgLo", min(fr["Mg"]), ".0f")
    put("LineFracMgHi", max(fr["Mg"]), ".0f")
    put("LineFracSiLo", min(fr["Si"]), ".0f")
    put("LineFracSiHi", max(fr["Si"]), ".0f")
    # secondary/primary with the CHIANTI spectrum
    mx = 0.0
    for T in [3, 5, 8, 15]:
        fm = xm.ForwardModel(secondary=True, e_max=15, de=0.005, solar=st, T_mk=T)
        for r in [{"Mg": 0.25, "Al": 0.30, "Ca": 0.40, "Ti": 0.05, "Fe": 0.60},
                  {"Mg": 0.20, "Al": 0.70, "Ca": 0.55, "Ti": 0.01, "Fe": 0.15}]:
            p, s, _ = fm.line_intensities(r, 45.0, return_parts=True)
            ka = [i for i, ln in enumerate(fm.lines) if ln.label == "Ka"]
            mx = max(mx, np.max(s[ka] / p[ka]))
    put("SecMaxPct", 100 * mx, ".1f")
    t = pd.read_parquet(solar.TEM)
    t = t[t.xrsb >= 1e-7]
    for y in [2024, 2025, 2026]:
        put(f"GoesTmed{['Zero','One','Two','Three','Four','Five','Six'][y - 2020]}", t[t.index.year == y].temperature.median(), ".1f")
    ms = pd.read_parquet(os.path.join(A, "maps", "mapping_set.parquet"), columns=["goes_T_mk"])
    put("MapTmed", ms.goes_T_mk.median(), ".1f")
    put("MapTlo", ms.goes_T_mk.quantile(0.05), ".1f")
    put("MapThi", ms.goes_T_mk.quantile(0.95), ".1f")


def synthetic_numbers():
    f = os.path.join(A, "synth", "summary_gen5.0_fit5.0_final.csv")
    if not os.path.exists(f):
        return
    s = pd.read_csv(f)
    put("SynMaxBiasMgAl", 100 * max(s.Mg_bias_med.abs().max(), s.Al_bias_med.abs().max()), ".0f")
    pulls = pd.concat([s.Mg_pull_std, s.Al_pull_std])
    put("SynPullLo", pulls.min(), ".2f")
    put("SynPullHi", pulls.max(), ".2f")
    with open(os.path.join(P, "tables", "synthetic.tex"), "w") as fo:
        fo.write("\\begin{tabular}{lrrrr}\n\\toprule\nSi K$\\alpha$ counts & " +
                 " & ".join(f"{int(x):,}" for x in s.si_counts) + "\\\\\n\\midrule\n")
        for e in ["Mg", "Al", "Ca"]:
            fo.write(f"{e}/Si bias / scatter (\\%) & " + " & ".join(
                f"{100*b:+.0f} / {100*c:.0f}" for b, c in zip(s[f'{e}_bias_med'], s[f'{e}_scatter_rmad'])) + "\\\\\n")
            fo.write(f"{e}/Si pull width & " + " & ".join(f"{x:.2f}" for x in s[f'{e}_pull_std']) + "\\\\\n")
        fo.write("Quality-0 fraction & " + " & ".join(f"{x:.2f}" for x in s.frac_quality0) + "\\\\\n")
        fo.write("\\bottomrule\n\\end{tabular}\n")
    for tag, nm in [("gen6.5_fit5.0_final", "Hot"), ("gen4.0_fit5.0_final", "Cold")]:
        f2 = os.path.join(A, "synth", f"summary_{tag}.csv")
        if os.path.exists(f2):
            s2 = pd.read_csv(f2)
            for e in ["Mg", "Al"]:
                lo, hi = 100 * s2[f"{e}_bias_med"].min(), 100 * s2[f"{e}_bias_med"].max()
                put(f"Syn{nm}{e}", f"{lo:+.0f} to {hi:+.0f}".replace("-", "$-$"))


def validation_numbers():
    V = os.path.join(A, "validation")
    s = pd.read_csv(os.path.join(V, "sites_stats.csv"))
    g = pd.read_csv(os.path.join(V, "grs_stats.csv"))
    names = {1.0: "One", 2.0: "Two", 5.0: "Five"}
    for _, r in s.iterrows():
        k = f"Site{r.element}{names[r.radius_deg]}"
        put(k + "R", r.pearson_r, ".2f")
        N[k + "CI"] = ci(r.pearson_ci95)
        put(k + "Rho", r.spearman_rho, ".2f")
        put(k + "RhoP", r.spearman_p, ".2f")
        put(k + "LooLo", r.loo_min, ".2f")
        put(k + "LooHi", r.loo_max, ".2f")
        put(k + "Rmse", r.rmse, ".3f")
        put(k + "Bias", r.bias, "+.3f")
    for _, r in g.iterrows():
        k = f"Grs{r.element}"
        put(k + "R", r.pearson_r, ".2f")
        N[k + "CI"] = ci(r.pearson_ci95)
        put(k + "Rho", r.spearman_rho, ".2f")
        put(k + "N", comma(r.n))
    v1 = pd.read_csv(os.path.join(V, "grs_stats_v1_prereg.csv"))
    r1 = v1[v1.element == "Al"].iloc[0]
    put("GrsAlRvOne", r1.pearson_r, ".2f")
    N["GrsAlCIvOne"] = ci(r1.pearson_ci95)
    for tag, nm in [("_min1000", "MinK"), ("_min3000", "MinKKK")]:
        gg = pd.read_csv(os.path.join(V, f"grs_stats{tag}.csv"))
        for e in ["Mg", "Al"]:
            put(f"Grs{e}R{nm}", gg[gg.element == e].iloc[0].pearson_r, ".2f")
        sm = json.load(open(os.path.join(A, "maps", f"summary_2deg{tag}.json")))
        put(f"NMapped{nm}", comma(sm["stages"]["mapping_set"]))
    # site table
    t = pd.read_csv(os.path.join(V, "sites_table.csv"))
    t2 = t[t.radius_deg == 2.0].pivot_table(index="site", columns="element",
                                            values=["ref_ratio", "class_ratio", "class_se", "n_footprints"], aggfunc="first")
    order = ["11", "12", "14", "15", "16", "17", "L16", "L20", "L24"]
    label = {s: (f"Apollo {s}" if not s.startswith("L") else f"Luna {s[1:]}") for s in order}
    with open(os.path.join(P, "tables", "sites.tex"), "w") as f:
        f.write("\\begin{tabular}{lrrrrrr}\n\\toprule\nSite & Lat & Lon & \\multicolumn{2}{c}{Mg/Si} & \\multicolumn{2}{c}{Al/Si}\\\\\n"
                " & (\\textdegree) & (\\textdegree E) & soil & CLASS & soil & CLASS\\\\\n\\midrule\n")
        for s_ in order:
            if s_ not in t2.index:
                continue
            row = t[(t.site == s_)].iloc[0]
            f.write(f"{label[s_]} & {row.lat:.2f} & {row.lon:.2f} & {t2.loc[s_, ('ref_ratio','Mg')]:.3f} & "
                    f"{t2.loc[s_, ('class_ratio','Mg')]:.3f}$\\pm${t2.loc[s_, ('class_se','Mg')]:.3f} & "
                    f"{t2.loc[s_, ('ref_ratio','Al')]:.3f} & {t2.loc[s_, ('class_ratio','Al')]:.3f}$\\pm${t2.loc[s_, ('class_se','Al')]:.3f}\\\\\n")
        f.write("\\bottomrule\n\\end{tabular}\n")
    nfp = t[(t.radius_deg == 2.0)].n_footprints
    put("SiteFpMin", comma(nfp.min()))
    put("SiteFpMax", comma(nfp.max()))


def tsens_numbers():
    base_rows, rows = [], []
    for sc, nm in [(0.8, "Lo"), (1.2, "Hi")]:
        diffs = {"Mg": [], "Al": []}
        for m in ["2024_12", "2025_06", "2025_12"]:
            f = os.path.join(A, "months", f"{m}_Tx{sc}", "fits.parquet")
            if not os.path.exists(f):
                continue
            a = pd.read_parquet(f)
            b = pd.read_parquet(os.path.join(A, "months", m, "fits.parquet"))
            j = a.merge(b, on="source_file", suffixes=("_s", "_b"))
            ok = (j.quality_s == 0) & (j.quality_b == 0) & (j.goes_xrsb_b >= 1e-7)
            for e in ["Mg", "Al"]:
                ok2 = ok & (j[f"r_{e}_b"] > 1e-3) & (j[f"r_{e}_s"] > 1e-3)
                diffs[e].append(np.log(j.loc[ok2, f"r_{e}_s"] / j.loc[ok2, f"r_{e}_b"]))
        for e in ["Mg", "Al"]:
            if diffs[e]:
                d = pd.concat(diffs[e])
                put(f"Tsens{e}{nm}", 100 * (np.exp(d.median()) - 1), "+.0f")
                N[f"Tsens{e}{nm}"] = N[f"Tsens{e}{nm}"].replace("-", "$-$")


def main():
    for fn in [header_numbers, stage_numbers, calib_numbers, background_numbers, solar_numbers,
               synthetic_numbers, validation_numbers, tsens_numbers]:
        try:
            fn()
        except Exception as ex:  # noqa: BLE001
            print(f"WARNING {fn.__name__}: {ex}")
    with open(os.path.join(A, "headline_numbers.json"), "w") as f:
        json.dump(N, f, indent=1, sort_keys=True)
    with open(os.path.join(P, "numbers.tex"), "w") as f:
        f.write("% AUTO-GENERATED by src/asr/make_numbers.py -- do not edit\n")
        for k, v in sorted(N.items()):
            f.write(f"\\newcommand{{\\{k}}}{{{v}}}\n")
    print(json.dumps(N, indent=1, sort_keys=True))


if __name__ == "__main__":
    main()
