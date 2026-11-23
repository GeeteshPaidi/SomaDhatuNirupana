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
# Primary product = v5 (revision 1): v3 (instrumental amplitude fixed to 0 on the dayside)
# + CLASS ARF + TEMP-binned detector calibration + 8-17 keV proxy band.
PT = "_v5"
PM = "months_v5"
# Archives of the primary mapping set (the 2022-2023 archives form a separate epoch, epoch_numbers()).
PRIMARY_MONTHS = "2020_12:2021_02,2022_03,2024_03,2024_11:2026_05"


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
    # A bound that rounds to -0.00 is shown with three decimals so the sign stays meaningful.
    f = lambda x: f"{x:.3f}" if x != 0 and abs(x) < 0.005 else f"{x:.2f}"
    return f"[{f(v[0])}, {f(v[1])}]".replace("-", "$-$")


def sgn(x):
    """Signed integer percent for LaTeX: 0 has no sign, minus as $-$."""
    v = int(round(x))
    return "0" if v == 0 else (f"+{v}" if v > 0 else f"$-${-v}")


def header_numbers():
    h = pd.read_csv(os.path.join(A, "header_survey.csv"))
    put("NHeaderSample", comma(len(h)))
    put("PctPartialExp", 100 * np.mean(h.EXPOSURE < 8.0), ".1f")
    put("PctFullExp", 100 * np.mean(h.EXPOSURE >= 8.0), ".1f")
    put("PctAllSCD", 100 * np.mean(h.n_scd == 16), ".0f")
    put("PctNotAllSCD", 100 * np.mean(h.n_scd < 16), ".0f")


def stage_numbers():
    s = json.load(open(os.path.join(A, "maps", f"summary_2deg{PT}.json")))
    st = s["stages"]
    put("NRawFiles", comma(st["raw_L1_files"]))
    put("NDayside", comma(st["dayside_valid_geometry"]))
    put("NGoesOK", comma(st["goes_T_reliable_B1"]))
    put("NConverged", comma(st["fit_converged"]))
    put("NMapped", comma(st["mapping_set"]))
    put("PctMapped", 100 * st["mapping_set"] / st["dayside_valid_geometry"], ".1f")
    from maps import parse_months
    keep = parse_months(PRIMARY_MONTHS)
    months = [d for d in glob.glob(os.path.join(A, PM, "*"))
              if len(os.path.basename(d)) == 7 and keep(os.path.basename(d))]
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
    ms = pd.read_parquet(os.path.join(A, "maps", f"mapping_set{PT}.parquet"), columns=["month"])
    nmap = ms["month"].value_counts()
    for d in sorted(months):
        cj = json.load(open(os.path.join(d, "calib.json")))
        tag = os.path.basename(d)
        tc = cj.get("temp_cal") or []
        fw = [x["fwhm_1.49keV_eV"] for x in tc] or [cj["fwhm_1.49keV_eV"]]
        src = cj.get("bkg_template_source", tag)
        rows.append({"month": tag.replace("_", "-"), "files": cj["n_files"], "dayside": cj["n_dayside"],
                     "mapped": int(nmap.get(tag, 0)), "gain": 1000 * cj["gain"], "offset": 1000 * cj["offset"],
                     "fwhm": cj["fwhm_1.49keV_eV"], "fwlo": min(fw), "fwhi": max(fw),
                     "tlo": tc[0]["t_lo"] if tc else np.nan, "thi": tc[-1]["t_hi"] if tc else np.nan,
                     "bkg": "own" if src == tag else pd.Timestamp(src.replace("_", "-") + "-01").strftime("%b %Y")})
    t = pd.DataFrame(rows)
    # per-month quality-pass fraction (mapped / dayside)
    t["frac"] = 100 * t.mapped / t.dayside
    lo, hi = t.loc[t.frac.idxmin()], t.loc[t.frac.idxmax()]
    mname = lambda m: pd.Timestamp(m + "-01").strftime("%B %Y")  # noqa: E731
    put("PassMinPct", lo.frac, ".0f")
    put("PassMinMonth", mname(lo.month))
    put("PassMaxPct", hi.frac, ".0f")
    put("PassMaxMonth", mname(hi.month))
    put("PassMedPct", t.frac.median(), ".0f")
    os.makedirs(os.path.join(P, "tables"), exist_ok=True)
    with open(os.path.join(P, "tables", "months.tex"), "w") as f:
        f.write("\\begin{tabular}{lrrrrrrrl}\n\\toprule\nMonth & Files & Dayside & Mapped & $g$ (eV/ch) & $o$ (eV) & "
                "TEMP (\\textdegree C) & FWHM$_{1.49}$ (eV) & Umbra template\\\\\n\\midrule\n")
        for r in rows:
            tr = "--" if not np.isfinite(r["tlo"]) else f"{r['tlo']:.0f} to {r['thi']:.0f}"
            fr = f"{r['fwlo']:.0f}" if round(r["fwlo"]) == round(r["fwhi"]) else f"{r['fwlo']:.0f}--{r['fwhi']:.0f}"
            f.write(f"{r['month']} & {comma(r['files'])} & {comma(r['dayside'])} & {comma(r['mapped'])} & "
                    f"{r['gain']:.3f} & {r['offset']:.1f} & {tr} & {fr} & {r['bkg']}\\\\\n")
        f.write("\\bottomrule\n\\end{tabular}\n")


def calib_numbers():
    cj = json.load(open(os.path.join(A, PM, "2024_12", "calib.json")))
    put("NStackDec", comma(cj["n_stack"]))
    put("GainDec", 1000 * cj["gain"], ".2f")
    put("GainPctDec", 100 * (1000 * cj["gain"] / 13.5 - 1), ".1f")  # vs header GAIN = 13.5 eV/ch
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
    ms = pd.read_parquet(os.path.join(A, "maps", f"mapping_set{PT}.parquet"), columns=["goes_T_mk"])
    put("MapTmed", ms.goes_T_mk.median(), ".1f")
    put("MapTlo", ms.goes_T_mk.quantile(0.05), ".1f")
    put("MapThi", ms.goes_T_mk.quantile(0.95), ".1f")


def synthetic_numbers():
    f = os.path.join(A, "synth", f"summary_gen5.0_fit5.0{PT}.csv")
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
                f"{sgn(100*b)} / {100*c:.0f}" for b, c in zip(s[f'{e}_bias_med'], s[f'{e}_scatter_rmad'])) + "\\\\\n")
            fo.write(f"{e}/Si pull width & " + " & ".join(f"{x:.2f}" for x in s[f'{e}_pull_std']) + "\\\\\n")
        fo.write("Quality-0 fraction & " + " & ".join(f"{x:.2f}" for x in s.frac_quality0) + "\\\\\n")
        fo.write("\\bottomrule\n\\end{tabular}\n")
    for tag, nm in [(f"gen6.5_fit5.0{PT}", "Hot"), (f"gen4.0_fit5.0{PT}", "Cold")]:
        f2 = os.path.join(A, "synth", f"summary_{tag}.csv")
        if os.path.exists(f2):
            s2 = pd.read_csv(f2)
            for e in ["Mg", "Al"]:
                lo, hi = 100 * s2[f"{e}_bias_med"].min(), 100 * s2[f"{e}_bias_med"].max()
                put(f"Syn{nm}{e}", f"{sgn(lo)} to {sgn(hi)}")


def validation_numbers():
    V = os.path.join(A, "validation")
    s = pd.read_csv(os.path.join(V, f"sites_stats{PT}.csv"))
    g = pd.read_csv(os.path.join(V, f"grs_stats{PT}.csv"))
    names = {1.0: "One", 2.0: "Two", 5.0: "Five"}
    for _, r in s.iterrows():
        k = f"Site{r.element}{names[r.radius_deg]}"
        put(k + "R", r.pearson_r, ".2f")
        N[k + "CI"] = ci(r.pearson_ci95)
        put(k + "Rho", r.spearman_rho, ".2f")
        N[k + "RhoP"] = ("$<$0.001" if r.spearman_p < 0.001 else f"{r.spearman_p:.3f}" if r.spearman_p < 0.01 else f"{r.spearman_p:.2f}")
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
    if os.path.exists(os.path.join(V, "grs_stats_v1_prereg.csv")):
        v1 = pd.read_csv(os.path.join(V, "grs_stats_v1_prereg.csv"))
        r1 = v1[v1.element == "Al"].iloc[0]
        put("GrsAlRvOne", r1.pearson_r, ".2f")
        N["GrsAlCIvOne"] = ci(r1.pearson_ci95)
    for tag, nm in [("_min1000", "MinK"), ("_min3000", "MinKKK")]:
        if not os.path.exists(os.path.join(V, f"grs_stats{PT}{tag}.csv")):
            continue
        gg = pd.read_csv(os.path.join(V, f"grs_stats{PT}{tag}.csv"))
        for e in ["Mg", "Al"]:
            put(f"Grs{e}R{nm}", gg[gg.element == e].iloc[0].pearson_r, ".2f")
        sm = json.load(open(os.path.join(A, "maps", f"summary_2deg{PT}{tag}.json")))
        put(f"NMapped{nm}", comma(sm["stages"]["mapping_set"]))
    # site table
    t = pd.read_csv(os.path.join(V, f"sites_table{PT}.csv"))
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
            f = os.path.join(A, PM, f"{m}_Tx{sc}", "fits.parquet")
            if not os.path.exists(f):
                continue
            a = pd.read_parquet(f)
            b = pd.read_parquet(os.path.join(A, PM, m, "fits.parquet"))
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


def diag_numbers():
    V = os.path.join(A, "validation")
    d = json.load(open(os.path.join(V, f"score_diagnostics{PT}.json")))
    put("GrsSitesMg", d["grs_vs_sites_r"]["Mg"], ".2f")
    put("GrsSitesAl", d["grs_vs_sites_r"]["Al"], ".2f")
    for e in ["Mg", "Al"]:
        put(f"Rel{e}", d[f"class_reliability_{e}"], ".2f")
        put(f"RelGrs{e}", d[f"grs_reliability_{e}"], ".2f")
        put(f"Ceil{e}", d[f"class_grs_ceiling_r_{e}"], ".2f")
        if f"old_pipeline_vs_grs_r_{e}" in d:
            put(f"OldGrs{e}", d[f"old_pipeline_vs_grs_r_{e}"], ".2f")
    m = json.load(open(os.path.join(V, f"alt_median_posthoc{PT}.json")))
    put("MedGrsMg", m["median_Mg_grs_r"], ".2f")
    put("MedGrsAl", m["median_Al_grs_r"], ".2f")
    put("MedRelAl", m["median_Al_reliability"], ".2f")
    # First-stage gamma = 0 model (v3fix): controlled reference for the free-gamma alternative.
    g3p = os.path.join(V, "grs_stats_v3fix.csv")
    if os.path.exists(g3p):
        g3 = pd.read_csv(g3p)
        for e in ["Mg", "Al"]:
            put(f"Grs{e}RvThree", g3[g3.element == e].iloc[0].pearson_r, ".2f")
    # Earlier-model (v2, gamma0) numbers: kept only when their inputs exist on this machine.
    if not os.path.exists(os.path.join(V, "score_diagnostics.json")):
        return
    d2 = json.load(open(os.path.join(V, "score_diagnostics.json")))
    put("RelAlvTwo", d2["class_reliability_Al"], ".2f")
    put("RelMgvTwo", d2["class_reliability_Mg"], ".2f")
    put("CeilAlvTwo", d2["class_grs_ceiling_r_Al"], ".2f")
    # v2 (free gamma) primary results, reported as the alternative model
    sv2 = pd.read_csv(os.path.join(V, "sites_stats.csv"))
    gv2 = pd.read_csv(os.path.join(V, "grs_stats.csv"))
    for e in ["Mg", "Al"]:
        put(f"Grs{e}RvTwo", gv2[gv2.element == e].iloc[0].pearson_r, ".2f")
        N[f"Grs{e}CIvTwo"] = ci(gv2[gv2.element == e].iloc[0].pearson_ci95)
        r = sv2[(sv2.radius_deg == 2.0) & (sv2.element == e)].iloc[0]
        put(f"Site{e}TwoRvTwo", r.pearson_r, ".2f")
        N[f"Site{e}TwoCIvTwo"] = ci(r.pearson_ci95)
    smv2 = json.load(open(os.path.join(A, "maps", "summary_2deg.json")))
    put("NMappedvTwo", comma(smv2["stages"]["mapping_set"]))
    put("MedAlvTwo", smv2["mapping_set_ratio_percentiles"]["Al"][2], ".2f")
    g0 = json.load(open(os.path.join(V, "gamma0_posthoc.json")))
    put("GzN", comma(g0["n"]))
    put("GzPct", 100 * g0["n"] / int(N["NMapped"].replace("{,}", "")), ".0f")
    for e in ["Mg", "Al"]:
        put(f"Gz{e}R", g0[e]["grs_r"], ".2f")
        N[f"Gz{e}CI"] = ci(g0[e]["grs_ci95"])
        put(f"Gz{e}Rel", g0[e]["reliability"], ".2f")
        put(f"Gz{e}Med", g0[e]["median"], ".2f")
    ss = pd.read_csv(os.path.join(V, "sites_stats_gamma0_posthoc.csv"))
    for e in ["Mg", "Al"]:
        r = ss[(ss.radius_deg == 2.0) & (ss.element == e)].iloc[0]
        put(f"GzSite{e}R", r.pearson_r, ".2f")
        N[f"GzSite{e}CI"] = ci(r.pearson_ci95)
    for e, rad in [("Al", 1.0), ("Al", 5.0), ("Mg", 1.0), ("Mg", 5.0)]:
        r = ss[(ss.radius_deg == rad) & (ss.element == e)].iloc[0]
        put(f"GzSite{e}R{'One' if rad == 1 else 'Five'}", r.pearson_r, ".2f")


XSM_MONTHS = ["2024_12", "2025_06", "2025_12"]


def xsm_numbers():
    X = os.path.join(A, "xsm")
    frames, nmin = [], 0
    for m in XSM_MONTHS:
        f = os.path.join(X, f"class_refits_{m}.parquet")
        if os.path.exists(f):
            d = pd.read_parquet(f)
            d["month"] = m
            frames.append(d)
            nmin += json.load(open(os.path.join(X, f"summary_{m}.json")))["n_xsm_minutes_used"]
    j = pd.concat(frames)
    xf = pd.concat([pd.read_parquet(os.path.join(X, f"xsm_fits_{m}.parquet")) for m in XSM_MONTHS
                    if os.path.exists(os.path.join(X, f"xsm_fits_{m}.parquet"))])
    put("XsmChiOne", np.median(xf.chi2_1t / xf.dof_1t), ".1f")
    put("XsmChiTwo", np.median(xf.chi2_2t / xf.dof_2t), ".1f")
    put("XsmNMonths", len(frames))
    put("XsmNClass", comma(len(j)))
    put("XsmNMin", comma(nmin))
    q = j.T1_xsm / j.goes_T_mk
    put("XsmTRatio", np.median(q), ".2f")
    put("XsmTRatioLo", np.percentile(q, 16), ".2f")
    put("XsmTRatioHi", np.percentile(q, 84), ".2f")
    put("XsmTmed", np.median(j.T1_xsm), ".1f")
    put("XsmTgoesMed", np.median(j.goes_T_mk), ".1f")
    for n, nm in [("xsm1t", "One"), ("xsm2t", "Two")]:
        ok = (j.quality_goes == 0) & (j[f"quality_{n}"] == 0)
        put(f"Xsm{nm}N", comma(ok.sum()))
        for e in ["Mg", "Al"]:
            r = 100 * (j.loc[ok, f"r_{e}_{n}"] / j.loc[ok, f"r_{e}_goes"] - 1)
            for p, s in [(16, "Lo"), (50, ""), (84, "Hi")]:
                put(f"Xsm{nm}{e}{s}", np.percentile(r, p), "+.0f")
                N[f"Xsm{nm}{e}{s}"] = N[f"Xsm{nm}{e}{s}"].replace("-", "$-$")
            for m, s in zip(XSM_MONTHS, "ABC"):  # A = 2024-12, B = 2025-06, C = 2025-12
                k = ok & (j.month == m)
                v = np.median(j.loc[k, f"r_{e}_{n}"] / j.loc[k, f"r_{e}_goes"]) * 100 - 100
                key = f"Xsm{nm}{e}M{s}"
                put(key, v, "+.0f")
                N[key] = N[key].replace("-", "$-$")
    # R2 Major 3: split by activity regime (flare enhancement = GOES XRS-B / 12-h running 10th percentile)
    regimes = [("Near baseline ($<1.5$)", 0, 1.5, "Lo"), ("Moderate enhancement (1.5--3)", 1.5, 3.0, "Mid"),
               (r"Strong enhancement ($\geq3$)", 3.0, np.inf, "Hi")]
    sg = lambda v: "0" if round(v) == 0 else f"{v:+.0f}".replace("-", "$-$")  # noqa: E731
    rows = []
    for lab, lo, hi, key in regimes:
        inr = (j.flare_enh >= lo) & (j.flare_enh < hi)
        cells = []
        # common subset: good fits with GOES, XSM 1T and XSM 2T spectra, so both columns share one N
        both = inr & (j.quality_goes == 0) & (j.quality_xsm1t == 0) & (j.quality_xsm2t == 0)
        for n, nm in [("xsm1t", "One"), ("xsm2t", "Two")]:
            ok = both
            put(f"XsmReg{key}{nm}N", comma(ok.sum()))
            for e in ["Mg", "Al"]:
                r = 100 * (j.loc[ok, f"r_{e}_{n}"] / j.loc[ok, f"r_{e}_goes"] - 1)
                v = np.median(r) if ok.sum() >= 10 else np.nan
                N[f"XsmReg{key}{nm}{e}"] = sg(v) if np.isfinite(v) else "--"
                cells.append((v, np.percentile(r, 16) if ok.sum() >= 10 else np.nan,
                              np.percentile(r, 84) if ok.sum() >= 10 else np.nan))
        ok1 = both
        tr = np.median(j.T1_xsm[ok1] / j.goes_T_mk[ok1]) if ok1.sum() >= 10 else np.nan
        put(f"XsmReg{key}TRatio", tr, ".2f")
        rows.append((lab, int(ok1.sum()), tr, np.median(j.goes_T_mk[ok1]) if ok1.sum() else np.nan, cells))
    put("XsmNBoth", comma(int(((j.quality_goes == 0) & (j.quality_xsm1t == 0) & (j.quality_xsm2t == 0)).sum())))
    fmt = lambda c: "--" if not np.isfinite(c[0]) else f"{sg(c[0])} ({sg(c[1])}, {sg(c[2])})"  # noqa: E731
    with open(os.path.join(P, "tables", "xsm_regimes.tex"), "w") as f:
        f.write("\\begin{tabular}{lrcccccc}\n\\toprule\n"
                " & & & & \\multicolumn{2}{c}{XSM isothermal} & \\multicolumn{2}{c}{XSM two-temperature}\\\\\n"
                "Regime ($F$) & $N$ & $T_\\mathrm{GOES}$ (MK) & $T_\\mathrm{XSM}/T_\\mathrm{GOES}$ & Mg/Si & Al/Si & Mg/Si & Al/Si\\\\\n"
                "\\midrule\n")
        for lab, n, tr, tg, c in rows:
            f.write(f"{lab} & {comma(n)} & {tg:.1f} & {tr:.2f} & " + " & ".join(fmt(x) for x in c) + "\\\\\n")
        f.write("\\bottomrule\n\\end{tabular}\n")


def epoch_numbers():
    """R1.2: per-year pass fractions for 2020-2023 and the two-epoch comparison."""
    import maps as mp
    # same criteria as the mapping set: quality 0, GOES >= B1, ratios not at a bound, >= 300 counts
    df = mp.load_fits(os.path.join(A, PM), months="2020_12:2023_12")
    t = df.groupby("month").agg(n=("use", "size"), q0=("use", "sum")).reset_index()
    t["year"] = t.month.str[:4]
    put("NArchEarly", len(t))
    for y, g in t.groupby("year"):
        nm = {"2020": "Twenty", "2021": "TwentyOne", "2022": "TwentyTwo", "2023": "TwentyThree"}[y]
        put(f"PassYr{nm}", 100 * g.q0.sum() / g.n.sum(), ".0f")
        put(f"NSpecYr{nm}", comma(g.n.sum()))
    t["frac"] = 100 * t.q0 / t.n
    for y in ["2021", "2022", "2023"]:
        g = t[t.year == y]
        nm = {"2021": "TwentyOne", "2022": "TwentyTwo", "2023": "TwentyThree"}[y]
        put(f"PassYr{nm}Lo", g.frac.min(), ".0f")
        put(f"PassYr{nm}Hi", g.frac.max(), ".0f")
    fr = lambda q, x: "0" if q == 0 else ("$<$0.1" if x < 0.05 else f"{x:.1f}")  # noqa: E731
    with open(os.path.join(P, "tables", "early_months.tex"), "w") as f:
        f.write("\\begin{tabular}{lrrr}\n\\toprule\nMonth & Dayside spectra & Pass & Fraction (\\%)\\\\\n\\midrule\n")
        for r in t.itertuples():
            f.write(f"{r.month.replace('_', '-')} & {comma(r.n)} & {comma(r.q0)} & {fr(r.q0, r.frac)}\\\\\n")
        f.write("\\bottomrule\n\\end{tabular}\n")
    e = json.load(open(os.path.join(A, "epochs", "epoch_compare.json")))
    for k, nm in [("A", "A"), ("B", "B")]:
        ep = e["epochs"][k]
        put(f"Ep{nm}NMonths", ep["n_months"])
        put(f"Ep{nm}NFp", comma(ep["n_footprints"]))
        put(f"Ep{nm}Tmed", ep["goes_T_mk_median"], ".1f")
        put(f"Ep{nm}Fwhm", ep["fwhm_Al_eV_median"], ".0f")
    for el in ["Mg", "Al"]:
        for c in e["compare"][el]:
            b = {2: "Two", 6: "Six", 10: "Ten"}[c["block_deg"]]
            put(f"Ep{el}{b}R", c["pearson_r"], ".2f")
            put(f"Ep{el}{b}N", comma(c["n_cells"]))
            put(f"Ep{el}{b}Ratio", c["median_ratio_B_over_A"], ".3f")
    V = os.path.join(A, "validation")
    for tag, nm in [("_v5_eA", "A"), ("_v5_eB", "B")]:
        g = pd.read_csv(os.path.join(V, f"grs_stats{tag}.csv"))
        s = pd.read_csv(os.path.join(V, f"sites_stats{tag}.csv"))
        for el in ["Mg", "Al"]:
            put(f"Ep{nm}Grs{el}R", g[g.element == el].iloc[0].pearson_r, ".2f")
            r = s[(s.radius_deg == 2.0) & (s.element == el)].iloc[0]
            put(f"Ep{nm}Site{el}R", r.pearson_r, ".2f")
            N[f"Ep{nm}Site{el}CI"] = ci(r.pearson_ci95)


def new_sites_numbers():
    """R2 Major 2: Chang'E-3/5/6 and Chandrayaan-3 sites (added after the pre-declared validation)."""
    from scipy.stats import pearsonr
    n = pd.read_csv(os.path.join(ROOT, "revision", f"new_sites{PT}.csv"))
    n = n[n.radius_deg == 2.0]
    s = pd.read_csv(os.path.join(A, "validation", f"sites_table{PT}.csv"))
    s = s[s.radius_deg == 2.0]
    key = {"Chang'E-3": "Ng2023", "Chang'E-5": "Li2022CE5", "Chang'E-6": "Li2024CE6", "Chandrayaan-3": "Vadawale2024"}
    kind = {"Chang'E-3": "APXS", "Chang'E-5": "returned soil", "Chang'E-6": "returned soil", "Chandrayaan-3": "APXS"}
    for e in ["Mg", "Al"]:
        a, b = s[s.element == e], n[n.element == e]
        x = np.r_[a.class_ratio, b.class_ratio]
        y = np.r_[a.ref_ratio, b.ref_ratio]
        put(f"NewSite{e}R", pearsonr(b.class_ratio, b.ref_ratio)[0], ".2f")
        put(f"AllSite{e}R", pearsonr(x, y)[0], ".2f")
        put(f"AllSite{e}Scale", np.median(x / y), ".2f")
    put("NewSiteFpMin", comma(n.n_footprints.min()))
    put("NewSiteFpMax", comma(n.n_footprints.max()))
    with open(os.path.join(P, "tables", "new_sites.tex"), "w") as f:
        f.write("\\begin{tabular}{lrrlrrrr}\n\\toprule\nSite & Lat & Lon & Reference & \\multicolumn{2}{c}{Mg/Si} & "
                "\\multicolumn{2}{c}{Al/Si}\\\\\n & (\\textdegree) & (\\textdegree E) & & ref. & CLASS & ref. & CLASS\\\\\n\\midrule\n")
        for site in key:
            r = n[n.site == site].set_index("element")
            f.write(f"{site} & {r.lat.iloc[0]:.2f} & {r.lon.iloc[0]:.2f} & {kind[site]}, \\citet{{{key[site]}}} & "
                    f"{r.loc['Mg', 'ref_ratio']:.3f} & {r.loc['Mg', 'class_ratio']:.3f}$\\pm${r.loc['Mg', 'class_se']:.3f} & "
                    f"{r.loc['Al', 'ref_ratio']:.3f} & {r.loc['Al', 'class_ratio']:.3f}$\\pm${r.loc['Al', 'class_se']:.3f}\\\\\n")
        f.write("\\bottomrule\n\\end{tabular}\n")


def goes1s_numbers():
    """R1.3: exposure-matched GOES 1 s (flx1s, rescaled to avg1m per minute) vs nearest 1 min."""
    j = pd.read_parquet(os.path.join(A, "months_v5_g1s", "compare.parquet"))
    put("GosN", comma(len(j)))
    dT = 100 * (j.goes_T_mk_s / j.goes_T_mk_b - 1)
    for p_, nm in [(50, "Med"), (5, "Lo"), (95, "Hi")]:
        put(f"GosT{nm}", np.percentile(dT, p_), "+.1f")
    big = dT.abs() > 10
    put("GosBigPct", 100 * big.mean(), ".0f")
    for e in ["Mg", "Al"]:
        d = 100 * (j[f"r_{e}_s"] / j[f"r_{e}_b"] - 1)
        for sub, sn in [(slice(None), ""), (big, "Big")]:
            x = d[sub]
            for p_, nm in [(50, "Med"), (16, "Lo"), (84, "Hi")]:
                put(f"Gos{e}{sn}{nm}", np.percentile(x, p_), "+.1f")
    for k in list(N):
        if k.startswith("Gos") and isinstance(N[k], str) and N[k].startswith("-"):
            N[k] = "$-$" + N[k][1:]


def main():
    base = os.path.join(ROOT, "products", "headline_numbers.json")
    old = json.load(open(base)) if os.path.exists(base) else {}
    for fn in [header_numbers, stage_numbers, calib_numbers, background_numbers, solar_numbers,
               synthetic_numbers, validation_numbers, tsens_numbers, diag_numbers, xsm_numbers, epoch_numbers, goes1s_numbers, new_sites_numbers]:
        try:
            fn()
        except Exception as ex:  # noqa: BLE001
            print(f"WARNING {fn.__name__}: {ex}")
    # Values whose inputs are not available here keep their previous (v3) value, and are
    # listed so that every carried-over number is reviewed before submission.
    stale = sorted(k for k in old if k not in N)
    for k in stale:
        N[k] = old[k]
    with open(os.path.join(A, "headline_numbers_stale.txt"), "w") as f:
        f.write("\n".join(stale) + "\n")
    print(f"{len(N) - len(stale)} numbers recomputed, {len(stale)} carried over from {base}:")
    print(" ".join(stale))
    with open(os.path.join(A, "headline_numbers.json"), "w") as f:
        json.dump(N, f, indent=1, sort_keys=True)
    with open(os.path.join(P, "numbers.tex"), "w") as f:
        f.write("% AUTO-GENERATED by src/asr/make_numbers.py -- do not edit\n")
        for k, v in sorted(N.items()):
            v = "0" if str(v) in ("$-$0", "-0", "+0") else v
            f.write(f"\\newcommand{{\\{k}}}{{{v}}}\n")


if __name__ == "__main__":
    main()
