"""Revision 1 text-only edits (R1-1, 4, 7, 9, 10, 11; R2-1 citations, R2 minor 1, 2, 4). Run once."""
import os
import sys

S = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "paper_asr", "sections")
E = {
    "introduction.tex": [
        (r"The fluorescence yield of each element depends on the solar flux above its absorption edge, so a change in the solar spectrum can change line ratios such as Mg/Si or Ca/Si even when the surface composition is unchanged.",
         r"The fluorescence line flux of each element depends on the solar flux above its K absorption edge, so a change in the solar spectrum can change line-flux ratios such as Mg/Si or Ca/Si even when the surface composition is unchanged."),
        (r"maps of X-ray line-flux ratios have also been presented \citep{Kumar2025}.",
         r"maps of X-ray line-flux ratios have also been presented \citep{Kumar2026}. Global Mg/Si and Al/Si maps were first derived from the Chang'E-2 X-ray spectrometer \citep{Dong2016}, and oxide maps calibrated against returned samples, including those of Chang'E-5, have been derived from Kaguya Multiband Imager reflectance \citep{Wang2021,Zhang2023}."),
    ],
    "model.tex": [
        (r"Minor oxides are fixed at typical soil values (Na$_2$O 0.5, K$_2$O 0.1, MnO 0.15, Cr$_2$O$_3$ 0.2 wt\%).",
         r"Minor oxides are fixed at typical soil values (Na$_2$O 0.5, K$_2$O 0.1, MnO 0.15, Cr$_2$O$_3$ 0.2 wt\%), close to the average soils of the Apollo and Luna landing sites (Na$_2$O 0.31--0.70, K$_2$O 0.04--0.55, MnO 0.13--0.30, Cr$_2$O$_3$ 0.23--0.36 wt\%; \citealp{Heiken1991}, Table 7.15). Varying them over these ranges changes $C_\mathrm{Si}$ in \cref{eq:closure} by less than 1\%."),
        (r"Instead, the scattered component carries its own free amplitude $B$.",
         r"Instead, the scattered component carries its own free amplitude $A_\mathrm{sc}$ (\cref{eq:model})."),
        (r"We interpret it as solar X-ray fluorescence of aluminium in the instrument structure. If present on the dayside at the nightside level, it would be comparable in magnitude to the lunar Al line.",
         r"We interpret it as solar X-ray fluorescence of aluminium in the instrument (filters and structure). In sunlit nightside spectra it is therefore significant; were it present on the dayside at the same level, it would be comparable in magnitude to the lunar Al line."),
        (r"The geometry of \cref{fig:background} suggests that the instrumental component requires the Sun to be behind the spacecraft: it is already small at \texttt{SOLARANG} $\approx90$--93\textdegree, and we therefore expect it to be negligible on the dayside, where the Sun is in front of the nadir-pointing instrument.",
         r"The geometry of \cref{fig:background} indicates when the component arises. For dayside nadir-pointing observations (\texttt{SOLARANG} $<90^\circ$), the Sun is above the plane of the instrument's Moon-facing (aperture) side, i.e.\ behind the instrument relative to its viewing direction, and that side, including the outer filter, is in the shadow of the spacecraft. Only when the Sun is below that plane (\texttt{SOLARANG} $>90^\circ$) and the spacecraft is still sunlit ($\lesssim\SI{109}{\degree}$) can sunlight reach the aperture side directly. The observed onset near 93\textdegree{} and the disappearance at the shadow boundary agree with this geometry. We therefore expect the component to be negligible on the dayside, although it is significant in sunlit nightside spectra."),
    ],
    "fitting.tex": [
        (r"A\Big[\sum_\ell \frac{P_\ell}{P_\mathrm{ref}} R_{\ell k} + B\,\frac{(\mathcal{F}\Sigma)_k}{\Sigma_\mathrm{ref}}\Big]",
         r"A\Big[\sum_\ell \frac{P_\ell}{P_\mathrm{ref}} R_{\ell k} + A_\mathrm{sc}\,\frac{(\mathcal{F}\Sigma)_k}{\Sigma_\mathrm{ref}}\Big]"),
        (r"The parameter vector is $\theta=(A,r_\mathrm{Mg},r_\mathrm{Al},r_\mathrm{Ca},B,\beta)$.",
         r"Here $A$ is the overall amplitude, $A_\mathrm{sc}$ the relative amplitude of the scattered solar spectrum (\cref{eq:scatter}) and $\beta$ the particle-background scale; the GOES flux classes (B1, B3) quoted elsewhere are unrelated. The parameter vector is $\theta=(A,r_\mathrm{Mg},r_\mathrm{Al},r_\mathrm{Ca},A_\mathrm{sc},\beta)$."),
        ("with a trust-region reflective least-squares solver applied to the signed deviance residuals, whose sum of squares equals $C$.",
         "with a trust-region reflective least-squares solver applied to the signed deviance residuals\n\\begin{equation}\n\\delta_k = \\mathrm{sign}(d_k-m_k)\\,\\sqrt{2\\big[m_k-d_k+d_k\\ln(d_k/m_k)\\big]},\n\\label{eq:devres}\n\\end{equation}\nwhose sum of squares equals $C$ (the logarithmic term is zero for $d_k=0$). The residuals plotted in \\cref{fig:fits} are these $\\delta_k$."),
        (r"core block $(A,r_\mathrm{Mg},r_\mathrm{Al},B,\beta)$",
         r"core block $(A,r_\mathrm{Mg},r_\mathrm{Al},A_\mathrm{sc},\beta)$"),
        (r"The products are Mg/Si and Al/Si weight ratios with per-footprint uncertainties.",
         r"The products are Mg/Si and Al/Si elemental weight ratios, $r_i=C_i/C_\mathrm{Si}$, with per-footprint uncertainties. They are not line-flux ratios: the fit converts the observed line fluxes into concentrations through the forward model of \cref{sec:model}, which includes the solar spectrum of each exposure, the incidence geometry, matrix effects and the instrument response. Line-flux ratios are not used at any stage of the mapping or the validation."),
    ],
    "discussion.tex": [
        (r"Refitting real data with the GOES temperature scaled by 0.8 and 1.2 changes the median Mg/Si by \TsensMgLo\% and \TsensMgHi\%, and Al/Si by \TsensAlLo\% and \TsensAlHi\%.",
         r"Refitting real data with the GOES temperature scaled by 0.8 and 1.2 changes the median Mg/Si by \TsensMgLo\% and \TsensMgHi\%, and Al/Si by \TsensAlLo\% and \TsensAlHi\%."
         "\n\n"
         r"Mg/Si is about three times as sensitive as Al/Si because of where the K edges lie relative to the solar lines. Photons above the Si K edge (\SI{1.84}{keV}) excite all three elements, so changes in that part of the spectrum largely cancel in both ratios. The ratios therefore respond mainly to photons that excite the numerator element but not Si. For Al this band is narrow (\SIrange{1.56}{1.84}{keV}) and contains only weak lines (e.g.\ Al\,\textsc{xii}). For Mg it is wider (\SIrange{1.30}{1.84}{keV}) and contains the strong Mg\,\textsc{xi} and Mg\,\textsc{xii} lines, whose emission relative to the Si\,\textsc{xiii} and S\,\textsc{xv} complexes changes with temperature."),
        (r"The most likely explanation is solar X-rays reaching aluminium in the instrument (filters, frames, housing; \citealp{Pillai2021}) when the Sun is behind the spacecraft. On the dayside the Sun is in front of the nadir-pointing instrument, and we expect the component to be negligible.",
         r"The most likely explanation is solar X-rays reaching aluminium in the instrument (filters, frames, housing; \citealp{Pillai2021}) when the Sun is below the plane of the instrument's Moon-facing side and can illuminate it directly (\cref{sec:background}). On the dayside the Sun is behind the instrument relative to its viewing direction, that side is shadowed, and we expect the component to be negligible. It is significant only in sunlit nightside spectra. It matters for dayside Al/Si because nightside spectra are used to build the background model, as follows."),
        (r"\item \emph{Relative, not absolute.} The products are weight ratios to Si. No absolute calibration is attempted, and the amplitude $A$ absorbs the solar flux, area and geometry.",
         r"\item \emph{Ratios, not absolute abundances.} The products are elemental weight ratios to Si derived through the forward model; they are neither line-flux ratios nor absolute concentrations. No absolute flux calibration is used: the amplitude $A$ absorbs the solar flux, area and geometry. Errors in the relative efficiency of the Mg, Al and Si lines would rescale the ratios (item 3)."),
    ],
    "conclusions.tex": [
        (r"\item CLASS nightside spectra contain a solar-driven instrumental Al K$\alpha$ component when the Sun illuminates the spacecraft from behind. How it is treated in the dayside background model affects Al/Si.",
         r"\item CLASS nightside spectra contain a solar-driven instrumental Al K$\alpha$ component while the spacecraft is sunlit and the Sun can illuminate the instrument's Moon-facing side. It is expected to be negligible on the dayside, but its treatment in the background model affects Al/Si."),
    ],
    "abstract.tex": [
        (r"an instrumental Al K$\alpha$ component that appears when the spacecraft is sunlit from behind and whose treatment",
         r"an instrumental Al K$\alpha$ component that appears in sunlit nightside spectra and whose treatment"),
    ],
    "validation.tex": [
        (r"The history of the primary model, including an earlier model whose validation motivated it, is given in \cref{sec:alt}.",
         r"The history of the primary model, including an earlier model whose validation motivated it, is given in \cref{sec:alt}. Both comparisons use elemental weight ratios on both sides: the CLASS values are concentrations derived through the forward model (\cref{sec:products}), and the reference oxide abundances are converted with stoichiometric factors (\cref{sec:refdata}). No numerical equality between line-flux ratios and abundance ratios is assumed."),
        (r"\caption{Landing-site validation at 2\textdegree{} radius: returned-soil weight ratios and CLASS weighted means $\pm$ standard error.}",
         r"\caption{Landing-site validation at 2\textdegree{} radius: returned-soil element weight ratios and CLASS element weight ratios (weighted means $\pm$ standard error).}"),
    ],
}

for fn, reps in E.items():
    p = os.path.join(S, fn)
    t = open(p, encoding="utf-8", newline="").read()
    crlf = "\r\n" in t
    t = t.replace("\r\n", "\n")
    for old, new in reps:
        if t.count(old) != 1:
            sys.exit(f"{fn}: expected 1 match, found {t.count(old)}: {old[:70]}")
        t = t.replace(old, new)
    if crlf:
        t = t.replace("\n", "\r\n")
    open(p, "w", encoding="utf-8", newline="").write(t)
    print("edited", fn, len(reps))
