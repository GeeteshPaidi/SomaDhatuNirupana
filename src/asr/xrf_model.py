"""
Simplified fundamental-parameter forward model for CLASS lunar XRF spectra
(ASR rework, PLAN.md section 4).

Physics (thick, flat, homogeneous target; photon counts per channel)
---------------------------------------------------------------------
Primary fluorescence of line l of analyte i (Sherman 1955, simplified):

    P_il = C_i * w_iK * p_il * Int_{E > E_K,i} S(E) tau_iK(E)
                   / [ mu(E)/sin(psi1) + mu(E_il)/sin(psi2) ] dE

  C_i      weight fraction of element i
  w_iK     K-shell fluorescence yield            (xraylib FluorYield)
  p_il     radiative rate of line l              (xraylib RadRate)
  tau_iK   K-shell photoionisation cross-section (xraylib CS_Photo_Partial, cm^2/g)
  mu(E)    matrix mass attenuation  sum_j C_j mu_j(E)  (xraylib CS_Total)
  S(E)     incident solar photon spectrum (solar_spectrum below)
  psi1     incidence elevation = 90 deg - SOLARANG; psi2 = 90 deg (nadir)

Secondary fluorescence of i excited by line m of element j (thick-target
formula of de Boer 1990 / Mantler; concentration-coupled through C_j and mu):

    S_il,jm = 0.5 * C_i w_iK p_il * C_j w_jK p_jm * tau_iK(E_jm)
              * Int S(E) tau_jK(E) / [mu(E)/sin psi1 + mu(E_il)/sin psi2]
                * [ sin psi1/mu(E) ln(1 + mu(E)/(sin psi1 mu(E_jm)))
                  + sin psi2/mu(E_il) ln(1 + mu(E_il)/(sin psi2 mu(E_jm))) ] dE
    (only for E_jm > E_K,i)

Scattered solar continuum (coherent + incoherent; angular dependence of the
differential cross-section absorbed into a free amplitude):

    Sc(E) = S(E) sigma_s(E) / [ mu(E)/sin psi1 + mu(E)/sin psi2 ]

Detector: each line is a Gaussian integrated over each channel's energy
bin (erf difference), truncated at +-5 sigma and renormalised to unit area,
sigma(E)^2 = sigma0^2 + F w E (Si Fano F = 0.115, w = 3.68 eV); optional
energy-dependent efficiency (unity by default: no CLASS efficiency data
ship with the L1 release).

Composition closure: elements enter as weight ratios r_i = C_i/C_Si; all
fitted elements plus fixed minor oxides are converted to oxides and closed
to 100 wt%, which fixes C_Si and stoichiometric O (O is not a free
parameter).
"""
import numpy as np
import xraylib
from scipy.special import erf

# Elements whose ratios to Si are fitted.
FIT_ELEMENTS = ["Mg", "Al", "Ca", "Ti", "Fe"]
LINE_ELEMENTS = ["Mg", "Al", "Si", "Ca", "Ti", "Fe"]
# Oxide stoichiometry: element -> (oxide mass / element mass), O mass per element mass.
_OXIDES = {"Si": ("SiO2", 1, 2), "Mg": ("MgO", 1, 1), "Al": ("Al2O3", 2, 3),
           "Ca": ("CaO", 1, 1), "Ti": ("TiO2", 1, 2), "Fe": ("FeO", 1, 1),
           "Na": ("Na2O", 2, 1), "K": ("K2O", 2, 1), "Mn": ("MnO", 1, 1), "Cr": ("Cr2O3", 2, 3)}
# Fixed minor oxides (wt fraction), typical lunar soil values; matrix only (no lines).
MINOR_OXIDES = {"Na": 0.005, "K": 0.001, "Mn": 0.0015, "Cr": 0.002}

M_O = xraylib.AtomicWeight(8)
FANO = 0.115
W_PAIR = 3.68e-3  # keV per electron-hole pair in Si
KEV_PER_MK = 0.08617


def oxygen_per_element(el):
    _, n_el, n_o = _OXIDES[el]
    return n_o * M_O / (n_el * xraylib.AtomicWeight(xraylib.SymbolToAtomicNumber(el)))


def composition_from_ratios(ratios):
    """ratios: dict El -> C_El/C_Si (weight) for FIT_ELEMENTS.
    Returns dict El -> weight fraction (incl. O and fixed minors), summing to 1."""
    minor_ox = sum(MINOR_OXIDES.values())
    oxide_per_si = 1 + oxygen_per_element("Si") + sum(
        ratios[e] * (1 + oxygen_per_element(e)) for e in FIT_ELEMENTS)
    c_si = (1 - minor_ox) / oxide_per_si
    comp = {"Si": c_si}
    for e in FIT_ELEMENTS:
        comp[e] = ratios[e] * c_si
    for e, ox in MINOR_OXIDES.items():
        comp[e] = ox / (1 + oxygen_per_element(e))
    comp["O"] = sum(comp[e] * oxygen_per_element(e) for e in list(comp))
    return comp


def solar_spectrum(E, kT=1.07):
    """Isothermal thermal-bremsstrahlung photon spectrum, Gaunt factor = 1:
    S(E) ~ exp(-E/kT)/E (arbitrary normalisation; only the shape matters,
    the amplitude is a free fit parameter). kT in keV; default 1.07 keV is
    the single-temperature fit to a C3 flare measured by Chandrayaan-1 XSM
    (Narendranath et al. 2013, quoted in arXiv:1305.6218)."""
    return np.exp(-E / kT) / E


class Line:
    def __init__(self, el, line_id, label):
        self.el, self.label = el, label
        self.Z = xraylib.SymbolToAtomicNumber(el)
        self.E = xraylib.LineEnergy(self.Z, line_id)
        self.Q = xraylib.FluorYield(self.Z, xraylib.K_SHELL) * xraylib.RadRate(self.Z, line_id)
        self.edge = xraylib.EdgeEnergy(self.Z, xraylib.K_SHELL)


def build_lines(elements=LINE_ELEMENTS):
    lines = []
    for el in elements:
        for lid, lab in [(xraylib.KA_LINE, "Ka"), (xraylib.KB_LINE, "Kb")]:
            try:
                lines.append(Line(el, lid, lab))
            except ValueError:  # Kb absent for Na/Mg in xraylib
                pass
    return lines


class ForwardModel:
    """Precomputes atomic data on an excitation grid; evaluates line
    intensities, scattered continuum and the folded channel spectrum."""

    def __init__(self, kT=1.07, e_max=20.0, de=0.01, secondary=True,
                 efficiency=None, lines=None):
        self.kT = kT
        self.secondary = secondary
        self.lines = lines if lines is not None else build_lines()
        self.efficiency = efficiency  # callable E -> efficiency, or None for unity
        self.elements = sorted(set(LINE_ELEMENTS) | set(MINOR_OXIDES) | {"O"})
        # Excitation grid starts just above the lowest edge (Mg K, 1.305 keV).
        self.Eg = np.arange(1.30, e_max, de) + de / 2
        self.de = de
        self.S = solar_spectrum(self.Eg, kT)
        Zs = {e: xraylib.SymbolToAtomicNumber(e) for e in self.elements}
        self.mu_g = {e: np.array([xraylib.CS_Total(Zs[e], x) for x in self.Eg]) for e in self.elements}
        self.sc_g = {e: np.array([xraylib.CS_Rayl(Zs[e], x) + xraylib.CS_Compt(Zs[e], x) for x in self.Eg])
                     for e in self.elements}
        self.tau_g = {}
        for e in LINE_ELEMENTS:
            Z = Zs[e]
            self.tau_g[e] = np.array([xraylib.CS_Photo_Partial(Z, xraylib.K_SHELL, x)
                                      if x > xraylib.EdgeEnergy(Z, xraylib.K_SHELL) else 0.0
                                      for x in self.Eg])
        # mu_j at each line energy, tau_iK at each enhancer line energy
        self.mu_line = {e: np.array([xraylib.CS_Total(Zs[e], ln.E) for ln in self.lines]) for e in self.elements}
        self.tau_at_line = {}
        for i, li in enumerate(self.lines):
            for m, lm in enumerate(self.lines):
                if lm.E > li.edge:
                    self.tau_at_line[(i, m)] = xraylib.CS_Photo_Partial(li.Z, xraylib.K_SHELL, lm.E)
        # Stacked arrays for the vectorised primary-only path used in fitting.
        self.MU_G = np.array([self.mu_g[e] for e in self.elements])          # (n_el, nE)
        self.MU_L = np.array([self.mu_line[e] for e in self.elements])       # (n_el, n_lines)
        self.S_TAU = np.array([self.S * self.tau_g[ln.el] for ln in self.lines])  # (n_lines, nE)
        self.Q = np.array([ln.Q for ln in self.lines])
        self.line_el = [ln.el for ln in self.lines]

    def comp_vector(self, comp):
        return np.array([comp[e] for e in self.elements])

    def primary_fast(self, ratios, s1, s2=1.0):
        """Vectorised primary intensities (no secondary); same result as
        line_intensities(..., secondary=False)."""
        comp = composition_from_ratios(ratios)
        c = self.comp_vector(comp)
        mu_g = c @ self.MU_G
        mu_l = c @ self.MU_L
        denom = mu_g[None, :] / s1 + mu_l[:, None] / s2
        integ = np.sum(self.S_TAU / denom, axis=1) * self.de
        conc = np.array([comp[e] for e in self.line_el])
        return conc * self.Q * integ, comp

    # -- matrix --------------------------------------------------------
    def _mu(self, comp):
        mu_g = sum(comp[e] * self.mu_g[e] for e in self.elements)
        mu_l = sum(comp[e] * self.mu_line[e] for e in self.elements)
        return mu_g, mu_l

    def line_intensities(self, ratios, psi1_deg, psi2_deg=90.0, return_parts=False):
        """Photon intensity of each line (arbitrary units, common scale)."""
        comp = composition_from_ratios(ratios)
        s1, s2 = np.sin(np.radians(psi1_deg)), np.sin(np.radians(psi2_deg))
        mu_g, mu_l = self._mu(comp)
        prim = np.zeros(len(self.lines))
        sec = np.zeros(len(self.lines))
        for i, li in enumerate(self.lines):
            denom = mu_g / s1 + mu_l[i] / s2
            prim[i] = comp[li.el] * li.Q * np.sum(self.S * self.tau_g[li.el] / denom) * self.de
            if not self.secondary:
                continue
            for m, lm in enumerate(self.lines):
                key = (i, m)
                if key not in self.tau_at_line or lm.el == li.el:
                    continue
                mu_jm = mu_l[m]
                geom = (s1 / mu_g) * np.log1p(mu_g / (s1 * mu_jm)) + \
                       (s2 / mu_l[i]) * np.log1p(mu_l[i] / (s2 * mu_jm))
                sec[i] += 0.5 * comp[li.el] * li.Q * comp[lm.el] * lm.Q * self.tau_at_line[key] * \
                    np.sum(self.S * self.tau_g[lm.el] / denom * geom) * self.de
        if return_parts:
            return prim, sec, comp
        return prim + sec

    def scatter_continuum(self, ratios, E, psi1_deg, psi2_deg=90.0):
        """Scattered solar photon spectrum per keV evaluated at energies E (keV).
        Uses the full solar spectrum (not only above the Mg edge)."""
        comp = composition_from_ratios(ratios)
        s1, s2 = np.sin(np.radians(psi1_deg)), np.sin(np.radians(psi2_deg))
        Zs = {e: xraylib.SymbolToAtomicNumber(e) for e in self.elements}
        key = tuple(np.round(E, 6))
        if getattr(self, "_sc_cache_key", None) != key:
            self._sc_cache_key = key
            self._sc_mu = {e: np.array([xraylib.CS_Total(Zs[e], x) for x in E]) for e in self.elements}
            self._sc_s = {e: np.array([xraylib.CS_Rayl(Zs[e], x) + xraylib.CS_Compt(Zs[e], x) for x in E])
                          for e in self.elements}
        mu = sum(comp[e] * self._sc_mu[e] for e in self.elements)
        sig = sum(comp[e] * self._sc_s[e] for e in self.elements)
        return solar_spectrum(E, self.kT) * sig / (mu / s1 + mu / s2)


class Detector:
    """Channel grid + Gaussian line response (+-5 sigma, erf-binned, renormalised)."""

    def __init__(self, gain=0.0135, offset=0.0, sigma0=0.045, n_channels=2048, n_sigma=5.0):
        self.gain, self.offset, self.sigma0 = gain, offset, sigma0
        self.n_channels, self.n_sigma = n_channels, n_sigma
        ch = np.arange(n_channels)
        self.E_center = gain * ch + offset
        self.E_lo = self.E_center - gain / 2
        self.E_hi = self.E_center + gain / 2

    def sigma(self, E):
        return np.sqrt(self.sigma0 ** 2 + FANO * W_PAIR * E)

    def line_response(self, E_line):
        """Fraction of a line's photons landing in each channel (sums to 1
        over the full +-n_sigma window when it lies inside the channel grid)."""
        s = self.sigma(E_line)
        lo = np.clip((self.E_lo - E_line) / s, -self.n_sigma, self.n_sigma)
        hi = np.clip((self.E_hi - E_line) / s, -self.n_sigma, self.n_sigma)
        frac = 0.5 * (erf(hi / np.sqrt(2)) - erf(lo / np.sqrt(2)))
        return frac / erf(self.n_sigma / np.sqrt(2))

    def response_matrix(self, lines):
        return np.array([self.line_response(ln.E) for ln in lines])
