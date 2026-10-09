"""
Poisson (Cash) fit of raw CLASS counts with the simplified FP model
(ASR rework, PLAN.md section 4; TRACKER tasks 1.5, 1.6).

Model counts in channel k (fit range E_lo..E_hi, default 1.15-7.5 keV):

    m_k = A * eff(E_k) * [ sum_l I_l(r, psi1) R_lk + B * Sc_k(r, psi1) ]
          + beta * b_k * t * n_SCD/16

  r = (Mg, Al, Ca, Ti, Fe)/Si weight ratios (C_Si fixed by oxide closure,
      O stoichiometric -> not free), A >= 0 overall amplitude,
  B >= 0 scattered-solar-continuum amplitude (relative to A),
  b_k nightside background rate template (counts/s, 16-SCD equivalent),
  beta >= 0 background scale, t exposure (s).

Minimises the Cash statistic C = 2 sum[m - d + d ln(d/m)] via
scipy least_squares on signed deviance residuals (sum of squares = C).
Uncertainties: inverse Fisher information F = J_m^T diag(1/m) J_m.
Goodness of fit: (C - E[C]) / sqrt(Var[C]) with exact per-channel Poisson
expectations (Kaastra 2017).
"""
import numpy as np
import xraylib
from scipy import optimize
from scipy.stats import poisson

from xrf_model import FIT_ELEMENTS, Detector, ForwardModel, composition_from_ratios, solar_spectrum

PARAMS = ["A"] + [f"r_{e}" for e in FIT_ELEMENTS] + ["B", "beta", "gamma"]
RATIO_UPPER = 5.0

# Quality-flag bits
Q_NOCONV = 1
Q_LOWCOUNTS = 2
Q_BADGOF = 4
Q_ATBOUND = 8  # Mg/Si or Al/Si at a bound (per-element flags in atbound_<El>)
Q_ILLCOND = 16
Q_BADGEOM = 32
Q_COVFAIL = 64


def _cstat_terms(d, m):
    with np.errstate(divide="ignore", invalid="ignore"):
        t = m - d + np.where(d > 0, d * np.log(d / m), 0.0)
    return 2.0 * t


class _CExpect:
    """Tabulated exact E[C_k] and Var[C_k] for a Poisson channel of mean m."""

    def __init__(self):
        grid = np.concatenate([np.logspace(-4, 1, 120), np.linspace(10.5, 200, 120)])
        ev, vv = [], []
        for m in grid:
            kmax = int(m + 12 * np.sqrt(m) + 25)
            k = np.arange(kmax + 1)
            p = poisson.pmf(k, m)
            c = _cstat_terms(k.astype(float), np.full(k.shape, m))
            e = np.sum(p * c)
            ev.append(e)
            vv.append(np.sum(p * (c - e) ** 2))
        self.lg = np.log(grid)
        self.e = np.array(ev)
        self.v = np.array(vv)

    def __call__(self, m):
        lm = np.log(np.clip(m, 1e-4, 200))
        e = np.interp(lm, self.lg, self.e)
        v = np.interp(lm, self.lg, self.v)
        big = m > 200  # asymptotic chi^2 limit
        e[big], v[big] = 1.0, 2.0
        return e, v


_CEXP = None


def cstat_expectation(m):
    global _CEXP
    if _CEXP is None:
        _CEXP = _CExpect()
    return _CEXP(m)


class SpectrumFitter:
    def __init__(self, fm: ForwardModel, det: Detector, bkg_rate, e_lo=1.15, e_hi=7.5,
                 min_counts=300, gof_z_max=3.0, cond_max=1e4, inst_rate=None, hi_band=(8.0, 20.0)):
        """bkg_rate: particle-background template P (counts/s/16 SCD).
        inst_rate: solar-driven instrumental template T (same units, unit mean
        level), amplitude gamma free; None -> gamma fixed at 0.
        hi_band: background-dominated band added to the fit as ONE summed bin
        to pin beta (shape-independent); None to disable."""
        self.fm, self.det = fm, det
        Ec = det.E_center
        self.main_sel = (Ec >= e_lo) & (Ec <= e_hi)
        self.sel = self.main_sel.copy()
        if hi_band is not None:
            self.sel |= (Ec >= hi_band[0]) & (Ec < hi_band[1])
        self.main = self.main_sel[self.sel]  # mask within fitted channels
        self.n_main = int(self.main.sum())   # main channels come first (energy order)
        self.E = Ec[self.sel]
        self.use_inst = inst_rate is not None
        self.inst_rate = (np.asarray(inst_rate, float) if self.use_inst else np.zeros(det.n_channels))[self.sel]
        self.R = det.response_matrix(fm.lines)[:, self.sel]
        self.eff = fm.efficiency(self.E) if fm.efficiency is not None else np.ones_like(self.E)
        self.bkg_rate = np.asarray(bkg_rate, float)[self.sel]
        self.min_counts, self.gof_z_max, self.cond_max = min_counts, gof_z_max, cond_max
        # Normalise intensities so A ~ "Si Ka counts per unit" for a reference composition.
        ref = {"Mg": 0.25, "Al": 0.4, "Ca": 0.4, "Ti": 0.03, "Fe": 0.4}
        I, _ = fm.primary_fast(ref, 1.0)
        self.i_norm = I.max()
        # matrix attenuation / scattering cross-sections at channel energies
        Zs = [xraylib.SymbolToAtomicNumber(e) for e in fm.elements]
        self.MU_E = np.array([[xraylib.CS_Total(Z, x) for x in self.E] for Z in Zs])
        self.SIG_E = np.array([[xraylib.CS_Rayl(Z, x) + xraylib.CS_Compt(Z, x) for x in self.E] for Z in Zs])
        self.S_E = solar_spectrum(self.E, fm.kT)
        self.sc_norm = 1.0
        self.sc_norm = self._scatter(composition_from_ratios(ref), 1.0).max()

    # -- model ---------------------------------------------------------
    def _scatter(self, comp, s1):
        """Scattered solar continuum per channel, nadir emergence (sin psi2 = 1)."""
        c = self.fm.comp_vector(comp)
        mu = c @ self.MU_E
        sig = c @ self.SIG_E
        return self.S_E * sig / (mu / s1 + mu) * self.det.gain / self.sc_norm

    def source(self, theta, s1):
        A, r, B = theta[0], dict(zip(FIT_ELEMENTS, theta[1:6])), theta[6]
        I, comp = self.fm.primary_fast(r, s1)
        lines = (I / self.i_norm) @ self.R
        sc = self._scatter(comp, s1)
        return A * self.eff * (lines + B * sc)

    def _agg(self, v):
        """Main-window channels individually + hi band summed into one bin."""
        if self.n_main == v.size:
            return v
        return np.r_[v[self.main], v[~self.main].sum()]

    def background_counts(self, exposure, n_scd):
        f = exposure * n_scd / 16.0
        return np.vstack([self.bkg_rate * f, self.inst_rate * f])

    def model(self, theta, s1, bkg):
        """bkg: 2 x n array from background_counts (particle, instrumental)."""
        return self.source(theta, s1) + theta[7] * bkg[0] + theta[8] * bkg[1]

    # -- fit -----------------------------------------------------------
    def fit(self, counts, exposure, n_scd, solarang, x0=None):
        d = self._agg(np.asarray(counts, float)[self.sel])
        psi1 = 90.0 - solarang
        out = {"n_counts_fit": float(d[:self.n_main].sum()), "n_counts_hiband": float(d[self.n_main:].sum()),
               "psi1": psi1}
        flag = 0
        if not (0.5 < psi1 <= 90.0) or exposure <= 0 or n_scd <= 0:
            out.update(quality=Q_BADGEOM, success=False)
            return out
        s1 = np.sin(np.radians(psi1))
        bkg = self.background_counts(exposure, n_scd)
        if x0 is None:
            excess = max(d[:self.n_main].sum() - bkg[0][self.main].sum(), 1.0)
            x0 = np.array([excess, 0.25, 0.4, 0.4, 0.03, 0.4, 1.0, 1.0, 0.3 if self.use_inst else 0.0])
            # scale A so the line part carries ~70% of the excess
            m_line = self.source(np.r_[1.0, x0[1:6], 0.0], s1).sum()
            x0[0] = 0.7 * excess / max(m_line, 1e-12)
        lb = np.r_[0.0, np.zeros(5), 0.0, 0.0, 0.0]
        ub = np.r_[np.inf, np.full(5, RATIO_UPPER), 1e4, 10.0, 50.0 if self.use_inst else 1e-12]
        x0 = np.clip(x0, lb + 1e-9, ub - 1e-9)

        def resid(th):
            m = np.maximum(self._agg(self.model(th, s1, bkg)), 1e-12)
            return np.sign(d - m) * np.sqrt(np.maximum(_cstat_terms(d, m), 0.0))

        try:
            res = optimize.least_squares(resid, x0, bounds=(lb, ub), method="trf",
                                         x_scale="jac", max_nfev=400)
        except Exception as ex:  # noqa: BLE001
            out.update(quality=Q_NOCONV, success=False, error=str(ex))
            return out
        th = res.x
        m_full = np.maximum(self.model(th, s1, bkg), 1e-12)
        m = np.maximum(self._agg(m_full), 1e-12)
        C = float(np.sum(_cstat_terms(d, m)))
        # GoF over the main (1.15-7.5 keV) channels only
        Cm = float(np.sum(_cstat_terms(d[:self.n_main], m[:self.n_main])))
        eC, vC = cstat_expectation(m[:self.n_main])
        z = (Cm - eC.sum()) / np.sqrt(vC.sum())
        out.update(success=bool(res.success), nfev=int(res.nfev), cstat=C, cstat_expected=float(eC.sum()),
                   gof_z=float(z), cstat_main=Cm, dof=int(self.n_main + 1 - len(th)))
        for k, v in zip(PARAMS, th):
            out[k] = float(v)
        out["bkg_counts_model"] = float(th[7] * bkg[0][self.main].sum())
        out["inst_counts_model"] = float(th[8] * bkg[1][self.main].sum())
        out["src_counts_model"] = float(m_full[self.main].sum()) - out["bkg_counts_model"] - out["inst_counts_model"]

        # Fisher covariance (finite-difference model Jacobian)
        free = np.arange(th.size) if self.use_inst else np.arange(th.size - 1)
        J = np.empty((m.size, free.size))
        for jj, j in enumerate(free):
            h = 1e-4 * max(abs(th[j]), 1e-3)
            tp = th.copy()
            tp[j] += h
            J[:, jj] = (self._agg(self.model(tp, s1, bkg)) - m) / h
        F = J.T @ (J / m[:, None])
        sd = np.full(th.size, np.nan)
        try:
            cov = np.linalg.inv(F)
            sd[free] = np.sqrt(np.clip(np.diag(cov), 0, None))
            corr = cov / np.outer(sd[free], sd[free])
            cond = float(np.linalg.cond(corr))
        except np.linalg.LinAlgError:
            corr, cond = None, np.inf
            flag |= Q_COVFAIL
        for k, v in zip(PARAMS, sd):
            out[f"sd_{k}"] = float(v)
        out["cond"] = cond
        if corr is not None:
            out["corr"] = corr

        if not res.success:
            flag |= Q_NOCONV
        if d[:self.n_main].sum() < self.min_counts:
            flag |= Q_LOWCOUNTS
        if abs(z) > self.gof_z_max:
            flag |= Q_BADGOF
        at_bound = (th[1:6] <= 1e-6) | (th[1:6] >= RATIO_UPPER - 1e-6)
        for e, b in zip(FIT_ELEMENTS, at_bound):
            out[f"atbound_{e}"] = bool(b)
        if at_bound[0] or at_bound[1]:  # Mg or Al at a bound
            flag |= Q_ATBOUND
        if cond > self.cond_max:
            flag |= Q_ILLCOND
        out["quality"] = flag
        return out
