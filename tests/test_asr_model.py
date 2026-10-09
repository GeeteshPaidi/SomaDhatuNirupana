"""Unit tests for the ASR-rework forward model (TRACKER tasks 1.2, 1.3, 1.4)."""
import os
import sys

import numpy as np
import pytest
import xraylib
from scipy.integrate import quad

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src", "asr"))
import xrf_model as xm  # noqa: E402

# Typical mare-ish ratios (wt): Mg/Si 0.25, Al/Si 0.30, Ca/Si 0.40, Ti/Si 0.05, Fe/Si 0.60
RATIOS = {"Mg": 0.25, "Al": 0.30, "Ca": 0.40, "Ti": 0.05, "Fe": 0.60}


def test_composition_closure():
    comp = xm.composition_from_ratios(RATIOS)
    assert abs(sum(comp.values()) - 1) < 1e-12
    for e, r in RATIOS.items():
        assert comp[e] / comp["Si"] == pytest.approx(r)
    # O in lunar soils is ~42-45 wt%
    assert 0.40 < comp["O"] < 0.46


def test_line_constant_matches_xraylib():
    """Q * tau_K must equal xraylib's fluorescence-line cross-section (non-cascade)."""
    fm = xm.ForwardModel(secondary=False)
    for ln in fm.lines:
        E = 8.0
        ours = ln.Q * xraylib.CS_Photo_Partial(ln.Z, xraylib.K_SHELL, E)
        lid = xraylib.KA_LINE if ln.label == "Ka" else xraylib.KB_LINE
        ref = xraylib.CS_FluorLine(ln.Z, lid, E)
        assert ours == pytest.approx(ref, rel=0.05), ln.el


def test_primary_si_against_independent_quadrature():
    """Hand-built Sherman integral for Si Ka with scipy.quad vs the gridded model."""
    psi1 = 40.0
    fm = xm.ForwardModel(kT=1.07, secondary=False, de=0.002)
    comp = xm.composition_from_ratios(RATIOS)
    Zs = {e: xraylib.SymbolToAtomicNumber(e) for e in comp}
    si = [ln for ln in fm.lines if ln.el == "Si" and ln.label == "Ka"][0]
    mu = lambda E: sum(c * xraylib.CS_Total(Zs[e], E) for e, c in comp.items())  # noqa: E731
    mu_l = mu(si.E)
    s1 = np.sin(np.radians(psi1))
    w = xraylib.FluorYield(14, xraylib.K_SHELL) * xraylib.RadRate(14, xraylib.KA_LINE)
    f = lambda E: xm.solar_spectrum(E, 1.07) * xraylib.CS_Photo_Partial(14, xraylib.K_SHELL, E) / (mu(E) / s1 + mu_l)  # noqa: E731
    edge = xraylib.EdgeEnergy(14, xraylib.K_SHELL)
    ref = comp["Si"] * w * quad(f, edge, 20.0, limit=200, points=[2.47, 4.04, 4.97, 7.11])[0]
    idx = fm.lines.index(si)
    ours = fm.line_intensities(RATIOS, psi1)[idx]
    assert ours == pytest.approx(ref, rel=0.01)


def test_primary_scales_with_incidence():
    """Grazing incidence must reduce intensity; psi1 -> 90 deg maximal."""
    fm = xm.ForwardModel(secondary=False)
    i30 = fm.line_intensities(RATIOS, 30.0)
    i80 = fm.line_intensities(RATIOS, 80.0)
    i5 = fm.line_intensities(RATIOS, 5.0)
    assert np.all(i5 < i30) and np.all(i30 < i80)


def test_secondary_is_small_and_coupled():
    fm = xm.ForwardModel(secondary=True)
    prim, sec, _ = fm.line_intensities(RATIOS, 45.0, return_parts=True)
    frac = sec / prim
    names = [f"{ln.el}{ln.label}" for ln in fm.lines]
    d = dict(zip(names, frac))
    assert d["FeKa"] == 0  # nothing above Fe K edge among the lines
    assert 0 < d["MgKa"] < 0.3 and 0 < d["CaKa"] < 0.3
    # coupling: doubling Fe/Si must increase the Ca secondary term
    r2 = dict(RATIOS, Fe=1.2)
    _, sec2, _ = fm.line_intensities(r2, 45.0, return_parts=True)
    ca = names.index("CaKa")
    assert sec2[ca] > sec[ca]


@pytest.mark.parametrize("E_line", [1.2536, 1.7398, 3.69, 6.3995])
@pytest.mark.parametrize("sigma0", [0.03, 0.06])
def test_detector_line_area_conserved(E_line, sigma0):
    det = xm.Detector(gain=0.01363, offset=0.017, sigma0=sigma0)
    r = det.line_response(E_line)
    assert r.sum() == pytest.approx(1.0, abs=1e-9)
    assert np.all(r >= 0)
    # centroid of the binned response equals the line energy (to < 1% of a channel)
    assert np.sum(r * det.E_center) == pytest.approx(E_line, abs=0.01 * det.gain)
