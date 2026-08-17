"""
Figure 1 (paper Sec 2): Hybrid smoothing comparison, dark theme.

Shows the two-stage denoiser on one representative background-subtracted
dayside spectrum:
  - Original noisy (background-subtracted count rate)
  - Wavelet only (Stationary Wavelet Transform, sym8, level 2, soft threshold)
  - Wavelet + Savitzky-Golay (window 7, order 3)  -- the polished result fed
    to the FP fit.

Matches the clean dark-background style the user preferred. Reads a saved
example spectrum (raw/bgsub) from an .npz produced by
src/fitting/prepare_data.py (its --example_out), so the traces are the real
pipeline output, not a mock-up.

Usage:
    python fig1_preprocessing.py [example_npz] [channel_max]
Defaults: ../../outputs/figures/preprocessing_example.npz, 200
"""
import sys

import numpy as np
import pywt
from scipy.signal import savgol_filter

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def wavelet_denoise(rate, wavelet="sym8", level=2):
    """SWT soft-threshold denoise WITHOUT the Savitzky-Golay stage."""
    n = len(rate)
    coeffs = pywt.swt(rate, wavelet, level=level)
    out = []
    for cA, cD in coeffs:
        sigma = np.median(np.abs(cD - np.median(cD))) / 0.6745 if len(cD) else 0.0
        thresh = sigma * np.sqrt(2 * np.log(max(n, 2)))
        cD_soft = np.sign(cD) * np.maximum(np.abs(cD) - thresh, 0.0)
        out.append((cA, cD_soft))
    rec = pywt.iswt(out, wavelet)[:n]
    return rec


def main():
    npz = sys.argv[1] if len(sys.argv) > 1 else "../../outputs/figures/preprocessing_example.npz"
    ch_max = int(sys.argv[2]) if len(sys.argv) > 2 else 200

    d = np.load(npz)
    noisy = d["bgsub"].astype(float)           # background-subtracted, still noisy
    wav = np.clip(wavelet_denoise(noisy), 0.0, None)
    wav_sg = np.clip(savgol_filter(wav, window_length=7, polyorder=3), 0.0, None)

    ch = np.arange(len(noisy))
    m = ch < ch_max

    # Standard publication style: white background, black axes, restrained
    # colours, no in-axes title (the LaTeX \caption carries the description).
    # Figure is sized to the actual IEEEtran single-column print width
    # (3.5 in) so font sizes set here ARE the printed point sizes -- no
    # shrink-to-fit blur when \includegraphics[width=\linewidth] scales it.
    plt.rcParams.update({
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
        "font.size": 9,
        "axes.linewidth": 0.6,
    })

    fig, ax = plt.subplots(figsize=(3.5, 2.55), dpi=600)
    ax.plot(ch[m], noisy[m], color="0.6", linewidth=0.6, alpha=0.9,
            label="Original (background-subtracted)")
    ax.plot(ch[m], wav[m], color="#1f77b4", linewidth=0.8,
            label="Wavelet denoised")
    ax.plot(ch[m], wav_sg[m], color="#d62728", linewidth=1.1, linestyle="--",
            label="Wavelet + Savitzky--Golay")

    ax.set_xlim(0, ch_max)
    ax.set_ylim(bottom=0)
    ax.set_xlabel("Channel", fontsize=8.5)
    ax.set_ylabel("Count rate (counts s$^{-1}$)", fontsize=8.5)
    ax.tick_params(labelsize=8)
    ax.grid(True, color="0.85", linewidth=0.4)
    ax.set_axisbelow(True)
    ax.legend(frameon=True, framealpha=0.95, edgecolor="0.7", fontsize=6.5,
              loc="upper right", handlelength=1.8)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)

    fig.tight_layout()
    import os
    os.makedirs("../../outputs/figures/final", exist_ok=True)
    fig.savefig("../../outputs/figures/final/fig1_preprocessing.pdf", bbox_inches="tight")
    fig.savefig("../../outputs/figures/final/fig1_preprocessing.png", dpi=600, bbox_inches="tight")
    print("Saved fig1_preprocessing.pdf/.png to outputs/figures/final/")
    print(f"noisy sum={noisy.sum():.2f}  wav sum={wav.sum():.2f}  wav+sg sum={wav_sg.sum():.2f}")


if __name__ == "__main__":
    main()
