"""
Well Log Analysis and Synthetic Seismogram Generation
-------------------------------------------------------
Module 1 — Well Log Petrophysical Analysis
  - Log loading and QC
  - Multi-track well log plot
  - Lithology flag from GR cutoff
  - Density porosity (PHID)
  - Neutron-Density crossplot
  - Water saturation (Archie's equation)

Module 2 — Synthetic Seismogram
  - Acoustic impedance from DT + RHOB
  - Reflection coefficient series
  - Ricker wavelet (bruges)
  - Convolution → synthetic trace
  - Wavelet amplitude spectrum (FFT)
  - Side-by-side log + synthetic display
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from matplotlib.colors import ListedColormap
import lasio
import bruges
from scipy.signal import welch
import warnings
import os

warnings.filterwarnings("ignore")

OUTPUT_DIR = "outputs"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ═══════════════════════════════════════════════════════════════════════════════
# MODULE 1 — WELL LOG ANALYSIS
# ═══════════════════════════════════════════════════════════════════════════════

# ── Load LAS ──────────────────────────────────────────────────────────────────
las  = lasio.read("sample.las")
df   = las.df().reset_index()
df.columns = [c.upper() for c in df.columns]
df   = df.rename(columns={"INDEX": "DEPTH"})
df   = df.replace(-999.25, np.nan)

print(f"Well: {las.well['WELL'].value}")
print(f"Depth range: {df['DEPTH'].min():.1f} – {df['DEPTH'].max():.1f} m")
print(f"Curves: {list(df.columns)}\n")

# ── Lithology flag from GR cutoff ─────────────────────────────────────────────
GR_SHALE_CUTOFF = 65   # API
GR_CARB_CUTOFF  = 25   # API — below this = carbonate

def litho_flag(gr):
    if gr >= GR_SHALE_CUTOFF: return 0   # shale
    elif gr <= GR_CARB_CUTOFF: return 2  # carbonate
    else: return 1                        # sand

df["LITHO"] = df["GR"].apply(litho_flag)
litho_colors = {0: "#A0785A", 1: "#F5E642", 2: "#A8D5E2"}
litho_labels = {0: "Shale", 1: "Sand", 2: "Carbonate"}

# ── Petrophysical Calculations ────────────────────────────────────────────────
# Matrix and fluid densities
RHO_MA  = 2.65   # g/cc quartz matrix
RHO_FL  = 1.00   # g/cc fresh water
RHO_SH  = 2.68   # g/cc shale matrix

# Density Porosity
df["PHID"] = (RHO_MA - df["RHOB"]) / (RHO_MA - RHO_FL)
df["PHID"] = df["PHID"].clip(0, 0.45)

# Average porosity (neutron-density)
df["PHIA"] = ((df["PHID"] + df["NPHI"]) / 2).clip(0, 0.45)

# Vshale from GR (linear)
GR_CLEAN = df["GR"].quantile(0.05)
GR_SHALE = df["GR"].quantile(0.95)
df["VSH"] = ((df["GR"] - GR_CLEAN) / (GR_SHALE - GR_CLEAN)).clip(0, 1)

# Water saturation — Archie's equation: Sw = (a / (phi^m × Rt)) ^ (1/n)
a, m, n = 1.0, 2.0, 2.0
RW = 0.05   # ohm.m formation water resistivity
df["SW"] = ((a * RW) / (df["PHIA"] ** m * df["ILD"])) ** (1/n)
df["SW"] = df["SW"].clip(0, 1)
df["SHC"] = 1 - df["SW"]   # hydrocarbon saturation

print("── Petrophysical Summary ──")
print(f"  Mean PHID : {df['PHID'].mean():.3f}")
print(f"  Mean PHIA : {df['PHIA'].mean():.3f}")
print(f"  Mean VSH  : {df['VSH'].mean():.3f}")
print(f"  Mean SW   : {df['SW'].mean():.3f}")
print()

# ── Plot 1: Multi-Track Well Log ──────────────────────────────────────────────
fig = plt.figure(figsize=(18, 14))
gs  = gridspec.GridSpec(1, 8, wspace=0.08)

depth = df["DEPTH"].values

tracks = [
    ("LITHO",  "Lithology",   None,        None),
    ("GR",     "GR (API)",    0,    150),
    ("SP",     "SP (mV)",     -120, 20),
    ("ILD",    "ILD (Ω·m)",   0.1,  200),
    ("RHOB",   "RHOB (g/cc)", 1.8,  2.9),
    ("NPHI",   "NPHI (v/v)",  0.45, -0.05),
    ("DT",     "DT (µs/ft)",  140,  40),
    ("PHIA",   "PHIA (v/v)",  0,    0.45),
]

for i, (curve, label, xmin, xmax) in enumerate(tracks):
    ax = fig.add_subplot(gs[i])

    if curve == "LITHO":
        for j in range(len(depth) - 1):
            fc = litho_colors[df["LITHO"].iloc[j]]
            ax.fill_betweenx([depth[j], depth[j+1]], 0, 1, color=fc)
        ax.set_xlim(0, 1)
        ax.set_xticks([])
        patches = [mpatches.Patch(color=litho_colors[k],
                                  label=litho_labels[k])
                   for k in litho_colors]
        ax.legend(handles=patches, loc="lower right",
                  fontsize=6, framealpha=0.7)
    else:
        color = {"GR":   "#2E7D32", "SP":  "#1565C0",
                 "ILD":  "#B71C1C", "RHOB":"#6A1B9A",
                 "NPHI": "#00838F", "DT":  "#E65100",
                 "PHIA": "#F9A825"}.get(curve, "black")
        ax.plot(df[curve], depth, color=color, linewidth=0.8)
        ax.set_xlim(xmin, xmax)
        ax.fill_betweenx(depth, xmin, df[curve],
                         alpha=0.15, color=color)

    ax.set_ylim(depth.max(), depth.min())
    ax.set_title(label, fontsize=8, fontweight="bold", pad=4)
    ax.tick_params(labelsize=7)
    ax.spines["right"].set_visible(False)
    if i > 0:
        ax.set_yticklabels([])
    else:
        ax.set_ylabel("Depth (m)", fontsize=9)

fig.suptitle("SYNTHETIC-01 — Multi-Track Well Log Display",
             fontsize=14, fontweight="bold", y=1.01)
plt.tight_layout()
fig.savefig(f"{OUTPUT_DIR}/well_log_tracks.png",
            dpi=180, bbox_inches="tight")
plt.close()
print("Saved: well_log_tracks.png")

# ── Plot 2: Neutron-Density Crossplot ─────────────────────────────────────────
fig, ax = plt.subplots(figsize=(7, 6))
lcolors = [litho_colors[l] for l in df["LITHO"]]
sc = ax.scatter(df["NPHI"], df["RHOB"], c=lcolors,
                s=18, alpha=0.7, edgecolors="none")
ax.set_xlabel("NPHI (v/v)", fontsize=11)
ax.set_ylabel("RHOB (g/cc)", fontsize=11)
ax.invert_yaxis()
ax.set_title("Neutron–Density Crossplot\n(Lithology Identification)",
             fontsize=12, fontweight="bold")
patches = [mpatches.Patch(color=litho_colors[k],
                          label=litho_labels[k])
           for k in litho_colors]
ax.legend(handles=patches, fontsize=9)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.tight_layout()
fig.savefig(f"{OUTPUT_DIR}/neutron_density_crossplot.png", dpi=180)
plt.close()
print("Saved: neutron_density_crossplot.png")

# ── Plot 3: Porosity and Saturation Summary ───────────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(13, 6))

# PHIA histogram
axes[0].hist(df["PHIA"].dropna(), bins=30,
             color="#F9A825", edgecolor="white", alpha=0.9)
axes[0].axvline(df["PHIA"].mean(), color="red", linestyle="--",
                linewidth=1.5, label=f"Mean={df['PHIA'].mean():.3f}")
axes[0].set_xlabel("Porosity (v/v)", fontsize=10)
axes[0].set_ylabel("Frequency", fontsize=10)
axes[0].set_title("Average Porosity (PHIA)", fontsize=11,
                  fontweight="bold")
axes[0].legend(fontsize=9)

# SW histogram
axes[1].hist(df["SW"].dropna(), bins=30,
             color="#1565C0", edgecolor="white", alpha=0.9)
axes[1].axvline(df["SW"].mean(), color="red", linestyle="--",
                linewidth=1.5, label=f"Mean={df['SW'].mean():.3f}")
axes[1].set_xlabel("Water Saturation (Sw)", fontsize=10)
axes[1].set_ylabel("Frequency", fontsize=10)
axes[1].set_title("Water Saturation (Archie)", fontsize=11,
                  fontweight="bold")
axes[1].legend(fontsize=9)

# VSH histogram
axes[2].hist(df["VSH"].dropna(), bins=30,
             color="#A0785A", edgecolor="white", alpha=0.9)
axes[2].axvline(df["VSH"].mean(), color="red", linestyle="--",
                linewidth=1.5, label=f"Mean={df['VSH'].mean():.3f}")
axes[2].set_xlabel("Vshale", fontsize=10)
axes[2].set_ylabel("Frequency", fontsize=10)
axes[2].set_title("Shale Volume (Vsh)", fontsize=11,
                  fontweight="bold")
axes[2].legend(fontsize=9)

for ax in axes:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

fig.suptitle("Petrophysical Parameter Distributions",
             fontsize=13, fontweight="bold")
plt.tight_layout()
fig.savefig(f"{OUTPUT_DIR}/petrophysical_distributions.png", dpi=180)
plt.close()
print("Saved: petrophysical_distributions.png")

# ═══════════════════════════════════════════════════════════════════════════════
# MODULE 2 — SYNTHETIC SEISMOGRAM
# ═══════════════════════════════════════════════════════════════════════════════

print("\n── Synthetic Seismogram ──")

# ── Acoustic Impedance ────────────────────────────────────────────────────────
# Convert DT (µs/ft) to velocity (m/s)
df["VP"]  = (1e6 / df["DT"]) * 0.3048   # µs/ft → m/s
df["AI"]  = df["VP"] * df["RHOB"]        # acoustic impedance

# ── Reflection Coefficients ───────────────────────────────────────────────────
ai     = df["AI"].ffill().bfill().values
rc     = np.zeros(len(ai))
rc[1:] = (ai[1:] - ai[:-1]) / (ai[1:] + ai[:-1])

# ── Ricker Wavelet (bruges) ───────────────────────────────────────────────────
dt_s   = 0.001            # 1 ms sample interval
f_dom  = 40               # 40 Hz dominant frequency
_wav   = bruges.filters.ricker(duration=0.128, dt=dt_s, f=f_dom)
wav    = np.array(_wav.amplitude).flatten()
t_wav  = np.array(_wav.time).flatten() * 1000   # ms

print(f"  Wavelet: Ricker {f_dom} Hz, {len(wav)} samples")
print(f"  AI range: {ai.min():.0f} – {ai.max():.0f} (g/cc·m/s)")
print(f"  Non-zero RC: {np.sum(np.abs(rc) > 0.01)}")

# ── Convolve ──────────────────────────────────────────────────────────────────
synthetic = np.convolve(rc, wav, mode="same")
synthetic = synthetic / np.max(np.abs(synthetic))   # normalise

# ── Time axis (TWT) ───────────────────────────────────────────────────────────
depth_step_m = float(depth[1] - depth[0])
twt = np.cumsum(np.ones(len(depth)) * depth_step_m * 2 / df["VP"].mean()) * 1000  # ms

# ── Plot 4: Synthetic Seismogram ──────────────────────────────────────────────
fig = plt.figure(figsize=(16, 13))
gs  = gridspec.GridSpec(1, 5, wspace=0.08)

# Track 1 — GR
ax1 = fig.add_subplot(gs[0])
ax1.plot(df["GR"], depth, color="#2E7D32", linewidth=0.9)
ax1.fill_betweenx(depth, 0, df["GR"], alpha=0.2, color="#2E7D32")
ax1.set_xlim(0, 150)
ax1.set_ylim(depth.max(), depth.min())
ax1.set_xlabel("GR (API)", fontsize=9)
ax1.set_ylabel("Depth (m)", fontsize=10)
ax1.set_title("GR", fontsize=10, fontweight="bold")
ax1.tick_params(labelsize=8)
ax1.spines["right"].set_visible(False)

# Track 2 — Acoustic Impedance
ax2 = fig.add_subplot(gs[1])
ax2.plot(df["AI"], depth, color="#6A1B9A", linewidth=0.9)
ax2.fill_betweenx(depth, df["AI"].min(), df["AI"],
                  alpha=0.2, color="#6A1B9A")
ax2.set_ylim(depth.max(), depth.min())
ax2.set_xlabel("AI (g/cc·m/s)", fontsize=9)
ax2.set_title("Acoustic\nImpedance", fontsize=10, fontweight="bold")
ax2.set_yticklabels([])
ax2.tick_params(labelsize=8)
ax2.spines["right"].set_visible(False)

# Track 3 — Reflection Coefficients
ax3 = fig.add_subplot(gs[2])
ax3.plot(rc, depth, color="#B71C1C", linewidth=0.6, alpha=0.8)
ax3.axvline(0, color="black", linewidth=0.5, alpha=0.5)
ax3.set_xlim(-0.25, 0.25)
ax3.set_ylim(depth.max(), depth.min())
ax3.set_xlabel("RC", fontsize=9)
ax3.set_title("Reflection\nCoefficients", fontsize=10, fontweight="bold")
ax3.set_yticklabels([])
ax3.tick_params(labelsize=8)
ax3.spines["right"].set_visible(False)

# Track 4 — Synthetic Seismogram (wiggle + fill)
ax4 = fig.add_subplot(gs[3])
ax4.plot(synthetic, depth, color="#01579B", linewidth=0.8)
ax4.fill_betweenx(depth, 0, synthetic,
                  where=(synthetic > 0), color="#01579B", alpha=0.5)
ax4.fill_betweenx(depth, 0, synthetic,
                  where=(synthetic < 0), color="#E53935", alpha=0.4)
ax4.axvline(0, color="black", linewidth=0.5, alpha=0.5)
ax4.set_xlim(-1.5, 1.5)
ax4.set_ylim(depth.max(), depth.min())
ax4.set_xlabel("Amplitude", fontsize=9)
ax4.set_title("Synthetic\nSeismogram", fontsize=10, fontweight="bold")
ax4.set_yticklabels([])
ax4.tick_params(labelsize=8)
ax4.spines["right"].set_visible(False)

# Track 5 — Ricker Wavelet
ax5 = fig.add_subplot(gs[4])
ax5.plot(wav, np.linspace(depth.min(), depth.max(), len(wav)),
         color="#E65100", linewidth=1.5)
ax5.fill_betweenx(np.linspace(depth.min(), depth.max(), len(wav)),
                  0, wav, alpha=0.3, color="#E65100")
ax5.axvline(0, color="black", linewidth=0.5, alpha=0.5)
ax5.set_xlabel("Amplitude", fontsize=9)
ax5.set_title(f"Ricker Wavelet\n({f_dom} Hz)", fontsize=10,
              fontweight="bold")
ax5.set_yticklabels([])
ax5.tick_params(labelsize=8)
ax5.spines["right"].set_visible(False)

fig.suptitle("SYNTHETIC-01 — Synthetic Seismogram Generation\n"
             "(GR → Acoustic Impedance → RC → Convolution with Ricker Wavelet)",
             fontsize=13, fontweight="bold", y=1.02)
plt.tight_layout()
fig.savefig(f"{OUTPUT_DIR}/synthetic_seismogram.png",
            dpi=180, bbox_inches="tight")
plt.close()
print("Saved: synthetic_seismogram.png")

# ── Plot 5: Wavelet Amplitude Spectrum ────────────────────────────────────────
N      = len(wav)
freqs  = np.fft.rfftfreq(N, d=dt_s)
amp    = np.abs(np.fft.rfft(wav))
amp   /= amp.max()

fig, axes = plt.subplots(1, 2, figsize=(12, 4))

axes[0].plot(t_wav - t_wav.mean(), wav, color="#E65100", linewidth=2)
axes[0].fill_between(t_wav - t_wav.mean(), 0, wav,
                     alpha=0.3, color="#E65100")
axes[0].axhline(0, color="black", linewidth=0.5)
axes[0].set_xlabel("Time (ms)", fontsize=10)
axes[0].set_ylabel("Amplitude", fontsize=10)
axes[0].set_title(f"Ricker Wavelet — Time Domain ({f_dom} Hz)",
                  fontsize=11, fontweight="bold")
axes[0].spines["top"].set_visible(False)
axes[0].spines["right"].set_visible(False)

axes[1].plot(freqs, amp, color="#1565C0", linewidth=2)
axes[1].fill_between(freqs, 0, amp, alpha=0.3, color="#1565C0")
axes[1].axvline(f_dom, color="red", linestyle="--",
                linewidth=1.5, label=f"Dominant freq = {f_dom} Hz")
axes[1].set_xlabel("Frequency (Hz)", fontsize=10)
axes[1].set_ylabel("Normalised Amplitude", fontsize=10)
axes[1].set_title("Wavelet Amplitude Spectrum — Frequency Domain",
                  fontsize=11, fontweight="bold")
axes[1].legend(fontsize=9)
axes[1].set_xlim(0, 150)
axes[1].spines["top"].set_visible(False)
axes[1].spines["right"].set_visible(False)

plt.tight_layout()
fig.savefig(f"{OUTPUT_DIR}/wavelet_spectrum.png", dpi=180)
plt.close()
print("Saved: wavelet_spectrum.png")

print(f"\n✓ All outputs saved to: {OUTPUT_DIR}")
print(f"  Module 1 — Well Log: 3 plots")
print(f"  Module 2 — Synthetic Seismic: 2 plots")
