"""
Synthetic LAS File Generator
------------------------------
Generates a realistic synthetic well log LAS file
simulating a sedimentary sequence with alternating
sand, shale, and carbonate intervals.

Curves generated:
  DEPTH — Measured depth (m)
  GR    — Gamma Ray (API)
  SP    — Spontaneous Potential (mV)
  ILD   — Deep Resistivity (ohm.m)
  RHOB  — Bulk Density (g/cc)
  NPHI  — Neutron Porosity (v/v)
  DT    — Sonic Transit Time (µs/ft)
  CALI  — Caliper (inches)
"""

import numpy as np
import lasio
import pandas as pd

np.random.seed(21)

# ── Depth ─────────────────────────────────────────────────────────────────────
depth_top    = 1500.0   # m
depth_bottom = 2100.0   # m
depth_step   = 0.5      # m
depth        = np.arange(depth_top, depth_bottom, depth_step)
n            = len(depth)

# ── Lithology model ───────────────────────────────────────────────────────────
# 0 = shale, 1 = sand, 2 = carbonate
litho = np.zeros(n, dtype=int)

# Define intervals (index based)
intervals = [
    (0,   80,  0),   # shale
    (80,  160, 1),   # sand
    (160, 220, 0),   # shale
    (220, 280, 2),   # carbonate
    (280, 360, 0),   # shale
    (360, 460, 1),   # sand (reservoir)
    (460, 520, 0),   # shale
    (520, 580, 2),   # carbonate
    (580, 660, 0),   # shale
    (660, 740, 1),   # sand
    (740, 800, 0),   # shale
    (800, 860, 2),   # carbonate
    (860, 940, 0),   # shale
    (940, 1000, 1),  # sand
    (1000, 1060, 0), # shale
    (1060, 1120, 2), # carbonate
                     # remaining = shale
]
for start, end, ltype in intervals:
    end = min(end, n)
    litho[start:end] = ltype

# ── Log responses per lithology ───────────────────────────────────────────────
# GR (API): shale=80-120, sand=15-45, carbonate=10-30
GR = np.where(litho == 0,
              np.random.uniform(80, 120, n),
     np.where(litho == 1,
              np.random.uniform(15, 45, n),
              np.random.uniform(10, 30, n)))

# SP (mV): shale≈0, sand negative deflection, carbonate slight negative
SP = np.where(litho == 0,
              np.random.uniform(-5, 5, n),
     np.where(litho == 1,
              np.random.uniform(-80, -40, n),
              np.random.uniform(-20, -5, n)))

# ILD (ohm.m): shale=1-5, sand=10-80, carbonate=20-150
ILD = np.where(litho == 0,
               np.random.uniform(1, 5, n),
      np.where(litho == 1,
               np.random.uniform(10, 80, n),
               np.random.uniform(20, 150, n)))

# RHOB (g/cc): shale=2.4-2.65, sand=2.1-2.45, carbonate=2.5-2.71
RHOB = np.where(litho == 0,
                np.random.uniform(2.40, 2.65, n),
       np.where(litho == 1,
                np.random.uniform(2.10, 2.45, n),
                np.random.uniform(2.50, 2.71, n)))

# NPHI (v/v): shale=0.25-0.45, sand=0.10-0.30, carbonate=0.05-0.20
NPHI = np.where(litho == 0,
                np.random.uniform(0.25, 0.45, n),
       np.where(litho == 1,
                np.random.uniform(0.10, 0.30, n),
                np.random.uniform(0.05, 0.20, n)))

# DT (µs/ft): shale=90-120, sand=55-90, carbonate=45-65
DT = np.where(litho == 0,
              np.random.uniform(90, 120, n),
     np.where(litho == 1,
              np.random.uniform(55, 90, n),
              np.random.uniform(45, 65, n)))

# CALI (inches): mostly 8.5", washes out in shale
CALI = np.where(litho == 0,
                np.random.uniform(8.8, 10.5, n),
                np.random.uniform(8.4, 8.7, n))

# Add noise
for arr in [GR, SP, ILD, RHOB, NPHI, DT, CALI]:
    arr += np.random.normal(0, arr.std() * 0.03, n)

GR   = np.clip(np.round(GR,   2),  0,   200)
SP   = np.round(SP,   2)
ILD  = np.clip(np.round(ILD,  3),  0.1, 500)
RHOB = np.clip(np.round(RHOB, 4),  1.8, 2.9)
NPHI = np.clip(np.round(NPHI, 4),  0,   0.6)
DT   = np.clip(np.round(DT,   2),  40,  140)
CALI = np.clip(np.round(CALI, 2),  6,   16)

# ── Write LAS ─────────────────────────────────────────────────────────────────
las = lasio.LASFile()

las.well["WELL"] = lasio.HeaderItem("WELL", value="SYNTHETIC-01",
                                     descr="Synthetic Well")
las.well["FLD"]  = lasio.HeaderItem("FLD",  value="DEMO FIELD")
las.well["LOC"]  = lasio.HeaderItem("LOC",  value="SYNTHETIC LOCATION")
las.well["COMP"] = lasio.HeaderItem("COMP", value="DEMO")
las.well["DATE"] = lasio.HeaderItem("DATE", value="2026")
las.well["STRT"] = lasio.HeaderItem("STRT", value=depth_top,
                                     unit="M")
las.well["STOP"] = lasio.HeaderItem("STOP", value=depth[-1],
                                     unit="M")
las.well["STEP"] = lasio.HeaderItem("STEP", value=depth_step,
                                     unit="M")
las.well["NULL"] = lasio.HeaderItem("NULL", value=-999.25)

las.append_curve("DEPTH", depth, unit="M",   descr="Measured Depth")
las.append_curve("GR",    GR,    unit="API", descr="Gamma Ray")
las.append_curve("SP",    SP,    unit="MV",  descr="Spontaneous Potential")
las.append_curve("ILD",   ILD,   unit="OHMM",descr="Deep Resistivity")
las.append_curve("RHOB",  RHOB,  unit="G/CC",descr="Bulk Density")
las.append_curve("NPHI",  NPHI,  unit="V/V", descr="Neutron Porosity")
las.append_curve("DT",    DT,    unit="US/F",descr="Sonic Transit Time")
las.append_curve("CALI",  CALI,  unit="IN",  descr="Caliper")

las.write("sample.las")
print(f"LAS file written: {len(depth)} samples, "
      f"{depth_top}–{depth[-1]} m")
print(f"Lithology: {np.sum(litho==0)} shale | "
      f"{np.sum(litho==1)} sand | "
      f"{np.sum(litho==2)} carbonate samples")
