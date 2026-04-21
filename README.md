# petrophysics-synthetic-seismic-python

![Python](https://img.shields.io/badge/Python-3.8+-blue)
![License](https://img.shields.io/badge/License-MIT-green)
![Status](https://img.shields.io/badge/Status-Complete-brightgreen)

A Python workflow for subsurface characterisation using well log petrophysical evaluation and seismic forward modelling. Designed for geoscientists working with LAS-format borehole data who need to derive rock physical properties and understand their seismic expression.

---

## Overview

Borehole well logs record the physical properties of subsurface formations at centimetre-scale resolution. This project processes those logs to extract meaningful petrophysical parameters — lithology, porosity, shale volume, water saturation — and then uses the derived acoustic properties to generate a synthetic seismogram, bridging the gap between borehole measurements and surface seismic data.

Developed as part of applied geophysics training in subsurface data integration.

---

## What This Does

**Well Log Analysis**
- Reads standard LAS format borehole logs (GR, SP, ILD, RHOB, NPHI, DT, CALI)
- Flags lithology using GR cutoff thresholds
- Computes density porosity, average porosity, shale volume, and water saturation
- Generates multi-track log display and neutron-density crossplot

**Synthetic Seismogram**
- Converts sonic log (DT) to P-wave velocity
- Computes acoustic impedance (AI) from velocity and bulk density
- Derives reflection coefficient series at layer boundaries
- Generates a Ricker wavelet using the bruges geophysics library
- Convolves reflection coefficients with wavelet to produce synthetic seismic trace
- Displays wavelet in time and frequency domain

---

## Why Synthetic Seismograms Matter

A synthetic seismogram directly connects well log data to seismic reflection data. It allows geoscientists to:
- Identify which geological boundaries produce seismic reflections
- Calibrate surface seismic data to borehole measurements
- Understand the seismic expression of specific lithologies
- Validate seismic interpretation at well locations

---

## Petrophysical Equations

**Density Porosity:**
$$\phi_D = \frac{\rho_{ma} - \rho_b}{\rho_{ma} - \rho_{fl}}$$

**Vshale (Linear GR method):**
$$V_{sh} = \frac{GR - GR_{clean}}{GR_{shale} - GR_{clean}}$$

**Water Saturation (Archie's Equation):**
$$S_w = \left(\frac{a \cdot R_w}{\phi^m \cdot R_t}\right)^{1/n}$$

**Acoustic Impedance:**
$$AI = V_p \times \rho_b$$

**Reflection Coefficient:**
$$RC_i = \frac{AI_{i+1} - AI_i}{AI_{i+1} + AI_i}$$

---

## Well Log Dataset

Synthetic LAS file simulating a 600 m sedimentary section (1500–2100 m depth) with three lithological units:

| Curve | Unit | Description |
|---|---|---|
| GR | API | Gamma Ray |
| SP | mV | Spontaneous Potential |
| ILD | Ω·m | Deep Resistivity |
| RHOB | g/cc | Bulk Density |
| NPHI | v/v | Neutron Porosity |
| DT | µs/ft | Sonic Transit Time |
| CALI | inches | Caliper |

---

## Project Structure

```
petrophysics-synthetic-seismic-python/
│
├── generate_las.py          # Synthetic LAS generator
├── sample.las               # Well log input file
├── well_log_seismic.py      # Main petrophysics + seismic script
└── README.md
```

---

## Dependencies

| Library | Purpose |
|---|---|
| `lasio` | Industry-standard LAS file reading |
| `numpy` | Numerical computation, convolution, FFT |
| `pandas` | Log curve data handling |
| `matplotlib` | Well log visualisation, seismic display |
| `bruges` | Geophysics-specific wavelet generation |
| `scipy` | Signal processing |

```bash
pip install lasio numpy pandas matplotlib bruges scipy
```

---

## Usage

```bash
# Generate LAS file
python generate_las.py

# Run petrophysical analysis and synthetic seismogram
python well_log_seismic.py
```

---

## Key Outputs

**Multi-Track Log Display** — Eight-track visualisation showing lithology panel, GR, SP, ILD, RHOB, NPHI, DT, and computed porosity in a standard well log format.

**Neutron-Density Crossplot** — Classic lithology identification plot showing sand, shale, and carbonate clusters using NPHI vs RHOB.

**Synthetic Seismogram** — Five-panel figure showing GR, acoustic impedance, reflection coefficients, synthetic seismic trace (wiggle + fill), and Ricker wavelet side by side.

**Wavelet Spectrum** — Time-domain wavelet shape and frequency-domain amplitude spectrum confirming dominant frequency at 40 Hz.

---

## Research Context

This workflow is directly applicable to subsurface characterisation in sedimentary basins — including the Bengal Basin where fluvial and deltaic stratigraphy controls aquifer architecture and hydrocarbon prospectivity. Well log analysis and seismic forward modelling are foundational tools for both groundwater and petroleum geoscience investigations.

---

## Author

**Banani Jana**
M.Sc. Applied Geology (2nd Year), Presidency University, Kolkata
Email: bananijana2002@gmail.com

---

## License

MIT License
