# Scanning Millimetre-Wave Radar Foliage Penetration

This repository contains the data processing and analysis code for my undergraduate thesis at the University of Sydney, which experimentally evaluates a scanning 120 GHz FMCW radar system for foliage penetration sensing in outdoor vegetated environments, benchmarked directly against a Livox Mid-70 LiDAR.

## Overview

Four experiments were conducted across March and April 2026 at the University of Sydney:

| # | Experiment | Description |
|---|-----------|-------------|
| 1 | **Control Scan** | Free-space corner reflector baseline |
| 2 | **Plane Tree Canopy Scan** | In-canopy attenuation and volumetric reconstruction |
| 3 | **Reed Bed Scan** | Corner reflector detection through dense reed canopy |
| 4 | **Cherry Ballart Shrub Scan** | Detection through two faces of differing canopy density |

## Hardware

| Device | Specification |
|--------|--------------|
| Silicon Radar TRA_120_002 | 120 GHz FMCW radar transceiver |
| Livox Mid-70 | Solid-state LiDAR (905 nm) |
| Emlid Reach RS2 | Dual RTK GNSS receivers for spatial registration |

## Repository Structure

Code is organised by experiment. Each folder contains the analysis scripts for that experiment.

```
.
├── Control Scan/
├── New Plane tree/
├── Reed Bed/
├── Cherry Ballart/
└── Figures/
```

## Getting Started

### Prerequisites

List any dependencies needed to run the code (e.g. MATLAB version, Python packages).

### Usage

Open the relevant experiment folder and run the main analysis script. Processed scan data and figures will be generated in the respective output directories.

## Author

Raquel Kampel  
Bachelor of Engineering (Honours) — University of Sydney, 2026  

## License

This repository is for academic purposes. Please contact the author before reusing any code or data.
