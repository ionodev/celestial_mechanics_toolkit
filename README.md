# Celestial Mechanics Toolkit

Orbital mechanics scripts from the FYS-3000 (Space Mission Design) course
project: atmospheric drag, TLE-based satellite tracking, and interplanetary
transfer trajectory design.

(Constellation cost/link budget scripts live in
[GIROS](../GIROS); GPS signal reception lives in [gps_sdr](../gps_sdr).)

## Layout

| Folder | Description |
|---|---|
| [`celestial_mechanics/`](celestial_mechanics) | Atmospheric drag on a cubesat, drag-driven decay height from a curve fit, and TLE-based satellite pass tracking. |
| [`lambert_problem/`](lambert_problem) | Lambert-problem porkchop plots for an Earth-to-comet transfer (log-scale and linear-scale Δv color mapping). |

## Data requirements

Small reference data (`MSIS.dat`, `TLE.txt`) is included directly in
`celestial_mechanics/`.

## Requirements

Scripts use `numpy`, `scipy`, and `matplotlib`; `pykep` and `tqdm` (Lambert
problem); `sgp4` and `astropy` (TLE tracking).
