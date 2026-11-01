#!/usr/bin/env bash
# Reproduces every numerical result in Sections 6 and 7 of the report.
# Scripts must run in this order: the figure scripts read what the solvers write.
set -e
cd "$(dirname "$0")"

echo "=============================================================="
echo " 1/6  Core thermal model (both interpass policies)   ~2-4 min"
echo "=============================================================="
python3 src/01_thermal_model.py

echo "=============================================================="
echo " 2/6  Mesh and time-step convergence                 ~3 min"
echo "=============================================================="
python3 src/02_convergence_mesh_timestep.py

echo "=============================================================="
echo " 3/6  Differential convergence test                  ~10-20 min"
echo "=============================================================="
python3 src/03_convergence_differential.py

echo "=============================================================="
echo " 4/6  Cost model sensitivity                         ~1 s"
echo "=============================================================="
python3 src/04_cost_sensitivity.py

echo "=============================================================="
echo " 5/6  Thermal figures"
echo "=============================================================="
python3 src/05_figures_thermal.py

echo "=============================================================="
echo " 6/6  Convergence and tornado figures"
echo "=============================================================="
python3 src/06_figures_convergence_tornado.py

echo
echo "Done. Data in results/data/, figures in results/figures/."
