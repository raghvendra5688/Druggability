#!/bin/bash
# Wrapper that sets the conda env's libstdc++ before running the Python script.
# Usage: nohup bash scripts/run_panviral_predictions.sh > logs/panviral_predictions.log 2>&1 &

EFFICHEM=/export/home/rmall/micromamba/envs/effichem
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

export LD_LIBRARY_PATH="${EFFICHEM}/lib:${LD_LIBRARY_PATH}"

exec "${EFFICHEM}/bin/python" "${SCRIPT_DIR}/panviral_drugpurpose_predictions.py"
