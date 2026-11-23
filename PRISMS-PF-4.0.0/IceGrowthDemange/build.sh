#!/bin/bash
# Convenience build wrapper for PRISMS-PF 4.0.0 on ACES.
# Set PRISMS_PF_INSTALL and DEAL_II_INSTALL to existing built/install locations.
# Use a compatible PRISMS-PF 4.0.0 and deal.II installation.
set -euo pipefail
module purge
module load foss/2023b CMake/3.27.6
: "${PRISMS_PF_INSTALL:?Set PRISMS_PF_INSTALL to your PRISMS-PF 4.0.0 install directory}"
: "${DEAL_II_INSTALL:?Set DEAL_II_INSTALL to your deal.II install directory}"
APP="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
cmake -S "$APP" -B "$APP/build" -DCMAKE_BUILD_TYPE=Release \
  -DPRISMS_PF_DIR="$PRISMS_PF_INSTALL/lib64/cmake/PRISMS-PF" \
  -Ddeal.II_DIR="$DEAL_II_INSTALL/lib/cmake/deal.II"
cmake --build "$APP/build" --parallel 4
