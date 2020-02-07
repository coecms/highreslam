#!/bin/bash
#  Copyright 2020 Scott Wales
#
#  \author  Scott Wales <scott.wales@unimelb.edu.au>
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.

set -eu

RUNID="$1"
CYLC_DIR="~/cylc-run/$RUNID"
BUILD_DIR="$CYLC_DIR/share/fcm_make"

test -d "$BUILD_DIR"

SCRIPT_DIR=$(dirname "$(readlink -f "$0")")

module load intel-compiler
module load openmpi/4.0.2
module load netcdf/4.7.1

mkdir -p build
pushd build

mpif90 -O3 -g -traceback -c "$BUILD_DIR/extract/shumlib/shum_constants/src/f_shum_conversions_mod.f90"
mpif90 -O3 -g -traceback -c "$BUILD_DIR/extract/shumlib/shum_spiral_search/src/f_shum_spiral_search.f90"
mpif90 -O3 -g -traceback -c "$SCRIPT_DIR/offline_spiral_circle.f90"
mpif90 "$SCRIPT_DIR/main.f90" f_shum_conversions_mod.o f_shum_spiral_search.o offline_spiral_circle.o -o offline_spiral_circle -lnetcdff

popd

qsub -v PROJECT,RUNID "$SCRIPT_DIR/run_offline_spiral_circle.pbs"
