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

SCRIPTDIR=$(dirname "$(readlink -f "$0")")

pushd "$BUILD_DIR"

cp "$SCRIPT_DIR/write_spiral_circle_mod.f90" extract/um/src/utility/qxreconf
cp "$SCRIPT_DIR/gather_saw.f90" extract/um/src/utility/qxreconf
patch -p0 < "$SCRIPT_DIR/write_spiral_circle.patch"

popd
