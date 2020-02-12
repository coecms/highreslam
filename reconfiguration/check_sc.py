#!/usr/bin/env python
# Copyright 2020 Scott Wales
# author: Scott Wales <scott.wales@unimelb.edu.au>
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import argparse
import xarray
import numpy

parser = argparse.ArgumentParser()
parser.add_argument('input', help="Input File")
args = parser.parse_args()

ds = xarray.open_dataset(args.input)
output = ds.output

if numpy.any(output < 1):
    print(output)
    raise Exception("Outputs below 1")

if numpy.any(output > ds.sizes['flat']):
    print(output)
    raise Exception("Outputs above grid size")
