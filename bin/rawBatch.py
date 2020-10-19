#!/usr/bin/python3
#
# Run pre-processing step on a batch of fastq files, then output results to a single JSON file and/or insert into database
# This is essentially just a command-line interface to rawBatch function in lib/httrpl.py module
#
# Usage:
#  rawBatch.py config.json
#

import sys
import os

# Make sure a parameter was specified
if len(sys.argv) < 2:
    sys.exit("No config file specified.")

config_file = sys.argv[1]

# Setup path to lib/ - assumes this script is in bin/ subdir next to lib/
TOP = '/'.join(os.path.abspath(__file__).split('/')[0:-2])
TOP = TOP+'/'
LIB = TOP+'lib'
if not LIB in sys.path:
    sys.path.insert(0,LIB)
os.environ['PYTHONPATH']=LIB

# Import Main API:
from httrpl import *

# Run everything from master function
rawBatch(config_file)
