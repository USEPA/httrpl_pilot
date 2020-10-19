#!/usr/bin/python3
#
# Run alignment and probe counting step on a batch of fastq files, then output results to a TSV file and/or database
# This is essentially just a command-line interface to countBatch function in lib/httrpl.py module
#
# Usage:
#  countBatch.py count_config.json
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

# Check versions:
envReport()

# Run everything from master function
countBatch(config_file)
