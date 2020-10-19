#!/usr/bin/python3
#
# Tests that all python modules and dependencies can be loaded in the current environment
#
# Usage:
#  python3 testInstall.py
#

import sys, os

# Setup path to lib/ - assumes this script is in bin/ subdir next to lib/
TOP = '/'.join(os.path.abspath(__file__).split('/')[0:-2])
TOP = TOP+'/'
LIB = TOP+'lib'
if not LIB in sys.path:
    sys.path.insert(0,LIB)
os.environ['PYTHONPATH']=LIB

# Import all httrpl libs to make sure everything loads properly
from httrpl import *
from httrplcore import *
from db.mongo import *
from db.raw import *
from gexp.fastq import *
from gexp.biospyder2 import *

stdoutlog = PipelineLogger(out=sys.stdout, dbg=True, strict=False)
envReport(log=stdoutlog)
