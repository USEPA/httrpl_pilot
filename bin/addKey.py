#!/usr/bin/python3
#
# Add/update a DB login key to the user's private keychain before working with a new database
# This is just a command-line interface to run addKeychainEntry function in fully interactive mode
#
# Usage:
#  addKey.py
#

import sys
import os

# Setup path to lib/ - assumes this script is in bin/ subdir next to lib/
TOP = '/'.join(os.path.abspath(__file__).split('/')[0:-2])
TOP = TOP+'/'
LIB = TOP+'lib'
if not LIB in sys.path:
    sys.path.insert(0,LIB)
os.environ['PYTHONPATH']=LIB

# Import Main API:
from httrpl import *

# Run addKeychainEntry interactively:
addKeychainEntry(interactive=True)
