#!/usr/bin/python
"""Shell script wrapper for `__init__.download()`.
"""
from __init__ import *

if __name__ == "__main__":
  print sys.argv
  download(**dict([s.split('=') for s in sys.argv[1:]]))
