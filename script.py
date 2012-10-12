#!/usr/bin/python
"""Shell script wrapper for `__init__.download()`.
EXAMPLE USE:
  python script.py gse_id=GSE15745 platform_id=GPL6104 outdir=$HOME/Desktop
"""
from __init__ import *

if __name__ == "__main__":
  print sys.argv
  download(**dict([s.split('=') for s in sys.argv[1:]]))
