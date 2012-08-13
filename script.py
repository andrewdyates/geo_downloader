#!/usr/bin/python
"""Save study to disk as ordered and aligned tab delimited files.

EXAMPLE USE:
python script.py gse_id=GSE15745 platform_id=GPL6104 outdir=$HOME/Desktop

raw
normed
samples
probes
"""
import os
import numpy as np

from geo_api import *
from quantile_normalize import *

OUT_FNAMES = {
  'raw_fname': "%s.raw.tab",
  'normed_fname': "%s.normed.tab",
  'samples_fname': "%s.samples.tab",
  'probes_fname': "%s.probes.tab",
}

def main(gse_id=None, outdir=None, platform_id=None):
  """Save study information to disk."""
  assert gse_id is not None
    
  # Configure output and download cache
  if outdir is None:
    outdir = os.getcwd()
    print "Warning: outdir not specified. Set outdir to current working directory %s." % (outdir)
  if "CACHE_DIR" not in os.environ:
    print "Warning: os enviroment variable CACHE_DIR not set. Setting CACHE_DIR to `outdir` %s" % (outdir)
    os.environ["CACHE_DIR"] = outdir

  # Generate output file paths.
  for k,v in OUT_FNAMES.items():
    OUT_FNAMES[k] = os.path.join(outdir, v%gse_id)

  # Fetch study.
  g = GSE(gse_id, platform_id=platform_id)

  # Populate g and save raw study data.
  row_ids = []
  print "Saving raw study data to %s." % (OUT_FNAMES['raw_fname'])
  fp = open(OUT_FNAMES['raw_fname'], 'w')
  fp.write('#'); fp.write('\t'.join(g.col_titles)); fp.write('\n');
  for row in g.get_rows():
    fp.write('\t'.join(row)); fp.write('\n');
    row_ids.append(row[0])
  fp.close()
  print "Wrote %d rows of %d columns (+1 header, includes ID column)" % \
      (len(row_ids), len(row))

  # Save GPL platform definitions in study data row order.
  print "Saving probe row definitions in row order to %s." % (OUT_FNAMES['probes_fname'])
  print g.platform.col_titles
  fp = open(OUT_FNAMES['probes_fname'], 'w')
  fp.write('#ID_REF\t'); fp.write('\t'.join(g.platform.col_titles[1:])); fp.write('\n');
  for s in row_ids:
    d, row = g.platform.row_desc[s], [s]
    for k in g.platform.col_titles[1:]:
      if k in d:
        row.append(d[k])
      else:
        row.append("")
    fp.write('\t'.join(row)); fp.write('\n')
  fp.close()
  print "Wrote %d rows of %d columns (+1 header, includes ID column)" % \
      (len(row_ids), len(g.platform.col_titles))

  # Save sample meta in study data column order.
  print "Saving sample meta in column order to %s." % (OUT_FNAMES['samples_fname'])

  # Get all GSM attributes
  attrs = {}
  for s, gsm in g.samples.items():
    for k, v in gsm.attr.items():
      for x in v:
        attrs.setdefault(k, set()).add(x)

  for k, v in attrs.items():
    print k, v
    
  #fp = open(OUT_FNAMES['samples_fname'], 'w')

  
if __name__ == "__main__":
  main(**dict([s.split('=') for s in sys.argv[1:]]))
