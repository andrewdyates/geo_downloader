#!/usr/bin/python
"""Save study to disk as ordered and aligned tab delimited files.

EXAMPLE USE:

  python script.py gse_id=GSE15745 platform_id=GPL6104 outdir=$HOME/Desktop
  python script.py gse_id=GSE7307 outdir=$HOME/Desktop

See README.md for details.
"""
import os
from geo_api import *
from lab_util import *

OUT_FNAMES = {
  'data': "%s.data.tab",
  'probes': "%s.probes.tab",
  'samples': "%s.samples.tab",
  'varlist': "%s.varlist.txt",
  'gpl_brief': "%s.gpl_brief.txt",
}

def download(gse_id=None, platform_id=None, outdir=None):
  """Save study information to disk.

  Returns:
    {str: str} of {type: file_path} of files saved.
  """
  assert gse_id is not None
  gse_id = gse_id.upper()
  if platform_id:
    platform_id = platform_id.upper()
  out_fnames = OUT_FNAMES.copy()

  # Configure output destination, download cache, and download logging
  if outdir is None:
    outdir = os.getcwd()
    print "Warning: outdir not specified. Set outdir to current working directory %s." % (outdir)
  if not os.path.exists(outdir):
    print "Creating %s..." % (outdir)
    make_dir(outdir)
  if "CACHE_DIR" not in os.environ:
    print "Warning: os enviroment variable CACHE_DIR not set. Setting CACHE_DIR to `outdir` %s" % (outdir)
    os.environ["CACHE_DIR"] = outdir

  # Fetch study.
  g = GSE(gse_id, platform_id=platform_id)
  if platform_id is not None:
    assert g.platform.id == platform_id
    
  # Generate output file paths.
  for k,v in out_fnames.items():
    name = "%s_%s" % (gse_id, g.platform.id)
    out_fnames[k] = os.path.join(outdir, v % name)

  # Save untransformed series matrix study data as text. Populate geo_api object g.
  row_ids = []
  print "Saving combined series matrix study data to %s." % (out_fnames['data'])
  fp = open(out_fnames['data'], 'w')
  fp.write('#'); fp.write('\t'.join(g.col_titles)); fp.write('\n');
  for row in g.get_rows():
    fp.write('\t'.join(row)); fp.write('\n');
    row_ids.append(row[0])
  fp.close()
  print "Wrote %d rows of %d columns (+1 header, includes ID column)" % \
      (len(row_ids), len(row))

  # Save GPL platform definitions in study data row order.
  print "Saving probe row definitions in row order to %s." % (out_fnames['probes'])
  print g.platform.col_titles
  fp = open(out_fnames['probes'], 'w')
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
  
  # Save original GPL file header.
  fp = open(out_fnames['gpl_brief'], 'w')
  for line in GPL.fp_download_brief(g.platform.id):
    fp.write(line)
  fp.close()
  print "Wrote original GPL brief file from %s" % (g.platform.brief_url)

  # Save sample meta in study data column order.
  print "Saving sample meta in column order to %s." % (out_fnames['samples'])

  # Get all GSM attributes
  attrs = {}
  for s, gsm in g.samples.items():
    for k, v in gsm.attr.items():
      for x in v:
        attrs.setdefault(k, set()).add(x)
        
  # Select attributes with at least 2 values. Print others as comments.
  global_attrs = {}
  for k, v in attrs.items():
    if len(v) == 1:
      global_attrs[k] = v.pop()
  for k in global_attrs:
    del attrs[k]
  print "Selected %d sample attributes and %d global attributes." % (len(attrs), len(global_attrs))
  
  fp = open(out_fnames['samples'], 'w')
  # Write global attributes as comment headers.
  for k, v in global_attrs.items():
    fp.write("##%s=%s\n" % (k, v))
  # Write data column titles as sample headers.
  fp.write("GSM_ID\t"); fp.write("\t".join(g.col_titles[1:])); fp.write('\n')
  # Write attribute values in sample order
  for attr_name in attrs:
    fp.write("%s\t" % attr_name)
    # Write each attribute value in column order
    row = [",".join(g.samples[s].attr.get(attr_name, [])) for s in g.col_titles[1:]]
    fp.write('\t'.join(row)); fp.write('\n')
  fp.close()
  print "Wrote %d rows of %d columns (+%d headers, includes attr name column)" % \
      (len(attrs), len(g.col_titles)-1, len(global_attrs)+1)
  
  return out_fnames
