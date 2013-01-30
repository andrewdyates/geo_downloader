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
import numpy as np
import itertools

OUT_FNAMES = {
  'data': "%s.seriesmatrix.tab",
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
  fp.write('ID_REF\t'); fp.write('\t'.join(g.platform.col_titles[1:])); fp.write('\n');
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


  # ==============================
  # GSM ATTRIBUTE HANDLING
  # ==============================

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

  key_list = list(attrs.keys())
  new_names = remove_prefixes(key_list)
  attr_name_map = dict(zip(key_list, new_names))
  print "Attribute name map after removing prefixes:", \
      ",".join(["%s=>%s"%(k,v) for k,v in attr_name_map.items()])

  # Attempt to merge attributes that seem to be split across multiple attributes.
  # collect attributes that have at least 5 missing values
  attr_masks = {}
  for attr_name in attrs:
    mask = [not bool(g.samples[s].attr.get(attr_name, False)) for s in g.col_titles[1:]]
    if np.sum(mask) >= 5:
      attr_masks[attr_name] = np.array(mask, dtype=np.bool)

  # merge pairs of attributes with disjoint masks and that have a name-corrected prefix of three characters
  mask_sets = [] # list of sets of attr keys
  if len(attr_masks) >= 2:
  # get equiv classes of prefixes
    prefixes = {}
    for k, mask in attr_masks.items():
      if len(k) < 3: continue
      pfx = attr_name_map[k][:3]
      prefixes.setdefault(pfx, set()).add(k)
      
    # for each equiv class, merge disjoint sets until no merge happens
    for pfx, keys in prefixes.items():
      mask_set = set()
      while True:
        if len(keys) == 1:
          break
        merges = False
        for k1, k2 in itertools.combinations(keys,2):
          # are masks disjoint?
          mask_1 = attr_masks[k1]
          mask_2 = attr_masks[k2]
          mask_u = mask_1&mask_2
          # match; masks are exactly disjoint. Merge masks, remove k2 from key set, repeat.
          if np.sum(~mask_1) + np.sum(~mask_2) == np.sum(~mask_u):
            attr_masks[k1] = mask_u
            mask_set.add(k1)
            mask_set.add(k2)
            keys.remove(k2)
            merges = True
            break
        if not merges:
          print "Cannot merge keys:", keys
          break
      if mask_set:
        mask_sets.append(mask_set)
  # report merge mask findings
  print "Mask merge results:"
  ignore_set = set()
  replacement_map = {}
  
  for keys in mask_sets:
    k = keys.pop()
    mask = attr_masks[k]
    keys.add(k)
    print "missing value mask size: %d, attributes: %s" % (np.sum(mask), ", ".join(list(keys)))
    # for longest string in mask set, set to merged values, delete other members from attr list
    best_key = sorted(keys, cmp=lambda q,r: len(q)<len(r), reverse=True)[0]
    keys.remove(best_key)
    ignore_set.update(keys)
    keys.add(best_key)
    print "best key for this set: %s. added %d other key(s) to ignore_set" % (best_key, len(keys)-1)
    # merge values
    merged_values = ['']*(len(g.col_titles)-1)
    for attr_name in keys:
      row = [",".join(g.samples[s].attr.get(attr_name, [])) for s in g.col_titles[1:]]
      for i, v in enumerate(row):
        if v:
          merged_values[i] = v
    # update values for best key
    replacement_map[best_key] = merged_values

  # ------------------------------
  # Write *.samples.tab
  # ------------------------------
  
  fp = open(out_fnames['samples'], 'w')
  # Write global attributes as comment headers.
  for k, v in global_attrs.items():
    fp.write("##%s=%s\n" % (k, v))
  # Write data column titles as sample headers.
  fp.write("GSM_ID\t"); fp.write("\t".join(g.col_titles[1:])); fp.write('\n')
  # Write attribute values in sample order
  n_wrote, n_replaced, n_ignored = 0, 0, 0
  for attr_name in attrs:
    if attr_name in ignore_set:
      n_ignored += 1
      continue
    # get merged values
    if attr_name in replacement_map:
      n_replaced += 1
      row = replacement_map[attr_name]
    else:
      row = [",".join(g.samples[s].attr.get(attr_name, [])) for s in g.col_titles[1:]]
    # write prefix-truncated name
    fp.write("%s\t" % attr_name_map[attr_name])
    # Write each attribute value in column order
    fp.write('\t'.join(row)); fp.write('\n')
    n_wrote +=1 
  fp.close()
  print "Wrote %d rows of %d columns, ignored %d rows, replaced %d rows (+%d headers, +1 attr name first column)" % \
      (n_wrote, len(g.col_titles)-1, n_ignored, n_replaced, len(global_attrs)+1)
  
  return out_fnames

def remove_prefixes(names):
  # From remaining multi-value attributes, remove any prefixes from attr names that are:
  #   over five characters long
  #   are shared by at least 3 and more than 2 attribute names
  #   if removed, all attribute names are still unique
  names = np.array(names)
  min_n = 3
  N = 5
  while True:
    prefixes = {}
    for i, name in enumerate(names):
      if len(name) > N:
        pfx = name[:N]
        prefixes.setdefault(pfx, set()).add(i)
    # continue if any prefix has enough names
    promising_pfxs = [k for k, v in prefixes.items() if len(v) >= min_n and len(v) > 2]
    
    # if no promising prefixes, break
    if not promising_pfxs:
      break
    
    for pfx in promising_pfxs:
      # expand prefixes
      idxs = list(prefixes[pfx])
      name_select = names[idxs]
      expand = 0
      while True:
        # All names must have remaining characters
        if any([True for s in name_select if len(s) <= N+expand+1]):
          break
        # All last characters should match. If not, break. Else, expand and loop.
        c=name_select[0][N+expand]
        for s in name_select[1:]:
          if s[N+expand]!=c:
            break
        if s[N+expand]!=c:
          break
        expand += 1
      full_pfx = name_select[0][:N+expand]
      print "Found prefix '%s' for %d words." % (full_pfx, len(name_select))
      name_stems = [s[N+expand:] for s in name_select]
      name_stem_set = set(name_stems)
      # does removing this prefix result in a duplicate name?
      if len(name_select) != len(name_stem_set):
        print "Stems are not all unique. Skipping..."
        continue
      # are any stems in the set of remaining names?
      if set(names) & name_stem_set:
        print "Stems found in complete name list. Skipping..."
      # Ok, remove prefixes from names.
      names[idxs] = name_stems
      
  return names


