geo_downloader
==============
Download NCBI GEO study data as "clean" text files.

Example Use
-----------

    python script.py gse_id=GSE15745 platform_id=GPL6104 outdir=$HOME/Desktop

Script Parameters
-----------

* gse_id: [required] str of GEO GSE Study Accession ID (uppercase)
* platform_id: str of GEO GPL Platform Accession ID;
    required only when multiple platforms are associated with the same study
* outdir: str of path to output directory. Defaults to current working directory.
* normalize: [default=True] bool if to perform quantile normalization on data matrix. Disabling normalize also prevents parsing and loading the text matrix as a numpy object. Thus, `*.normed.tab`, `*.normed.masked.pkl`, and `*.varlist.txt` will also not be saved. Set `normalize=False` to conserve memory and processor time in a constrained environment like a web server.

About
-----------

GEO Downloader fetches "series matrix" data files from the [NCBI Gene Expression Omnibus (GEO)](http://www.ncbi.nlm.nih.gov/geo/)
public database as easy to work with plain text files. Matrices are aligned so that rows in
the probes (GPL platform definition) matrix align with corresponding rows in the data matrix, and
columns in the samples (GSM sample attribute and metadata) matrix align with corresponding columns
in the data matrix. The data matrix itself is saved in three forms: raw, normed, and numpy normed.
See Table 1.
    
GEO Downloader uses the GEO API and thus, by proxy, the CACHE_DIR operating system
environment variable. If CACHE_DIR is not set, GEO Downloader sets CACHE_DIR equal
to the output directory `outdir.` Any file downloaded from the Internet will be saved
into CACHE_DIR as a temporary file. 


### Figure 1: Matrix Alignments ###

             +---------+
             | samples |
    +--------+ --------+
    | probes |  data   |
    +--------+---------+


### Table 1: Output Files ###
_saved to `outdir`, full file name is generated using `gse_id` and `platform_id`_

| File Suffix | Matrix  | Type               | Notes |
| ----------- | ------  | ----               | ----- |
| *.raw.tab     | data    | tab delimited text | Concatenated series matrix data without transformation. |
| *.normed.tab  | data    | tab delimited text | Quantile normed 'raw.tab' in same text format. |
| *.normed.masked.pkl  | data    | pickled numpy.ma binary    | Quantile normed 'raw.tab' as binary; no row labels. Load back into a numpy.MaskedArray (ma) object using cPickle.load |
| *.probes.tab  | probes  | tab delimited text | Platform "GPL" probe definition per probe in data row order. |
| *.varlist.txt  | data  | text | List of probe (row variable) names in row order from top to bottom, one per line. |
| *.gpl_brief.txt  | probes  | SOFT text | GPL metadata file "Brief" as downloaded from the GEO website; no row descriptions. |
| *.samples.tab | samples | tab delimited text | All sample "GSM" attributes with at least two unique values in data column order. Attributes with only one value are included as '#' prefixed headers at the top of the file. |

