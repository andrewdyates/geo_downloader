geo_downloader
==============
Download NCBI GEO study data as "clean" text files.

Example Use
-----------

    python script.py gse_id=GSE15745 platform_id=GPL6104 outdir=$HOME/Desktop

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


| File Suffix | Matrix  | Type               | Notes |
| ----------- | ------  | ----               | ----- |
| *.raw.tab     | data    | tab delimited text | Concatenated series matrix data without transformation. |
| *.normed.tab  | data    | tab delimited text | Quantile normed 'raw.tab' in same text format. |
| *.normed.npy  | data    | numpy.ma binary    | Quantile normed 'raw.tab' as binary; no row labels. |
| *.probes.tab  | probes  | tab delimited text | Platform "GPL" probe definition per probe in data row order. |
| *.samples.tab | samples | tab delimited text | All sample "GSM" attributes with at least two unique values in data column order. Attributes with only one value are included as '#' prefixed headers at the top of the file. |

