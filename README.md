# QuantileNormalisation
# Description
This script takes as input files in tabular format and perfoms a quantile normalisation on the given data. Files have to contain 4 tab-seperated columns: chromosome, start position, stop position and raw reads. The files for quantile normalisation have to be giving as command line arguments, seperated by whitespace. The function "quantile_normalize” and parts of the code for visualisation were taken from a tutorial published on June 5, 2020 by cmdline (access dated: 08.10.2020 https://cmdlinetips.com/2020/06/computing-quantile-normalization-in-python/)

# Input data
The script takes a varible number of bedgraph files as input files (two or more). It produces one output file per input file, containing the same columns as the input file with the normalised instead of the raw data.
