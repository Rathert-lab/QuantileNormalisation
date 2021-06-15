# QuantileNormalisation
## Quantile Normalisation of .bedgraph files

Performs quantile normalization of ChIP-seq read data stored in bedgraph files using Python. The function "quantile_normalize” and parts of the code for visualisation were taken from a tutorial published on June 5, 2020 by cmdline (access dated: 08.10.2020 https://cmdlinetips.com/2020/06/computing-quantile-normalization-in-python/).

## Usage
The script requires an installation of Python to run and five different packages. These will be imported when running the script but might 
- sys
- matplotlib.pyplot
- numpy
- pandas
- seaborn

The files for quantile normalisation have to be giving as command line arguments, seperated by whitespace. The script can be executed from a Linux command terminal. Use the following command to run the script:

```python QuantileNormalisation.py input1 input2 input3
```

The number of input files is variable from two up to infinite files in principle. It was tested for up to eight input files. Using more files than 10 should be possible but will take very long to compute.
Output files will be gerneated automatically and saved as input1_normalised.tabular, input2_normalised.tabular, ... and so on.

There are to optional graphical outputs: density graphs or boxplots of the raw vs. quantile normalised data. By default the script produces box-plots, but the flag for density_graph in the script can be chnaged to "True" before running the script to produce density graphs as well.

## Input Data
This script takes as input files in tabular format, like .bedgraph files, and perfoms a quantile normalisation on the given data. Files have to contain 4 tab-seperated columns: chromosome, start position, stop position and raw reads. All files for normalization should contain the same number of datapoints (same bin size and same genome).
```
chromA  chromStartA  chromEndA  dataValueA
chromB  chromStartB  chromEndB  dataValueB
```


## License

MIT
