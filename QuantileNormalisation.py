#===============================================================================================#
#                                                                                               #
#   Title:       QuantileNormalisation of ChIP-Seq Data                                         #
#   Author:      Tabea Bauer                                                                    #
#   Email:       Tabea.Bauer@ibtb.uni-stuttgart.de                                              #
#   Description: This script takes as input data in bedgraph format and perfoms a quantile      # 
#                normalisation on the given data. The files for quantile normalisation have     #
#                to be giving as command line arguments. The function "quantile_normalize"      #
#                and parts of the code for visualisation were taken from a tutorial             #
#                published on June 5, 2020 by cmdline (access dated: 08.10.2020 ;               #
#                https://cmdlinetips.com/2020/06/computing-quantile-normalization-in-python/)   #
#   Data:        The script takes a varible number of bedgraph files as input files (two or     #
#                more). It produces one output file per input file, containing the same         #
#                columns as the input file with the normalised instead of the raw data.         #
#   Version:     V1.1                                                                           #
#   LastChange:  09.10.2020                                                                     #
#                                                                                               #
#===============================================================================================#  

# import all needed packages
import sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

# determine wheter density and/or boxplot graphs should be created
density_graph = False
boxplot_graph = True

# define empty lists needed for the script to rnu
reads, chrom, start, stop, files = ([] for i in range(5)) 

"""
In this section, the files given as command line arguments are read in
"""
for f in sys.argv[1:]:                                                  # perform same action for each given file
    files.append(f.rpartition('.')[0])                                  # write file names (without extension) into list
    list1 = np.loadtxt(f, dtype='str', delimiter="\t", usecols=[0,1,2]) # save the first 3 columns of the file (chrom, startposition and stopposition) into a list
    list2 = np.loadtxt(f, dtype='int', delimiter="\t", usecols=3)       # write raw data into a list
    reads.append(list2.tolist())                                        # append data of the files to create a list of lists containing all the raw reads of all files
    chrom.append(list1[:, 0].tolist())                                  # ---------------------"--------------------------------------- chromosome entries of all files
    start.append(list1[:, 1].tolist())                                  # ---------------------"--------------------------------------- startpositions of all files
    stop.append(list1[:, 2].tolist())                                   # ---------------------"--------------------------------------- stoppositions of all files
    list1 = []                                                          # empty the list for the next file
    list2 = []                                                          # empty the list for the next file

# save the created list of lists in arrays, each row represents data of one input file
input = np.array(reads)
inputchrom = np.array(chrom)
inputstart = np.array(start)
inputstop = np.array(stop)

# data is saved into Dataframes, so that is can be written to a file and save eacsier later, for this 
# the arrays need to be transformed, so that each column contains the data from one file, rather than each row

data = pd.DataFrame(input.transpose(), columns=files)
datachrom = pd.DataFrame(inputchrom.transpose(), columns=files)
datastart = pd.DataFrame(inputstart.transpose(), columns=files)
datastop = pd.DataFrame(inputstop.transpose(), columns=files)

# function that performs the quantile normalisation as taken from the tutorial by "cmdline"
def quantile_normalize(df):

    """
    input: dataframe with numerical columns
    output: dataframe with quantile normalized values
    """

    df_sorted = pd.DataFrame(np.sort(df.values,
                                     axis=0),
                             index=df.index,
                             columns=df.columns)
    df_mean = df_sorted.mean(axis=1)
    df_mean.index = np.arange(1, len(df_mean) + 1)
    df_qn =df.rank(method="min").stack().astype(int).map(df_mean).unstack()
    return(df_qn)

# perform normalisation of data with the quantile normalisation function
normed = quantile_normalize(data)

"""
Data visualisation
"""

# This paragraphs builds density graphs before and after normalisation

if density_graph is True:
    fig = plt.figure(figsize=[12.8, 4.8])
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    data.plot.density(ax=ax1,linewidth=4)
    ax1.set_title("Density plot before Quantile Normalization")
    ax1.set_xlabel("Measurement", size=14)
    ax1.set_ylabel("Density", size=14)

    normed.plot.density(ax=ax2,linewidth=4)
    ax2.set_title("Density plot after Quantile Normalization")
    ax2.set_xlabel("Measurement", size=14)
    ax2.set_ylabel("Density", size=14)
    plt.savefig('Density_plot_Quantile_Normalization_Pandas.png', dpi=150)
    #plt.show()

# This paragraph builds boxplots before and after normalisaiton

# calculate potions of whispers to get plot limits

def plotlimits(columnnames,values):

    """
    input: array containing the column names and Dataframe array containing the data (values)
    output: array containing the upper and lower limit for the plot
    """
    whiskerup=[]
    whiskerlo=[]
    for n in columnnames:
        IQR = np.percentile(values[n], [75]) - np.percentile(values[n], [25])
        whiskerup.append(np.percentile(values[n], [75]) + 1.6 * IQR)
        if (min(values[n]) < 0):
            whiskerlo.append(np.percentile(values[n], [25]) - 1.6 * IQR)
        else:
            whiskerlo.append(-50)
    limits = [min(whiskerlo),max(whiskerup)]
    return(limits)


ylim1 = plotlimits(files,data)
ylim2 = plotlimits(files,normed)

if boxplot_graph is True:
    fig = plt.figure(figsize=[12.8, 4.8])
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    sns.boxplot(data=data, ax=ax1, whis=1.5)
    # set x-axis label
    ax1.set_xlabel("Samples", size=14)
    # set y-axis label
    ax1.set_ylabel("Measurement", size=14)
    ax1.set_ylim(ylim1)
    ax1.set_title("Boxplot of raw data before Quantile Normalization")

    sns.boxplot(data=normed, ax=ax2, whis=1.5)
    # set x-axis label
    ax2.set_xlabel("Samples", size=14)
    # set y-axis label
    ax2.set_ylabel("Measurement", size=14)
    ax2.set_ylim(ylim2)
    ax2.set_title("Boxplot after Quantile Normalization")
    plt.savefig('Boxplot_Quantile_Normalization_Seaborn.png', dpi=150)

plt.show()

for file in files:
    f = open(file+"_normalised.tabular", "w")
    frames = [datachrom[file], datastart[file], datastop[file], normed[file]]
    output = pd.concat(frames,axis=1)
    output.to_csv(f, sep='\t', header=False, index=False, line_terminator='\n')

