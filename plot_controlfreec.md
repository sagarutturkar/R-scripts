## Background:
[Control-FREEC](http://boevalab.inf.ethz.ch/FREEC/index.html) is a tool for assessing copy number and allelic content 
using next generation sequencing data. The tool provides the script [makeGraph.R](https://github.com/BoevaLab/FREEC/blob/master/scripts/makeGraph.R) 
to visualize normalized copy number profile with predicted Copy Number Alterations (CNAs). However, this script generate different plots for different chromosomes
and the Y-scale could be different for each chromosome based on the distribution of the data.

## Motivation:
For the genome-wide CNA visualization, we needed a scaled plot where all chromosomes are plotted with uniform Y-scale and combined plot for
all given chromosomes. This script is inspired by the original [makeGraph.R](https://github.com/BoevaLab/FREEC/blob/master/scripts/makeGraph.R)
and attempts to achieve combined chromosome plot with uniform Y-scale.

