![1429-bacteria banner](https://github.com/user-attachments/assets/49bc62ec-13a4-45f9-8895-c225339edac7)

# Scientific Project: Difference in Microbiome Compositions

This repository houses the computational code developed for the scientific programming project (2425-MSB1015). The project aims to investigate the microbiome composition of individuals with Autism Spectrum Disorder (ASD) in comparison to controls.

The code (ScientificProjectV2.R) implements a pipeline that encompasses data processing, cleaning, and statistical analysis. By examining microbiome data at both the Phylum (higher-level) and Family (lower-level) taxonomic levels, this analysis provides a balanced perspective on potential differences in microbiome composition between these groups.

#### There are two datasets used:
#### GSE113690_Autism_16S_rRNA_OTU_assignment_and_abundance.csv

  From: NCBI study GSE113690

  Contains: 1322 OTUs, 254 samples, taxonomy info and count values 

  Link: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE113690 

#### SampleInfo.xlsx

From: Supplemental Data of GSE113690 paper

  Contains: Sample ID, Stage, Gender, Age, etc.

  Original file name: Table S1 Sample information.xlsx

  Link to paper: https://doi.org/10.1080/19490976.2020.1747329

## Prerequisites
This software has to be installed:

[R 4.4.1](https://cran.r-project.org/bin/windows/base/)

[RStudio](https://posit.co/download/rstudio-desktop/)

[Hands-On installation guide](https://rstudio-education.github.io/hopr/starting.html) 

## Installation
Clone the git repository: 
```
git clone https://github.com/MaxMourmans/Scientific-Project
```

## R packages/libraries to install
BiocManager is required to install the following packages:
- tidyverse
- glue
- gtsummary
- ggplot2
- phyloseq
- readxl
- ggpubr
- microbiome
- microViz

## Usage
Once the necessary prerequisites are installed, the project can be accessed in RStudio. The following general steps outline the data analysis process after importing the datasets:

### (1)	Data processing:
What’s done: Important information is extracted and set to the right format

Plots/table: Characteristics table based on sample information
### (2)	General data cleaning:
What’s done: Duplicates, missing values and mislabels are identified/corrected

Plots/table: none
### (3)	Statistical tests:
What’s done: Phyloseq object created, Alpha- and Beta diversity tests

Plots/table: Composition plots, Alpha/Beta boxplots, PCoA

## Warning 
##### The verification functions for correcting mislabels don't operate anymore, due to deletion of  the package [Taxize](https://cran.r-project.org/web/packages/taxize/index.html) 

The script can be simulated without impact on the statisitical results!

## License
See License for license rights and limitations [MIT](https://choosealicense.com/licenses/mit/)

## Contact

Max Mourmans - maxmourmans22@gmail.com

Project link - https://github.com/MaxMourmans/Scientific-Project 


