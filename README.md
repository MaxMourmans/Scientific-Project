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

## R packages/libraries to install
- tidyverse
- glue
- gtsummary
- ggplot2
- taxize
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

## Important when running 
#### DO NOT RUN THE SCRIPT IN ONE GO!

The code includes two verification functions: one for Phyla and one for Family. These functions can be used to validate the data at the specified taxonomic level. Be aware to run these 'for' loops solely!

```R
#Apply Verification Phyla Function
## WARNING: RUN THIS FOR-LOOP SEPERATELY AND GIVE INPUT WHERE NEEDED!
for (i in seq_along(unique_phyla)) {
  check_phylum(unique_phyla[i], i)
}
```

To run the Phyla verification, you'll be prompted for input. Simply enter '1' or '2' in the Console and press Enter. The correct sequence for this verification is: 2, 2, 2, 1.

```R
  status   rank     division scientificname               commonname    uid genus species subsp modificationdate
1 active phylum GNS bacteria  Chloroflexota green nonsulfur bacteria 200795                     2023/02/13 00:00
2 active  class GNS bacteria   Chloroflexia                           32061                     2018/10/26 00:00
1: _?
```

## License
See License for license rights and limitations [MIT](https://choosealicense.com/licenses/mit/)

## Contact

Max Mourmans - maxmourmans22@gmail.com

Project link - https://github.com/MaxMourmans/Scientific-Project 


