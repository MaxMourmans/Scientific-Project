##### SCIENTIFIC PROGRAMMING #####
# Subject: Meta-genomic sequencing of gut microbiome associated with ASD patients
# Author: Max Mourmans
# Date: 20-10-2024
# Info: Project code for Scientific Programming Project, MSB


####################
# Install packages #
####################
#Skip installation part, when already installed all packages below
if (!requireNamespace("BiocManager", quietly = TRUE)) {
 install.packages("BiocManager")
}
#Install R packages using BiocManager
BiocManager::install("tidyverse")
BiocManager::install("glue")
BiocManager::install("gtsummary")
BiocManager::install("ggplot2")
BiocManager::install("phyloseq")
BiocManager::install("readxl")
BiocManager::install("ggpubr")
BiocManager::install("microbiome")
BiocManager::install("microViz",
repos = c(davidbarnett = "https://david-barnett.r-universe.dev", getOption("repos"))
)
#TAXIZE DOESN'T WORK ANYMORE!
#BiocManager::install("taxize")

####################
#    R libraries   #
####################
library(tidyverse)
library(glue)
library(gtsummary)
library(ggplot2)
library(phyloseq)
library(readxl)
library(ggpubr)
library(microbiome)
library(microViz)
#TAXIZE DOESN'T WORK ANYMORE!
#library(taxize)

#If not working in GitHub cloned repository, uncomment and set working directory
setwd("C:/Set/Your/Working/Directory/Here")

#Load data set
raw_rRNA_OTU <- read.csv("16S_rRNA_OTU_assignment_and_abundance.csv")
metadata <- read_excel("SampleInfo.xlsx")

# Data Sorting ------------------------------------------------------------

#Check if there are any duplicate samples/missing values in metadata
#Give 'TRUE' when there are any duplicate samples
any(duplicated(metadata$'Sample ID' == TRUE))
#Give 'TRUE' when there are any missing values in metadata
any(missing(metadata) == TRUE)
#Set Sample IDs as row names
metadata <- metadata %>% 
  tibble::column_to_rownames("Sample ID")
#Create a new column with 3 age categories
metadata <- metadata %>%
  mutate(ageCat = case_when(
    Age >= 2 & Age <= 3  ~ "2-3",
    Age >= 4 & Age <= 6  ~ "4-6",
    Age >= 7 & Age <= 11 ~ "7-11",
    Age >= 12 & Age <= 13 ~ "12-13",
    TRUE ~ NA_character_  # if age does not fit into any category
  ))

#Reorder categories, to have them in the right order for characteristics table
metadata$ageCat <- factor(metadata$ageCat, levels = c("2-3", "4-6", "7-11", "12-13"))
#Create a characteristics table using metadata
  tbl_summary(
    metadata,
    include = c(Gender, Age, ageCat, Constipation),
    type = Age ~ "continuous2",
    statistic = 
      list(Age ~ c("{mean} ({sd})","{min} - {max}")),
    label = ageCat ~ "Age Categories",
    by = Stage,
    digits = Age ~ 1,
  ) |> 
  modify_header(label = "**Patients**") |> 
  bold_labels()

#Extract phylum taxa-level
extracted_phylum <- data.frame(gsub(".*p__(.*?)?;_.*", "\\1", raw_rRNA_OTU$taxonomy))
colnames(extracted_phylum)[1] <- "Phylum"
#Extract family taxa-level
extracted_family <- data.frame(gsub(".*_f__(.*?)?;_.*", "\\1", raw_rRNA_OTU$taxonomy))
colnames(extracted_family)[1] <- "Family"
#Drop taxonomy column
rRNA_OTU = raw_rRNA_OTU[,-2]
rRNA_OTU = cbind(rRNA_OTU, NewCol = extracted_phylum)
rRNA_OTU = cbind(rRNA_OTU, NewCol = extracted_family)
#Reorder phylum, family columns in data matrix
rRNA_OTU = rRNA_OTU[, c(1, 256:257, 2:255)]

# General Data Cleaning --------------------------------------------------------

# 1. Are there any sample duplicates?
#Give 'TRUE' when there are any duplicate samples in OTU data
any(duplicated(rRNA_OTU) == TRUE)
#Give 'TRUE' when there are any duplicate ID's in OTU data
any(duplicated(colnames(rRNA_OTU)) == TRUE)

# 2. Are there any missing values(NAs) in OTU data?
any(missing(rRNA_OTU) == TRUE)

# 3.1 Are there any mislabels/unclassified labels in Phylum labeling?
#Find all unique ID's in Phylum column
unique_phyla <- unique(rRNA_OTU$Phylum)
#If there are any phyla duplicates this statement is TRUE
any(duplicated(unique_phyla) == TRUE)

#Create data frame to allocate verification status of all phyla
verified_phyla <- data.frame(verified_status = rep(0, length(unique_phyla)))

#Verification Phyla Function:
check_phylum <- function(phylum, index) {
  result <- tryCatch({
    #Query classification hierarchy from NCBI database
    res <- classification(phylum, db = "ncbi")
    
    #Check if the result includes "phylum" as a rank
    if (!is.null(res[[1]]) && any(res[[1]]$rank == "phylum")) {
      verified_phyla$verified_status[index] <<- "Verified"  #Update verified_phyla with "Verified"
    } else {
      verified_phyla$verified_status[index] <<- "Not Verified"  #Update verified_phyla with "Not Verified"
    }
  }, error = function(e) {
    verified_phyla$verified_status[index] <<- "Not Verified"  #Handle errors by updating verified_phyla
  })
}

# #Apply Verification Phyla Function
# ## WARNING: DOESN'T WORK ANYMORE DUE TO TAXIZE
# for (i in seq_along(unique_phyla)) {
#   check_phylum(unique_phyla[i], i)
# }
# rm(i)
#Data frame with all phyla and verification status
verified_phyla <- cbind(data.frame(unique_phyla), verified_phyla)


#Correct not verified data by renaming the mislabels based on NCBI data
rRNA_OTU <- rRNA_OTU %>%
  mutate(Phylum = case_when(
    Phylum == "Saccharibacteria" ~ "Candidatus Saccharibacteria",
    Phylum == "SR1__Absconditabacteria_" ~ "Candidatus Absconditabacteria",
    TRUE ~ Phylum  # Keep the original name if no typo
  ))
#Remove unclassified phyla
rRNA_OTU <- rRNA_OTU %>%
  filter(!grepl("unclassified_k__norank", Phylum))

#Run again to verify the mislabel corrections by repeating the above steps
unique_phyla <- unique(rRNA_OTU$Phylum)
verified_phyla <- data.frame(verified_status = rep(0, length(unique_phyla)))
# ## WARNING: DOESN'T WORK ANYMORE DUE TO TAXIZE
# for (i in seq_along(unique_phyla)) {
#   check_phylum(unique_phyla[i], i) 
# }
#rm(i)
#Data frame with Phyla and UPDATED verification status
verified_phyla <- cbind(data.frame(unique_phyla), verified_phyla)



# 3.2 Are there any mislabels/unclassified labels in Family labeling?
#Find all unique ID's in Family column
unique_families <- unique(rRNA_OTU$Family)
#If there are any families duplicates this statement is TRUE
any(duplicated(unique_families) == TRUE)

#Create data frame to allocate verification status of all families
verified_families <- data.frame(verified_status = rep(0, length(unique_families)))

#Verification Family Function:
check_family <- function(family, index) {
  result <- tryCatch({
    #Query classification hierarchy from NCBI database
    res <- classification(family, db = "ncbi")
    
    #Check if the result includes "family" as a rank
    if (!is.null(res[[1]]) && any(res[[1]]$rank == "family")) {
      verified_families$verified_status[index] <<- "Verified"  #Update verified_families with "Verified"
    } else {
      verified_families$verified_status[index] <<- "Not Verified"  #Update verified_families with "Not Verified"
    }
  }, error = function(e) {
    verified_families$verified_status[index] <<- "Not Verified"  #Handle errors by updating verified_families
  })
}

#Apply Verification Family Function
# ## WARNING: DOESN'T WORK ANYMORE DUE TO TAXIZE
# for (i in seq_along(unique_families)) {
#   check_family(unique_families[i], i)
# }
#rm(i)
#Data frame with all families and verification status
verified_families <- cbind(data.frame(unique_families), verified_families)


#Remove unclassified phyla
rRNA_OTU <- rRNA_OTU %>%
  filter(!grepl("norank", Family))
rRNA_OTU <- rRNA_OTU %>%
  filter(!grepl("unclassified", Family))
rRNA_OTU <- rRNA_OTU %>%
  filter(!grepl("Unknown", Family))
#Remove unclassified phyla
rRNA_OTU <- rRNA_OTU %>%
  filter(!grepl("Mitochondria", Family))

#Correct mislabels
rRNA_OTU <- rRNA_OTU %>%
  mutate(Family = case_when(
    Family == "Rubrobacteriaceae" ~ "Rubrobacteraceae",
    Family == "Bacteroidales_S24-7_group" ~ "Muribaculaceae",
    Family == "Clostridiales_vadinBB60_group" ~ "Clostridiaceae",
    Family == "Family_XI_o__Clostridiales" ~ "Clostridiaceae",
    Family == "Family_XI" ~ "Clostridiaceae",
    Family == "Clostridiaceae_1" ~ "Clostridiaceae",    
    Family == "Family_XIII" ~ "Clostridiales Family XIII. Incertae Sedis",   
    Family == "Family_XI_o__Bacillales" ~ "Bacillaceae",
    Family == "PHOS-HE36" ~ "Ignavibacteriaceae",
    Family == "NS9_marine_group" ~ "Flavobacteriaceae",
    Family == "vadinBE97" ~ "Victivallaceae",
    Family == "Burkholderiales_Incertae_Sedis" ~ "Burkholderiaceae",
    Family == "Rickettsiales_Incertae_Sedis" ~ "Rickettsiaceae",
    Family == "Bacteroidales_Incertae_Sedis" ~ "Bacteroidaceae",
    TRUE ~ Family  # Keep the original name if no typo
  ))

#Run again to verify the mislabel corrections by repeating the above steps
unique_families <- unique(rRNA_OTU$Family)
verified_families <- data.frame(verified_status = rep(0, length(unique_families)))
### WARNING: DOESN'T WORK ANYMORE DUE TO TAXIZE
# for (i in seq_along(unique_families)) {
#   check_family(unique_families[i], i) 
# }
#rm(i)
#Family and UPDATED verification status
verified_families <- cbind(data.frame(unique_families), verified_families)



# Statistical tests -------------------------------------------------------

#Creating of a Phyloseq object using rRNA_OTU
otu_table <- rRNA_OTU
row.names(otu_table) <- rRNA_OTU$OTU
otu_table <- otu_table[,-c(1:3)]
OTU <- otu_table(as.matrix(otu_table), taxa_are_rows = TRUE)
#Make Taxa matrix using column Phylum from rRNA_OTU
tax_table <- data.frame(Phylum = rRNA_OTU[,c(2)], Family = rRNA_OTU[,c(3)] )
row.names(tax_table) <- rRNA_OTU[,1]
TAX <- tax_table(as.matrix(tax_table))
#Make Sample information table based on groups from meta_autism
SAMPLES <- sample_data(metadata)

#Use all these matrices to create phyloseq object: physeq
physeq <- phyloseq(OTU,TAX,SAMPLES)


#ALPHA DIVERSITY
#Create composition plot at phylum level
comp_barplot(ps = physeq, tax_level = "Phylum",             
             taxon_renamer = function(x) stringr::str_remove(x, " [ae]t rel."),
             palette = distinct_palette(n = 15, add = "grey90"),
             merge_other = FALSE, bar_outline_colour = "darkgrey"
             )+
          coord_flip() +
           facet_wrap("Stage", nrow = 1, scales = "free") +
           labs(title = "Compositions at Phylum level" ,x = NULL, y = NULL) +
           theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                 panel.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4"),
                 plot.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4"),
                 legend.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4"))
#Create composition plot at Family level
comp_barplot(ps = physeq, tax_level = "Family", n_taxa = 15, other_name = "Other",
             taxon_renamer = function(x) stringr::str_remove(x, " [ae]t rel."),
             palette = distinct_palette(n = 15, add = "grey90"),
             merge_other = FALSE, bar_outline_colour = "darkgrey"
             )+
           coord_flip() +
            facet_wrap("Stage", nrow = 1, scales = "free") +
            labs(title = "Compostions at Family level" , x = NULL, y = NULL) +
            theme(axis.text.y = element_blank(), 
                  axis.ticks.y = element_blank(),
                  panel.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4"),
                  plot.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4"),
                  legend.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4"))             


#Calculate alpha diversity on Phylum rank with Shannon index
a_div_phylum <- ps_calc_diversity(physeq, rank = "Phylum", index = "shannon")
#Merge sample data and Shannon scores together
a_div_phylum <- samdat_tbl(a_div_phylum)
#ggplot to show the alpha diversity on Phylum rank
ggplot(a_div_phylum, aes(x = Stage, y = shannon_Phylum, colour = Stage)) +
  geom_boxplot(width = 0.4) +
  scale_color_manual(values = c("#0099f8", "#2ecc71"))+
  geom_point(position = position_jitter(height = 0.2), alpha = 0.5) +
  theme(axis.ticks = element_blank(),
        axis.line = element_line(colour = "grey50"),
        panel.grid = element_line(color = "#b4aea9"),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(linetype = "dashed"),
        panel.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4"),
        plot.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4"),
        legend.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4"))+
  labs(title="Alpha Diversity Phylum", x="Stage", y="Shannon Index") +
  stat_compare_means(comparisons = list(c("TD", "Autism")),
                     method = "wilcox.test", 
                     label = "p.signif")
#Perform Wilcox sum rank test to check if there is significance
wilcox.test(shannon_Phylum ~ Stage, data = a_div_phylum)


#Repeat alpha diversity calculation but now on family level
a_div_fam <- ps_calc_diversity(physeq, rank = "Family", index = "shannon")
a_div_fam <- samdat_tbl(a_div_fam)
#ggplot to show the alpha diversity on Family rank
ggplot(a_div_fam, aes(x = Stage, y = shannon_Family, colour = Stage)) +
  geom_boxplot(width = 0.4) +
  scale_color_manual(values = c("#0099f8", "#2ecc71"))+
  geom_point(position = position_jitter(height = 0.2), alpha = 0.5) +
  theme(axis.ticks = element_blank(),
        axis.line = element_line(colour = "grey50"),
        panel.grid = element_line(color = "#b4aea9"),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(linetype = "dashed"),
        panel.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4"),
        plot.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4"),
        legend.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4"))+
  labs(title="Alpha Diversity Family", x="Stage", y="Shannon Index") +
  stat_compare_means(comparisons = list(c("TD", "Autism")),
                     method = "wilcox.test", 
                     label = "p.signif")
#Perform Wilcox sum rank test to check if there is significance
wilcox.test(shannon_Family ~ Stage, data = a_div_fam)


#Check the difference of alpha diversity between age categories on phylum rank 
a_div_age <- a_div_phylum
#Create two data frames, first contains only TD samples
a_div_td <- a_div_age %>%
  filter(grepl("TD", Stage))
#Second data frame contains only autistic samples
a_div_asd <- a_div_age %>%
  filter(grepl("Autism", Stage))%>%
  #Because there is 1 sample outlying sample based on age categories, this sample is filtered
  filter(!grepl(13, Age))
#ggplot to show the alpha diversity between age categories in TD
ggplot(a_div_td, aes(x = ageCat, y = shannon_Phylum, colour = ageCat)) +
  geom_boxplot(width = 0.4) +
  scale_color_manual(values = c("#0099f8", "#2ecc71", "#e74c3c"))+
  geom_point(position = position_jitter(height = 0.2), alpha = 0.5) +
  theme(axis.ticks = element_blank(),
        axis.line = element_line(colour = "grey50"),
        panel.grid = element_line(color = "#b4aea9"),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(linetype = "dashed"),
        panel.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4"),
        plot.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4"),
        legend.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4"))+
  labs(title="Alpha Diversity: Age in TD", x="Age", y="Shannon Index") + 
  stat_compare_means(comparisons = list(c("2-3", "4-6"), c("4-6", "7-11"), c("2-3", "7-11")),
                     method = "wilcox.test", 
                     label = "p.signif")
#Subset data to compare age categories with Wilcox test
subset_td1 <- subset(a_div_td, ageCat %in% c("2-3", "4-6"))
wilcox.test(shannon_Phylum ~ ageCat, data = subset_td1)
subset_td2 <- subset(a_div_td, ageCat %in% c("4-6", "7-11"))
wilcox.test(shannon_Phylum ~ ageCat, data = subset_td2)
subset_td3 <- subset(a_div_td, ageCat %in% c("2-3", "7-11"))
wilcox.test(shannon_Phylum ~ ageCat, data = subset_td3)


#ggplot to show the alpha diversity between age categories in Autism
ggplot(a_div_asd, aes(x = ageCat, y = shannon_Phylum, colour = ageCat)) +
  geom_boxplot(width = 0.4) +
  scale_color_manual(values = c("#0099f8", "#2ecc71", "#e74c3c"))+
  geom_point(position = position_jitter(height = 0.2), alpha = 0.5) +
  theme(axis.ticks = element_blank(),
        axis.line = element_line(colour = "grey50"),
        panel.grid = element_line(color = "#b4aea9"),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(linetype = "dashed"),
        panel.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4"),
        plot.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4"),
        legend.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4"))+
  labs(title="Alpha Diversity: Age in Autism", x="Age", y="Shannon Index") +
  stat_compare_means(comparisons = list(c("2-3", "4-6"), c("4-6", "7-11"), c("2-3", "7-11")),
                     method = "wilcox.test", 
                     label = "p.signif")
#Subset data to compare age categories with Wilcox test
subset_asd1 <- subset(a_div_asd, ageCat %in% c("2-3", "4-6"))
wilcox.test(shannon_Phylum ~ ageCat, data = subset_asd1)
subset_asd2 <- subset(a_div_asd, ageCat %in% c("4-6", "7-11"))
wilcox.test(shannon_Phylum ~ ageCat, data = subset_asd2)
subset_asd3 <- subset(a_div_asd, ageCat %in% c("2-3", "7-11"))
wilcox.test(shannon_Phylum ~ ageCat, data = subset_asd3)


#Calculate Firmicutes/Bateriodetes Ratio:
#Create a phyloseq object grouped on Phylum ranks
physeq_phylum <- tax_glom(physeq, taxrank = "Phylum")
#Check the taxonomy table to identify the phylum names
tax_table(physeq_phylum)
# Extract the abundance matrix
abundance_matrix <- otu_table(physeq_phylum)
#Extract taxonomy table
taxonomy <- tax_table(physeq_phylum)
#Get the indices for Firmicutes and Bacteroidetes
firmicutes_idx <- which(taxonomy[, "Phylum"] == "Firmicutes")
bacteroidetes_idx <- which(taxonomy[, "Phylum"] == "Bacteroidetes")
#Extract the abundance of Firmicutes and Bacteroidetes
firmicutes_abundance <- abundance_matrix[firmicutes_idx, ]
bacteroidetes_abundance <- abundance_matrix[bacteroidetes_idx, ]
#Calculate F/B ratio for each sample
fb_ratio <- data.frame(firmicutes_abundance / bacteroidetes_abundance)
rownames(fb_ratio) <- "FBratio"
#Transpose data frame
fb_ratio <- data.frame(t(fb_ratio))
fb_ratio$Stage <- metadata$Stage
#ggplot the F/B ratio between the two groups
ggplot(fb_ratio, aes(x = Stage, y = FBratio, colour = Stage)) +
  geom_boxplot(width = 0.4) +
  scale_color_manual(values = c("#0099f8", "#2ecc71"))+
  geom_point(position = position_jitter(height = 0.2), alpha = 0.5) +
  theme(axis.ticks = element_blank(),
        axis.line = element_line(colour = "grey50"),
        panel.grid = element_line(color = "#b4aea9"),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(linetype = "dashed"),
        panel.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4"),
        plot.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4"),
        legend.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4"))+
  labs(title="Firmicutes/Bacteriodetes ratio", x="Stage", y="F/B ratio")+ 
  stat_compare_means(comparisons = list(c("TD", "Autism")),
                     method = "wilcox.test", 
                     label = "p.signif")
#Wilcox test to prove significance between F/B ratios
wilcox.test(FBratio ~ Stage, data = fb_ratio)


#BETA DIVERSITY
#Filter on prevalence
filterPrev <- rRNA_OTU[,c(3:257)]
#Aggregating by Phylum
filterPrev <- filterPrev%>%
  group_by(Family) %>%
  summarise(across(everything(), sum))
#Count the number of samples (columns) where the OTU (row) is non-zero
prev_abun_stats <- data.frame(Family = filterPrev$Family, 
                         Prevalence = rowSums(filterPrev > 0)/length(filterPrev), 
                         Abundance =  rowSums(filterPrev[,c(2:255)]))

#Add 1 for visualization purposes on log-scale
prev_abun_stats$Abundance <- prev_abun_stats$Abundance + 1
#Visualize cut-offs in plot
prev_abun_stats %>%
  ggplot(aes(Abundance,Prevalence)) +
  geom_vline(xintercept = 800, color = "red", linetype = "dotted", linewidth = 0.8) +
  geom_hline(yintercept = 15 / 100, color = "red", linetype = "dotted", linewidth = 0.8) +
  geom_point(alpha = 0.5) +
  geom_rug(alpha = 0.1) +
  scale_x_log10(
    labels = scales::label_number(), name = "Total Abundance") +
  scale_y_continuous(
    labels = scales::label_percent(), breaks = scales::breaks_pretty(n = 5),
    name = "Prevalence (%)",) +
  theme(axis.ticks = element_blank(),
        axis.line = element_line(colour = "grey50"),
        panel.grid = element_line(color = "#b4aea9"),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(linetype = "dashed"),
        panel.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4"),
        plot.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4"),
        legend.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4")) +
  labs(title="Prevalence & Abundance filtering: Family level")


#Create PCoA plot after filtering on prevalence and abundance on family level
physeq %>%
  tax_filter(min_prevalence = 15/100, tax_level = "Family") %>%
  tax_filter(min_total_abundance = 800, tax_level = "Family") %>%
  dist_calc(dist = "bray") %>%
  ord_calc(method = "PCoA") %>%
  ord_plot(alpha = 0.6, size = 2, color = "Stage") +
  theme(axis.ticks = element_blank(),
        axis.line = element_line(colour = "grey50"),
        panel.grid = element_line(color = "#b4aea9"),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(linetype = "dashed"),
        panel.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4"),
        plot.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4"),
        legend.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4"))+
  coord_fixed(0.8) +
  stat_ellipse(aes(color = Stage)) +
  scale_color_brewer(palette = "Set1") +
  labs(title = "PCoA plot on Bray-Curtis distances")
#PERMANOVA, to check if there is significant difference between the two groups
physeq %>%
  tax_filter(min_prevalence = 15/100, tax_level = "Family") %>%
  tax_filter(min_total_abundance = 800, tax_level = "Family") %>%
  dist_calc(dist = "bray") %>%
  dist_permanova(variables = "Stage", n_perms = 99, seed = 22) %>%
  perm_get()

