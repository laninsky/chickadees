# This code creates PCA components 

# 1. Installing/loading packages
# install.packages("smartsnp")
# install.packages("GGally")
library(smartsnp)
library(tidyverse)
library(GGally)

# 2. Setwd
setwd("output/")

# 3. Loading in STRUCTURE  file
temp <- read_delim("../data/chickadee_singleton_filtered.stru",delim=" ",col_names = FALSE)

# 4. Reformatting the data
# 0 = homozogous "ref" (the first listed allele), 1 = heterozygous, 2 = homozygous "alt"
# Creating an output object
snp_input <- NULL
# Grabbing the sample names
samples <- unique(temp$X1)
# For each of the samples
for (i in samples) {
  print(i)
  # Temporary object to hold the data per sample
  temp_snps <- NULL
  # For each of the SNPs
  for (j in 2:dim(temp)[2]) {
    # Getting the first listed allele
    ref_snp <- temp[1,j]
    # Pulling out the SNPs for the sample of interest
    sample_snps <- temp[which(temp[,1]==i),j]
    # If that sample has missing data
    if (-9 %in% as.matrix(sample_snps)) {
      # Coding using smartsnp missing data terminology
      temp_snps <- c(temp_snps,9)
    } else {
      # If the two alleles are different
      if (length(unique(as.matrix(sample_snps[,1])))==2) {
        temp_snps <- c(temp_snps,1)
      } else {
        if (as.matrix(ref_snp) %in% as.matrix(sample_snps)) {
          temp_snps <- c(temp_snps,0) 
        } else {
          temp_snps <- c(temp_snps,2) 
        }
      }
    }
  }
  snp_input <- cbind(snp_input,temp_snps)
}

write.table(snp_input,"../data/PCA_SNP_input.txt", col.names = FALSE, row.names = FALSE)

# 5. Pulling together files to describe sample groups 
# Some minor tweaking of reference sample names based on tissue 
# number used in structure file vs catalog number in Table S1
samples[which(samples=="3474_")] <- "90612"
samples[which(samples=="6281_")] <- "95776"
samples[which(samples=="9898_")] <- "131638"
samples[which(samples=="7420_")] <- "92269"
samples[which(samples=="7421_")] <- "92270"

# Pulling in Table_S1 and dropping empty final row
tempS1 <- read.delim("../data/Table_S1.txt", sep="\t")
tempS1 <- tempS1[1:165,]
# Also will drop the sample excluded from Structure analyses due to low depth
tempS1 <- tempS1[-(which(tempS1$Catalog_number=="99788")),]

# Matching names in samples to those available in Table_S1
samples <- sub("_.*","",samples)

sample_groups <- NULL
# For each sample
for (i in samples) {
  # Pulling out the sampling period as a group based on Catalog_number and or Tissue_number
  sample_groups <- c(sample_groups,
                     as.matrix(tempS1[which(tempS1$Catalog_number==i),
                                      which(names(tempS1)=="Sampling_period")]))
  sample_groups <- c(sample_groups,
                     as.matrix(tempS1[which(tempS1$Tissue_number==i),
                                      which(names(tempS1)=="Sampling_period")]))
  
}

# 5. Running the PCA
pca_results <- smart_pca(snp_data="../data/PCA_SNP_input.txt",sample_group=sample_groups)

# Amount of variance explained:
pca_results$pca.eigenvalues

# 6. Exploring correlations between PC1 and STRUCTURE assignments
structure_assignments <- NULL
# For each sample
for (i in samples) {
  # Pulling out the structure assignements (to BC cluster) based on Catalog_number and or Tissue_number
  structure_assignments <- c(structure_assignments,
                     as.matrix(tempS1[which(tempS1$Catalog_number==i),
                                      which(names(tempS1)=="BC_genetic_cluster_assignment")]))
  structure_assignments <- c(structure_assignments,
                     as.matrix(tempS1[which(tempS1$Tissue_number==i),
                                      which(names(tempS1)=="BC_genetic_cluster_assignment")]))
  
}

# Pulling out the amount of missing data
missing_data <- colSums(snp_input==9)

data <- as.data.frame(cbind(missing_data,structure_assignments,pca_results$pca.sample_coordinates$PC1,pca_results$pca.sample_coordinates$PC2))
names(data) <- c("missing_data", "structure_assignments", "PC1", "PC2")

# Displaying the correlations
ggpairs(data)
# Saved as a pdf 12 by 8 inches

# Saving the PCA output
saveRDS(pca_results, file = "../data/pca_results.rds")

# Pulling in Table_S1 and dropping empty final row
tempS1 <- read.delim("../data/Table_S1.txt", sep="\t")
tempS1 <- tempS1[1:165,]

pc1scores <- cbind(samples,pca_results$pca.sample_coordinates)

PC_cat_no <- NULL

# Getting the catalog number that corresponds to each row in the PCA results
for (i in 1:dim(pc1scores)[1]) {
  whichcat <- c(which(tempS1$Catalog_number %in% pc1scores$samples[i]),
    which(tempS1$Tissue_number %in% pc1scores$samples[i]))
  tempPCcat_no <- tempS1$Catalog_number[whichcat]
  PC_cat_no <- c(PC_cat_no,tempPCcat_no)
}

pc1scores <- cbind(PC_cat_no, pc1scores)

# Joining this to Table S1
test <- full_join(x=tempS1, y=pc1scores, by = c("Catalog_number" = "PC_cat_no"))
# Dropping superfluous rows
test <- test %>% select(-samples,-Group,-Class,-PC2)

# Writing out Table S1
write.table(test,"../data/Table_S1.txt",quote = FALSE,row.names= FALSE,col.names = TRUE, sep = "\t")

sessionInfo()
#R version 4.1.1 (2021-08-10)
#Platform: x86_64-apple-darwin17.0 (64-bit)
#Running under: macOS Catalina 10.15.7
#
#Matrix products: default
#BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
#LAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylib
#
#locale:
#  [1] en_NZ.UTF-8/en_NZ.UTF-8/en_NZ.UTF-8/C/en_NZ.UTF-8/en_NZ.UTF-8
#
#attached base packages:
#  [1] stats     graphics  grDevices utils     datasets  methods   base     
#
#other attached packages:
#  [1] GGally_2.1.2      forcats_0.5.1     stringr_1.4.0     dplyr_1.0.7       purrr_0.3.4      
#[6] readr_2.1.2       tidyr_1.1.3       tibble_3.1.6      ggplot2_3.3.5     tidyverse_1.3.1  
#[11] smartsnp_1.1.0    data.table_1.14.2
#
#loaded via a namespace (and not attached):
#  [1] Rcpp_1.0.8         lubridate_1.7.10   lattice_0.20-44    prettyunits_1.1.1  assertthat_0.2.1  
#[6] digest_0.6.29      foreach_1.5.2      utf8_1.2.2         RSpectra_0.16-0    plyr_1.8.6        
#[11] R6_2.5.1           cellranger_1.1.0   backports_1.2.1    reprex_2.0.1       RcppZiggurat_0.1.6
#[16] httr_1.4.2         pillar_1.7.0       progress_1.2.2     rlang_1.0.1        readxl_1.3.1      
#[21] rstudioapi_0.13    Matrix_1.3-4       labeling_0.4.2     munsell_0.5.0      broom_0.7.9       
#[26] compiler_4.1.1     modelr_0.1.8       pkgconfig_2.0.3    tidyselect_1.1.1   codetools_0.2-18  
#[31] reshape_0.8.9      fansi_1.0.2        crayon_1.4.2       tzdb_0.2.0         dbplyr_2.1.1      
#[36] withr_2.4.3        grid_4.1.1         jsonlite_1.7.3     gtable_0.3.0       lifecycle_1.0.1   
#[41] DBI_1.1.1          magrittr_2.0.2     scales_1.1.1       Rfast_2.0.6        cli_3.1.1         
#[46] stringi_1.7.6      farver_2.1.0       fs_1.5.0           xml2_1.3.2         ellipsis_0.3.2    
#[51] generics_0.1.0     vctrs_0.3.8        RColorBrewer_1.1-2 iterators_1.0.14   tools_4.1.1       
#[56] glue_1.6.1         hms_1.1.1          parallel_4.1.1     colorspace_2.0-2   rvest_1.0.1       
#[61] haven_2.4.3 
