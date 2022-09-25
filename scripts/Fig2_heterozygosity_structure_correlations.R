# This code creates aspects of Fig 2 (heterozygosity calculations) 
#  of Alexander et al.

# 1. Loading necessary libraries
library(tidyverse)
library(rstatix)
library(geosphere)

# 2. Setwd
setwd("output/")

# 3. Reading in data (tab delimited), dropping last blank row
temp <- read_tsv("../data/Table_S1.txt")
temp <- temp[1:165,]
structure <- read_table("../data/chickadee_singleton_filtered.stru",col_names = FALSE)

# 4. First need to restrict to SNPs that are diagnostic between black-capped and Carolina
# Combining structure dataset with Table S1 (including structure assigments)
# Preparation to combine files
structure <- structure %>% mutate(Tissue_number=gsub(".*_","",X1),
                              Catalog_number=gsub("_.*","",X1))

structure <- structure %>% mutate(Tissue_number=ifelse(Tissue_number=="",Catalog_number,Tissue_number))
structure$Tissue_number[which(structure$Tissue_number==9898)] <- 29898
structure <- structure %>% mutate(Tissue_number=as.numeric(Tissue_number),
                              Catalog_number=as.numeric(Catalog_number))

# Combining structure file with table S1
combined_SNPs <- full_join(temp,structure,by="Tissue_number") %>% select(-Catalog_number.y)  %>% select(-X1)

# Getting a list of diagnostic SNPs
to_filter <- NULL
for (i in 2:8057) {
  SNP_col <- paste("X",i,sep="")
  print(paste("Up to",SNP_col))
  BC_SNPs <- names(table(combined_SNPs %>% filter(BC_genetic_cluster_assignment >= 0.99) %>% select(UQ(SNP_col))))
  BC_SNPs <- BC_SNPs[BC_SNPs!=-9]
  CC_SNPs <- names(table(combined_SNPs %>% filter(BC_genetic_cluster_assignment <= 0.01) %>% select(UQ(SNP_col))))
  CC_SNPs <- CC_SNPs[CC_SNPs!=-9]
  if (!(any(BC_SNPs %in% CC_SNPs))) {
    to_filter <- c(to_filter,SNP_col)
  }
}

# 5. Filtering dataset to these SNPs and calling sites as het or not
filtered_SNPs <- combined_SNPs %>% select(Catalog_number.x,UQ(to_filter))

# Need to filter out the low coverage sample not present in the Structure dataset
filtered_SNPs <- filtered_SNPs %>% filter(!(Catalog_number.x=="99788"))
  
struc_names <- c(unique(filtered_SNPs$Catalog_number.x))

het_record <- matrix(NA,ncol = (dim(filtered_SNPs)[2]-1),nrow=length(struc_names))

for (i in 1:length(struc_names)) {
  # For each of the samples
  sample <- struc_names[i]
  print(paste("Up to",sample,sep=" "))
  # For each of the SNPs
  for (j in 2:dim(filtered_SNPs)[2]) {
    # If neither of the alleles are -9 (our missing data value)
    if (!(any(filtered_SNPs[(which(filtered_SNPs$Catalog_number.x==sample)),j]==-9))) {
      # If only one allele state is present (e.g. homozygous)
      if(dim(unique(filtered_SNPs[(which(filtered_SNPs$Catalog_number.x==sample)),j]))[1]==1) {
        het_record[i,(j-1)] <- 0
      } else {
        het_record[i,(j-1)] <- 1
      }
    }
  }
}

het_est <- as_tibble(cbind(struc_names,rowMeans(het_record,na.rm=TRUE)))
names(het_est) <- c("Catalog_number","Het")

combined_het <- full_join(temp,het_est,by="Catalog_number")

# 7. Estimates of prevalence of hybrid individuals
# Maximum heterozygosity estimate among "unadmixed" individuals
combined_het %>% 
  filter(BC_genetic_cluster_assignment <= 0.01 | BC_genetic_cluster_assignment >= 0.99) %>% 
  select(Het) %>% max()
# 0.07692308

# Proportions of individuals with greater heterozygosity than that estimate
# for each study period. 92 Modern samples,68 Smithsonian
combined_het %>% filter(Sampling_period=="SMITHSONIAN") %>% filter(Het>0.07692308)
13/68
combined_het %>% filter(Sampling_period=="MODERN") %>% filter(Het>0.07692308)
23/92

# Lower levels of heterozygosity could also be evidence of more complex backcrossing
combined_het %>% filter(Het>0.07692308) %>% t_test(Het ~ Sampling_period, var.equal = F)
combined_het  %>% t_test(Het ~ Sampling_period, var.equal = F)

# Suggestive overall proportions of hybrids within the dataset have not shifted 
# much over time, however caution needed because full extent of hybrid zone has
# not been sampled in the modern dataset. However, also consistent with strong
# selection against hybrids previously observed

# Writing out Table S1 with Het calculations
write.table((combined_het  %>% select(-distance)),"../data/Table_S1.txt",quote = FALSE,row.names= FALSE,col.names = TRUE, sep = "\t")

sessionInfo()
#R version 4.2.1 (2022-06-23)
#Platform: x86_64-apple-darwin17.0 (64-bit)
#Running under: macOS Catalina 10.15.7

#Matrix products: default
#BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
#LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib

#locale:
#  [1] en_NZ.UTF-8/en_NZ.UTF-8/en_NZ.UTF-8/C/en_NZ.UTF-8/en_NZ.UTF-8

#attached base packages:
#  [1] stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#  [1] rstatix_0.7.0   forcats_0.5.1   stringr_1.4.0   dplyr_1.0.9     purrr_0.3.4    
#[6] readr_2.1.2     tidyr_1.2.0     tibble_3.1.7    ggplot2_3.3.6   tidyverse_1.3.1

#loaded via a namespace (and not attached):
#  [1] tidyselect_1.1.2 haven_2.5.0      carData_3.0-5    colorspace_2.0-3 vctrs_0.4.1     
#[6] generics_0.1.2   utf8_1.2.2       rlang_1.0.2      pillar_1.7.0     glue_1.6.2      
#[11] withr_2.5.0      DBI_1.1.3        bit64_4.0.5      dbplyr_2.2.0     modelr_0.1.8    
#[16] readxl_1.4.0     lifecycle_1.0.1  munsell_0.5.0    gtable_0.3.0     cellranger_1.1.0
#[21] rvest_1.0.2      labeling_0.4.2   tzdb_0.3.0       parallel_4.2.1   fansi_1.0.3     
#[26] broom_0.8.0      scales_1.2.0     backports_1.4.1  vroom_1.5.7      jsonlite_1.8.0  
#[31] abind_1.4-5      farver_2.1.0     fs_1.5.2         bit_4.0.4        hms_1.1.1       
#[36] digest_0.6.29    stringi_1.7.6    grid_4.2.1       cli_3.3.0        tools_4.2.1     
#[41] magrittr_2.0.3   crayon_1.5.1     car_3.1-0        pkgconfig_2.0.3  ellipsis_0.3.2  
#[46] xml2_1.3.3       reprex_2.0.1     lubridate_1.8.0  assertthat_0.2.1 httr_1.4.3      
#[51] rstudioapi_0.13  R6_2.5.1         compiler_4.2.1