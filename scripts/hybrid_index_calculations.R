# This code calculates hybrid indices for the the 
# main manuscript of Alexander et al.

# 1. Loading required libraries and scripts
#library(devtools)
#devtools::install_github("ribailey/gghybrid")
library(gghybrid)
library(tidyverse)
library(GGally)

# 2. Setwd
setwd("chickadee/output/")

# 3. Reading in data
temp <- read.data("../data/chickadee_singleton_filtered.stru", precol.headers=0,
          nprecol=1,NUMLOCI.autoAccept=FALSE,MARKERNAME=0,MISSINGVAL=-9,
          NUMINDS=164,NUMLOCI=8056,PLOIDY=2,POPID=0)

# 4. Prepping the data: note done this "backwards" with BC = 0, CC = 1
temp_cleaned <- data.prep(data=temp$data, loci=temp$loci, sourceAbsent=FALSE,alleles=temp$alleles,
          precols=temp$precols,S0=c("3474_", "9898_"),S1=c("7420_", "7421_"),POPID.name="INDLABEL")

# 5. Estimating the HI
hi_cals <- esth(data.prep.object=temp_cleaned$data.prep,read.data.precols=temp$precols,nitt=3000,burnin=1000) 

# 6. Pulling in Table_S1 and dropping empty final row to match HI up to other data
tempS1 <- read.delim("../data/Table_S1.txt", sep="\t")
tempS1 <- tempS1[1:165,]

samples <- hi_cals$hi$INDLABEL

# Matching names in samples to those available in Table_S1
samples <- sub("_.*","",samples)
samples[164] <- "29898"

hi_scores <- cbind(samples,hi_cals$hi$h_posterior_mode)
hi_scores <- data.frame(hi_scores)

hi_cat_no <- NULL

# Getting the catalog number that corresponds to each row in the PCA results
for (i in 1:dim(hi_scores)[1]) {
  whichcat <- c(which(tempS1$Catalog_number %in% hi_scores$samples[i]),
                which(tempS1$Tissue_number %in% hi_scores$samples[i]))
  temp_hi_cat_no <- tempS1$Catalog_number[whichcat]
  hi_cat_no <- c(hi_cat_no,temp_hi_cat_no)
}

hi_scores <- cbind(hi_cat_no, hi_scores)
names(hi_scores) <- c("hi_cat_no","Hybrid_Index")

# Joining this to Table S1
test <- full_join(x=tempS1, y=hi_scores, by = c("Catalog_number" = "hi_cat_no"))
# Dropping superfluous rows
test <- test %>% select(-samples)

# Writing out Table S1
write.table(test,"../data/Table_S1.txt",quote = FALSE,row.names= FALSE,col.names = TRUE, sep = "\t")

# Examining correlation between STRUCTURE, PC1 and Hybrid_Index
ggpairs(test %>% select(BC_genetic_cluster_assignment,PC1,Hybrid_Index))

sessionInfo()
#R version 4.2.1 (2022-06-23)
#Platform: x86_64-apple-darwin17.0 (64-bit)
#Running under: macOS Catalina 10.15.7
#
#Matrix products: default
#BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
#LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib
#
#locale:
#  [1] en_NZ.UTF-8/en_NZ.UTF-8/en_NZ.UTF-8/C/en_NZ.UTF-8/en_NZ.UTF-8
#
#attached base packages:
#  [1] stats     graphics  grDevices utils     datasets  methods   base     
#
#other attached packages:
#  [1] forcats_0.5.1     stringr_1.4.0     dplyr_1.0.9       purrr_0.3.4       readr_2.1.2      
#[6] tidyr_1.2.0       tibble_3.1.7      ggplot2_3.3.6     tidyverse_1.3.1   gghybrid_1.0.0.0 
#[11] data.table_1.14.2 devtools_2.4.3    usethis_2.1.6    
#
#loaded via a namespace (and not attached):
#  [1] lubridate_1.8.0   prettyunits_1.1.1 ps_1.7.1          assertthat_0.2.1  rprojroot_2.0.3  
#[6] utf8_1.2.2        truncnorm_1.0-8   cellranger_1.1.0  R6_2.5.1          backports_1.4.1  
#[11] reprex_2.0.1      httr_1.4.3        pillar_1.7.0      rlang_1.0.2       curl_4.3.2       
#[16] readxl_1.4.0      rstudioapi_0.13   callr_3.7.0       desc_1.4.1        munsell_0.5.0    
#[21] broom_0.8.0       compiler_4.2.1    modelr_0.1.8      pkgconfig_2.0.3   pkgbuild_1.3.1   
#[26] tidyselect_1.1.2  fansi_1.0.3       crayon_1.5.1      tzdb_0.3.0        dbplyr_2.2.0     
#[31] withr_2.5.0       brio_1.1.3        grid_4.2.1        jsonlite_1.8.0    gtable_0.3.0     
#[36] lifecycle_1.0.1   DBI_1.1.3         magrittr_2.0.3    scales_1.2.0      stringi_1.7.6    
#[41] cli_3.3.0         cachem_1.0.6      fs_1.5.2          remotes_2.4.2     testthat_3.1.4   
#[46] xml2_1.3.3        tester_0.1.7      ellipsis_0.3.2    generics_0.1.2    vctrs_0.4.1      
#[51] tools_4.2.1       glue_1.6.2        hms_1.1.1         processx_3.6.1    pkgload_1.2.4    
#[56] fastmap_1.1.0     colorspace_2.0-3  sessioninfo_1.2.2 rvest_1.0.2       memoise_2.0.1    
#[61] haven_2.5.0 
