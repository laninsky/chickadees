# This code creates Fig 2 (spatial interpolation of hybrid zone movement) in the 
# main manuscript of Alexander et al.

# 1. Loading required libraries and scripts
# library(devtools)
# Before the next step need to have openmp installed if on mac. Did this by:
# brew install llvm
# mkdir -p ~/.R
# vi ~/.R/Makevars # Then inserting the following
# C=/usr/local/opt/llvm/bin/clang
# CXX=/usr/local/opt/llvm/bin/clang++
# devtools::install_github("bcm-uga/TESS3_encho_sen")
library(tidyverse)
library(tess3r)
library(ggrepel)
source("http://membres-timc.imag.fr/Olivier.Francois/POPSutilities.R")
library(fields)

# 2. Setwd
setwd("chickadee/output/")

# 3. Reading in data (tab delimited), dropping row
# corresponding to sample 99788 with low coverage
temp <- read_tsv("../data/Table_S1.txt")
temp <- temp %>% filter(Catalog_number!="99788")

# 4. Creating variables with our variables of interest
mod_coords <- as.matrix(temp %>% filter(Sampling_period=="MODERN") %>% filter(Included_in_tess3r=="YES") %>% select(DecimalLongitude,DecimalLatitude))
  
hist_coords <- as.matrix(temp %>% filter(Sampling_period=="SMITHSONIAN") %>% filter(Included_in_tess3r=="YES") %>% select(DecimalLongitude,DecimalLatitude))

mod_q.matrix <- as.matrix(temp %>% filter(Sampling_period=="MODERN") %>% filter(Included_in_tess3r=="YES") %>% select(BC_genetic_cluster_assignment,CC_genetic_cluster_assignment))

hist_q.matrix <-  as.matrix(temp %>% filter(Sampling_period=="SMITHSONIAN") %>% filter(Included_in_tess3r=="YES") %>% select(BC_genetic_cluster_assignment,CC_genetic_cluster_assignment))

grid <- createGrid(min(mod_coords[,1],hist_coords[,1]),max(mod_coords[,1], hist_coords[,1]),min(mod_coords[,2], hist_coords[,2]),max(mod_coords[,2], hist_coords[,2]),2000,2000)

# 5. Creating the base maps of interpolated genome make up
maps(matrix = mod_q.matrix, mod_coords, grid, method = "max", colorGradientsList = list(c("gray95",brewer.pal(9,"Reds")),c("gray95",brewer.pal(9,"Blues"))))
# Manually exported as a *.png 2000*2000 pixels in size, Fig_2_modern_baselayer.png in output folder

maps(matrix = hist_q.matrix, hist_coords, grid, method = "max", colorGradientsList = list(c("gray95",brewer.pal(9,"Reds")),c("gray95",brewer.pal(9,"Blues"))))
# Manually exported as a *.png 2000*2000 pixels in size, Fig_2_historical_baselayer.png in output folder

# 6. Creating site-specific points to manually overlay on base maps
grid <- as.data.frame(grid)

modern_mymarkers <- temp %>% filter(Sampling_period=="MODERN") %>% filter(Included_in_tess3r=="YES") %>% group_by(Location_code, DecimalLongitude, DecimalLatitude) %>% 
  summarise(BCsum=sum(BC_genetic_cluster_assignment),CCsum=sum(CC_genetic_cluster_assignment),max_admixture=max(min(BC_genetic_cluster_assignment,CC_genetic_cluster_assignment))) %>%
  mutate(r=BCsum+CCsum,hybrid_status=ifelse(((BCsum/r)>=0.95 & max_admixture <= 0.05),"BC",ifelse(((CCsum/r)>=0.95 & max_admixture <= 0.05),"CC","Hybrid"))) 

ggplot(grid, aes(x = grid$V1, y = grid$V2)) + geom_point(modern_mymarkers, mapping=aes(x = DecimalLongitude, y = DecimalLatitude,fill = hybrid_status), shape=21,color = "#F2B01E",size=24, stroke = 6)+scale_fill_manual(values=c("#CE1B26", "#15326C","#9437FF")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) + theme(legend.position="none") + geom_label_repel(data=modern_mymarkers, aes(x = DecimalLongitude, y = DecimalLatitude, label = modern_mymarkers$Location_code, color= hybrid_status),fill="#F2B01E",segment.size=0,size=15, force=1.5, max.iter = 40000)+scale_color_manual(values=c("#CE1B26", "#15326C","#9437FF"))+
  xlim(min(grid$V1),max(grid$V1)) +
  ylim(min(grid$V2),max(grid$V2))
# Manually exported as a *.png 2000*2000 pixels in size, Fig_2_modern_overlay.png in output folder

historical_mymarkers <- temp %>% filter(Sampling_period=="SMITHSONIAN") %>% filter(Included_in_tess3r=="YES") %>% group_by(Location_code, DecimalLongitude, DecimalLatitude) %>% 
  summarise(BCsum=sum(BC_genetic_cluster_assignment),CCsum=sum(CC_genetic_cluster_assignment),max_admixture=max(min(BC_genetic_cluster_assignment,CC_genetic_cluster_assignment))) %>%
  mutate(r=BCsum+CCsum,hybrid_status=ifelse(((BCsum/r)>=0.95 & max_admixture <= 0.05),"BC",ifelse(((CCsum/r)>=0.95 & max_admixture <= 0.05),"CC","Hybrid"))) 

ggplot(grid, aes(x = grid$V1, y = grid$V2)) + geom_point(historical_mymarkers, mapping=aes(x = DecimalLongitude, y = DecimalLatitude,fill = hybrid_status), shape=21,color = "#F2B01E",size=24, stroke = 6)+scale_fill_manual(values=c("#CE1B26", "#15326C","#9437FF")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) + theme(legend.position="none") + geom_label_repel(data=historical_mymarkers, aes(x = DecimalLongitude, y = DecimalLatitude, label = historical_mymarkers$Location_code, color= hybrid_status),fill="#F2B01E",segment.size=0,size=15, force=1.5, max.iter = 40000)+scale_color_manual(values=c("#CE1B26", "#15326C","#9437FF"))+
  xlim(min(grid$V1),max(grid$V1)) +
  ylim(min(grid$V2),max(grid$V2))
# Manually exported as a *.png 2000*2000 pixels in size, Fig_2_historical_overlay.png in output folder

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
#  [1] RColorBrewer_1.1-3 fields_14.0        viridis_0.6.2      viridisLite_0.4.0 
#[5] spam_2.9-0         ggrepel_0.9.1      tess3r_1.1.0       forcats_0.5.1     
#[9] stringr_1.4.0      dplyr_1.0.9        purrr_0.3.4        readr_2.1.2       
#[13] tidyr_1.2.0        tibble_3.1.7       ggplot2_3.3.6      tidyverse_1.3.1   
#
#loaded via a namespace (and not attached):
#  [1] Rcpp_1.0.9          lubridate_1.8.0     lattice_0.20-45     digest_0.6.29      
#[5] assertthat_0.2.1    utf8_1.2.2          R6_2.5.1            cellranger_1.1.0   
#[9] backports_1.4.1     reprex_2.0.1        httr_1.4.3          pillar_1.7.0       
#[13] rlang_1.0.2         readxl_1.4.0        rstudioapi_0.13     Matrix_1.4-1       
#[17] labeling_0.4.2      RcppEigen_0.3.3.9.2 bit_4.0.4           munsell_0.5.0      
#[21] broom_0.8.0         compiler_4.2.1      modelr_0.1.8        pkgconfig_2.0.3    
#[25] tidyselect_1.1.2    gridExtra_2.3       fansi_1.0.3         crayon_1.5.1       
#[29] tzdb_0.3.0          dbplyr_2.2.0        withr_2.5.0         grid_4.2.1         
#[33] jsonlite_1.8.0      gtable_0.3.0        lifecycle_1.0.1     DBI_1.1.3          
#[37] magrittr_2.0.3      scales_1.2.0        cli_3.3.0           stringi_1.7.6      
#[41] vroom_1.5.7         farver_2.1.0        fs_1.5.2            xml2_1.3.3         
#[45] ellipsis_0.3.2      generics_0.1.2      vctrs_0.4.1         tools_4.2.1        
#[49] bit64_4.0.5         glue_1.6.2          maps_3.4.0          hms_1.1.1          
#[53] parallel_4.2.1      colorspace_2.0-3    rvest_1.0.2         dotCall64_1.0-1    
#[57] haven_2.5.0       
