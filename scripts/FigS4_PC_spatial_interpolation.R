# This code creates Fig. S4 (spatial interpolation of hybrid zone movement using HI) in the 
# supp mats of Alexander et al.

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
library(geosphere)
# Download package from https://cran.r-project.org/src/contrib/Archive/hzar/hzar_0.2-5.tar.gz
# install.packages("MCMCpack")
#install.packages("../../../../../Downloads/hzar_0.2-5.tar.gz", repos = NULL, type="source")
library(hzar)

# 2. Setwd
setwd("chickadees/output/")

# 3. Reading in data (tab delimited), dropping row
# corresponding to sample 99788 with low coverage
temp <- read_tsv("../data/Table_S1.txt")
temp <- temp %>% filter(Catalog_number!="99788")

# Converting PC1 to 0 to 1 scale (to use a q-score)
PC1_range <- max(temp$PC1)-min(temp$PC1)
temp <- temp %>% mutate(PC1_adjusted=((PC1+abs(min(temp$PC1)))/PC1_range),PC1_BC=1-PC1_adjusted)

# 4. Creating variables with our variables of interest
mod_coords <- as.matrix(temp %>% filter(Sampling_period=="MODERN") %>% filter(Included_in_tess3r=="YES") %>% dplyr::select(DecimalLongitude,DecimalLatitude))
  
hist_coords <- as.matrix(temp %>% filter(Sampling_period=="SMITHSONIAN") %>% filter(Included_in_tess3r=="YES") %>% dplyr::select(DecimalLongitude,DecimalLatitude))

mod_q.matrix <- as.matrix(temp %>% filter(Sampling_period=="MODERN") %>% filter(Included_in_tess3r=="YES") %>% dplyr::select(PC1_BC,PC1_adjusted))

hist_q.matrix <-  as.matrix(temp %>% filter(Sampling_period=="SMITHSONIAN") %>% filter(Included_in_tess3r=="YES") %>% dplyr::select(PC1_BC,PC1_adjusted))

grid <- createGrid(min(mod_coords[,1],hist_coords[,1]),max(mod_coords[,1], hist_coords[,1]),min(mod_coords[,2], hist_coords[,2]),max(mod_coords[,2], hist_coords[,2]),2000,2000)

# 5. Creating the base maps of interpolated genome make up
maps(matrix = mod_q.matrix, mod_coords, grid, method = "max", colorGradientsList = list(c("gray95",brewer.pal(9,"Reds")),c("gray95",brewer.pal(9,"Blues"))))
# Manually exported as a *.png 2000*2000 pixels in size, FigS4_PCA_modern_baselayer.png in output folder

maps(matrix = hist_q.matrix, hist_coords, grid, method = "max", colorGradientsList = list(c("gray95",brewer.pal(9,"Reds")),c("gray95",brewer.pal(9,"Blues"))))
# Manually exported as a *.png 2000*2000 pixels in size, FigS4_PCA_historical_baselayer.png in output folder

# 6. Creating site-specific points to manually overlay on base maps
grid <- as.data.frame(grid)

modern_mymarkers <- temp %>% filter(Sampling_period=="MODERN") %>% filter(Included_in_tess3r=="YES") %>% group_by(Location_code, DecimalLongitude, DecimalLatitude) %>% 
  summarise(BCsum=sum(PC1_BC),CCsum=sum(PC1_adjusted),max_admixture=max(min(PC1_BC,PC1_adjusted))) %>%
  mutate(r=BCsum+CCsum,hybrid_status=ifelse(((BCsum/r)>=0.95 & max_admixture <= 0.05),"BC",ifelse(((CCsum/r)>=0.95 & max_admixture <= 0.05),"CC","Hybrid"))) 

# All sampling locations look to be "hybrid" when using the transformed PC1 scores
ggplot(grid, aes(x = grid$V1, y = grid$V2)) + geom_point(modern_mymarkers, mapping=aes(x = DecimalLongitude, y = DecimalLatitude,fill = hybrid_status), shape=21,color = "#F2B01E",size=24, stroke = 6) + scale_fill_manual(values="#9437FF") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) + theme(legend.position="none") + geom_label_repel(data=modern_mymarkers, aes(x = DecimalLongitude, y = DecimalLatitude, label = modern_mymarkers$Location_code, color= hybrid_status),fill="#F2B01E",segment.size=0,size=15, force=1.5, max.iter = 40000)+scale_color_manual(values=c("#9437FF"))+
  xlim(min(grid$V1),max(grid$V1)) +
  ylim(min(grid$V2),max(grid$V2))
# Manually exported as a *.png 2000*2000 pixels in size, FigS4_PCA_modern_overlay.png in output folder

historical_mymarkers <- temp %>% filter(Sampling_period=="SMITHSONIAN") %>% filter(Included_in_tess3r=="YES") %>% group_by(Location_code, DecimalLongitude, DecimalLatitude) %>% 
  summarise(BCsum=sum(PC1_BC),CCsum=sum(PC1_adjusted),max_admixture=max(min(BC_genetic_cluster_assignment,CC_genetic_cluster_assignment))) %>%
  mutate(r=BCsum+CCsum,hybrid_status=ifelse(((BCsum/r)>=0.95 & max_admixture <= 0.05),"BC",ifelse(((CCsum/r)>=0.95 & max_admixture <= 0.05),"CC","Hybrid"))) 

# All sampling locations look to be "hybrid" when using the transformed PC1 scores
ggplot(grid, aes(x = grid$V1, y = grid$V2)) + geom_point(historical_mymarkers, mapping=aes(x = DecimalLongitude, y = DecimalLatitude,fill = hybrid_status), shape=21,color = "#F2B01E",size=24, stroke = 6)+scale_fill_manual(values=c("#9437FF")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) + theme(legend.position="none") + geom_label_repel(data=historical_mymarkers, aes(x = DecimalLongitude, y = DecimalLatitude, label = historical_mymarkers$Location_code, color= hybrid_status),fill="#F2B01E",segment.size=0,size=15, force=1.5, max.iter = 40000)+scale_color_manual(values=c("#9437FF"))+
  xlim(min(grid$V1),max(grid$V1)) +
  ylim(min(grid$V2),max(grid$V2))
# Manually exported as a *.png 2000*2000 pixels in size, FigS4_PCA_historical_overlay.png in output folder

# Now running hzar on things:
# 4. Creating variables with our data of interest
mod_names <- as.matrix(temp %>% filter(Sampling_period=="MODERN") %>% filter(Included_in_tess3r=="YES") %>% dplyr::select(Catalog_number))

hist_names <- as.matrix(temp %>% filter(Sampling_period=="SMITHSONIAN") %>% filter(Included_in_tess3r=="YES") %>% dplyr::select(Catalog_number))

mod_coords <- as.matrix(temp %>% filter(Sampling_period=="MODERN") %>% filter(Included_in_tess3r=="YES") %>% dplyr::select(DecimalLongitude,DecimalLatitude))

hist_coords <- as.matrix(temp %>% filter(Sampling_period=="SMITHSONIAN") %>% filter(Included_in_tess3r=="YES") %>% dplyr::select(DecimalLongitude,DecimalLatitude))

mod_q.matrix <- as.matrix(temp %>% filter(Sampling_period=="MODERN") %>% filter(Included_in_tess3r=="YES") %>% dplyr::select(PC1_BC,PC1_adjusted))

hist_q.matrix <-  as.matrix(temp %>% filter(Sampling_period=="SMITHSONIAN") %>% filter(Included_in_tess3r=="YES") %>% dplyr::select(PC1_BC,PC1_adjusted))

# 5. Calculating distance across the transect
# Obtaining the coordinates to create the approximate SW-NE beginning of transect line
min_long <- min(mod_coords[,1],hist_coords[,1])
max_long <- max(mod_coords[,1], hist_coords[,1])
min_lat <- min(mod_coords[,2], hist_coords[,2])
max_lat <- max(mod_coords[,2], hist_coords[,2])

# Centering the line on the SE corner of the plotting coordinates
long_diff <- max_long-min_long
lat_diff <- max_lat-min_lat

# Longitude/latitude of line as a matrix of 2 columns (first one is longitude, second is latitude) 
transect_center <- matrix(c(min_long,(min_lat-lat_diff),(max_long+long_diff),max_lat),ncol=2,byrow=T)

# Calculating distances 
mod_dist <- (dist2Line(mod_coords, transect_center, distfun=distGeo))[,1]
hist_dist <- (dist2Line(hist_coords, transect_center, distfun=distGeo))[,1]

# Setting minimum distance to zero
min_dist <- min(mod_dist,hist_dist)
mod_dist <- mod_dist-min_dist
hist_dist <- hist_dist-min_dist

# 6. Creating the modern dataframe and modeling
mod_df <- as.data.frame(cbind(mod_names,mod_dist,mod_q.matrix[,1],1))
names(mod_df) <- c("names","distance","frequency","sample_size")
mod_hzar <- hzar.doMolecularData1DPops(mod_df$distance,mod_df$frequency,mod_df$sample_size)

# Building the model, including upper constraints
mod_hzar_cline <- hzar.makeCline1DFreq(data = mod_hzar,scaling = "fixed", tails = "none", direction = NULL)
mod_hzar_cline <- hzar.model.addMaxCenter(mod_hzar_cline, max(mod_dist,hist_dist))
mod_hzar_cline <- hzar.model.addMaxWidth(mod_hzar_cline, max(mod_dist,hist_dist))

# Creating fit request model
mod_fitrequest <- hzar.first.fitRequest.old.ML(mod_hzar_cline,mod_hzar,verbose=TRUE)

# Fitting model
mod_fit <- hzar.doFit(mod_fitrequest)

# Checking fit
plot(mod_fit$mcmcRaw)
hzar.get.ML.cline(mod_fit)$param.free
hzar.getLLCutParam(mod_fit,params = c("center","width"))

# Repeating as often as necessary for stationarity
new_mod_fit <- hzar.next.fitRequest(oldFitRequest = mod_fit)
mod_fit <- hzar.doFit(new_mod_fit)
plot(mod_fit$mcmcRaw)
hzar.get.ML.cline(mod_fit)$param.free
hzar.getLLCutParam(mod_fit,params = c("center","width"))

# After being satisfied stationarity has been reached
# running second chain to verify convergence

# Creating fit request model
mod_fitrequest2 <- hzar.first.fitRequest.old.ML(mod_hzar_cline,mod_hzar,verbose=TRUE)

# Fitting model
mod_fit2 <- hzar.doFit(mod_fitrequest2)

# Checking fit
plot(mod_fit2$mcmcRaw)
hzar.get.ML.cline(mod_fit)$param.free
hzar.get.ML.cline(mod_fit2)$param.free
hzar.getLLCutParam(mod_fit,params = c("center","width"))
hzar.getLLCutParam(mod_fit2,params = c("center","width"))

# Summarizing results (from initial run) to plot
# obtaining the 95% credible interval
mod_conf <- hzar.getCredParamRed(mod_fit)
# recording the parameters of interest
mod_conf_95 <- matrix(NA,ncol=4,nrow=length(mod_conf$clines))
for (i in 1:length(mod_conf$clines)) {
  mod_conf_95[i,1] <- mod_conf$clines[[i]]$param.all$center
  mod_conf_95[i,2] <- mod_conf$clines[[i]]$param.all$width
  mod_conf_95[i,3] <- mod_conf$clines[[i]]$param.all$pMin
  mod_conf_95[i,4] <- mod_conf$clines[[i]]$param.all$pMax
}

# Summarizing the 95% CI for our params of interest
mod_center_min <- min(mod_conf_95[,1])
mod_center_max <- max(mod_conf_95[,1])
mod_width_min <- min(mod_conf_95[,2])
mod_width_max <- max(mod_conf_95[,2])

# Summarizing the best estimate of our parameters
mod_center <- hzar.get.ML.cline(mod_fit)$param.all$center
mod_width <- hzar.get.ML.cline(mod_fit)$param.all$width
mod_pMin <- hzar.get.ML.cline(mod_fit)$param.all$pMin
mod_pMax <- hzar.get.ML.cline(mod_fit)$param.all$pMax

# Generating the 95% CIs for the cline
est_freq <- function (x) {
  pMin + (pMax - pMin) * (1/(1 + exp(-((x - center) * 4/width))))
}

cline_conf <- matrix(NA,ncol=length(mod_conf$clines),nrow=ceiling(max(mod_dist,hist_dist)))

center <- mod_center
width <- mod_width
pMin <- mod_pMin
pMax <- mod_pMax

for (i in 1:length(mod_conf$clines)) {
  center <- mod_conf_95[i,1]
  width <- mod_conf_95[i,2]
  for (j in 1:ceiling(max(mod_dist,hist_dist))) {
    cline_conf[j,i] <- est_freq(j)
  }
}

# Taking the min, max, and mean at each distance
mod_conf_min_max_mean <- matrix(NA,ncol=4,nrow=ceiling(max(mod_dist,hist_dist)))
mod_conf_min_max_mean[,1] <- seq(1,ceiling(max(mod_dist,hist_dist)))

center <- mod_center
width <- mod_width
pMin <- mod_pMin
pMax <- mod_pMax

for (j in 1:ceiling(max(mod_dist,hist_dist))) {
  mod_conf_min_max_mean[j,2] <- min(cline_conf[j,])
  mod_conf_min_max_mean[j,3] <- max(cline_conf[j,])
  mod_conf_min_max_mean[j,4] <- est_freq(j)
}  

# 7. Creating the historical dataframe and modeling
hist_df <- as.data.frame(cbind(hist_names,hist_dist,hist_q.matrix[,1],1))
names(hist_df) <- c("names","distance","frequency","sample_size")
hist_hzar <- hzar.doMolecularData1DPops(hist_df$distance,hist_df$frequency,hist_df$sample_size)

# Building the model, including upper constraints
hist_hzar_cline <- hzar.makeCline1DFreq(data = hist_hzar,scaling = "fixed", tails = "none", direction = NULL)
hist_hzar_cline <- hzar.model.addMaxCenter(hist_hzar_cline, max(mod_dist,hist_dist))
hist_hzar_cline <- hzar.model.addMaxWidth(hist_hzar_cline, max(mod_dist,hist_dist))

# Creating fit request model
hist_fitrequest <- hzar.first.fitRequest.old.ML(hist_hzar_cline,hist_hzar,verbose=TRUE)

# Fitting model
hist_fit <- hzar.doFit(hist_fitrequest)

# Checking fit
plot(hist_fit$mcmcRaw)
hzar.get.ML.cline(hist_fit)$param.free
hzar.getLLCutParam(hist_fit,params = c("center","width"))

# Repeating as often as necessary for stationarity
new_hist_fit <- hzar.next.fitRequest(oldFitRequest = hist_fit)
hist_fit <- hzar.doFit(new_hist_fit)
plot(hist_fit$mcmcRaw)
hzar.get.ML.cline(hist_fit)$param.free
hzar.getLLCutParam(hist_fit,params = c("center","width"))

# After being satisfied stationarity has been reached
# running second chain to verify convergence

# Creating fit request model
hist_fitrequest2 <- hzar.first.fitRequest.old.ML(hist_hzar_cline,hist_hzar,verbose=TRUE)

# Fitting model
hist_fit2 <- hzar.doFit(hist_fitrequest2)

# Checking fit
plot(hist_fit2$mcmcRaw)
hzar.get.ML.cline(hist_fit)$param.free
hzar.get.ML.cline(hist_fit2)$param.free
hzar.getLLCutParam(hist_fit,params = c("center","width"))
hzar.getLLCutParam(hist_fit2,params = c("center","width"))

# Summarizing results (from initial run) to plot
# obtaining the 95% credible interval
hist_conf <- hzar.getCredParamRed(hist_fit)
# recording the parameters of interest
hist_conf_95 <- matrix(NA,ncol=4,nrow=length(hist_conf$clines))
for (i in 1:length(hist_conf$clines)) {
  hist_conf_95[i,1] <- hist_conf$clines[[i]]$param.all$center
  hist_conf_95[i,2] <- hist_conf$clines[[i]]$param.all$width
  hist_conf_95[i,3] <- hist_conf$clines[[i]]$param.all$pMin
  hist_conf_95[i,4] <- hist_conf$clines[[i]]$param.all$pMax
}

# Summarizing the 95% CI for our params of interest
hist_center_min <- min(hist_conf_95[,1])
hist_center_max <- max(hist_conf_95[,1])
hist_width_min <- min(hist_conf_95[,2])
hist_width_max <- max(hist_conf_95[,2])

# Summarizing the best estimate of our parameters
hist_center <- hzar.get.ML.cline(hist_fit)$param.all$center
hist_width <- hzar.get.ML.cline(hist_fit)$param.all$width
hist_pMin <- hzar.get.ML.cline(hist_fit)$param.all$pMin
hist_pMax <- hzar.get.ML.cline(hist_fit)$param.all$pMax

cline_conf <- matrix(NA,ncol=length(hist_conf$clines),nrow=ceiling(max(mod_dist,hist_dist)))

center <- hist_center
width <- hist_width
pMin <- hist_pMin
pMax <- hist_pMax

for (i in 1:length(hist_conf$clines)) {
  center <- hist_conf_95[i,1]
  width <- hist_conf_95[i,2]
  for (j in 1:ceiling(max(mod_dist,hist_dist))) {
    cline_conf[j,i] <- est_freq(j)
  }
}

# Taking the min, max, and mean at each distance
hist_conf_min_max_mean <- matrix(NA,ncol=4,nrow=ceiling(max(mod_dist,hist_dist)))
hist_conf_min_max_mean[,1] <- seq(1,ceiling(max(mod_dist,hist_dist)))

center <- hist_center
width <- hist_width
pMin <- hist_pMin
pMax <- hist_pMax

for (j in 1:ceiling(max(mod_dist,hist_dist))) {
  hist_conf_min_max_mean[j,2] <- min(cline_conf[j,])
  hist_conf_min_max_mean[j,3] <- max(cline_conf[j,])
  hist_conf_min_max_mean[j,4] <- est_freq(j)
}  

# 8. Plotting
# Getting the data frames together
mod_qmat_dist <- as.data.frame(cbind(mod_dist, mod_q.matrix[,1]))
hist_qmat_dist <- as.data.frame(cbind(hist_dist, hist_q.matrix[,1]))
names(mod_qmat_dist) <- c("dist","BC_cluster")
names(hist_qmat_dist) <- c("dist","BC_cluster")
hist_conf_min_max_mean <- as.data.frame(hist_conf_min_max_mean)
mod_conf_min_max_mean <- as.data.frame(mod_conf_min_max_mean)
names(hist_conf_min_max_mean) <- c("dist","min","max","mean")
names(mod_conf_min_max_mean) <- c("dist","min","max","mean")

# Writing out the dataframes in case of needing to reload
write.table(mod_conf_min_max_mean,"../data/hzar_mod_PCA.txt",quote=FALSE,col.names=TRUE,row.names=FALSE)
write.table(hist_conf_min_max_mean,"../data/hzar_hist_PCA.txt",quote=FALSE,col.names=TRUE,row.names=FALSE)

ggplot() +
  geom_ribbon(data=mod_conf_min_max_mean,aes(x=dist/1000,ymax = max, ymin = min),fill="grey50",color="black") +
  geom_line(data=mod_conf_min_max_mean,aes(x=dist/1000,y=mean),color="black", size=5) +
  geom_vline(xintercept = mod_center/1000, linetype="longdash", color = "black", size=7) +
  geom_vline(xintercept = mod_center_min/1000, color = "black", size=5) +
  geom_vline(xintercept = mod_center_max/1000, color = "black", size=5) +
  geom_point(data=mod_qmat_dist,aes(x = dist/1000, y = BC_cluster, fill=BC_cluster),shape=21, color="black", stroke=5,size=30) +
  scale_fill_gradient2(
    low = "#15326C", mid = "#9437FF", high="#CE1B26", midpoint = 0.5
  ) +
  theme_bw(base_size=70) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  theme(axis.title=element_text(size=80,face="bold")) +
  scale_y_continuous(name="Assignment to PC1 (normalized)") +
  scale_x_continuous(name="Distance along transect (km)") +
  theme(legend.position = "none")

# Saved manually as a plot 4000 pixels wide * 2000 pixels wall
# FigS4_PCA_modern_transect.png

ggplot() +
  geom_ribbon(data=hist_conf_min_max_mean,aes(x=dist/1000,ymax = max, ymin = min),fill="grey50",color="black") +
  geom_line(data=hist_conf_min_max_mean,aes(x=dist/1000,y=mean),color="black", size=5) +
  geom_vline(xintercept = hist_center/1000, linetype="longdash", color = "black", size=7) +
  geom_vline(xintercept = hist_center_min/1000, color = "black", size=5) +
  geom_vline(xintercept = hist_center_max/1000, color = "black", size=5) +
  geom_point(data=hist_qmat_dist,aes(x = dist/1000, y = BC_cluster, fill=BC_cluster),shape=21, color="black", stroke=5,size=30) +
  scale_fill_gradient2(
    low = "#15326C", mid = "#9437FF", high="#CE1B26", midpoint = 0.5
  ) +
  theme_bw(base_size=70) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  theme(axis.title=element_text(size=80,face="bold")) +
  scale_y_continuous(name="Assignment to PC1 (normalized)") +
  scale_x_continuous(name="Distance along transect (km)") +
  theme(legend.position = "none")

# Saved manually as a plot 4000 pixels wide * 2000 pixels wall
# FigS4_PCA_historical_transect.png

# 9. Printing to screen some parameters of interest to report in manuscript
results <- as_tibble(rbind(c("Center",hist_center/1000,mod_center/1000),
                           c("min 95% CI Center",hist_center_min/1000,mod_center_min/1000),
                           c("max 95% CI Center",hist_center_max/1000,mod_center_max/1000),
                           c("Width",hist_width/1000,mod_width/1000),
                           c("min 95% CI Width",hist_width_min/1000,mod_width_min/1000),
                           c("max 95% CI Width",hist_width_max/1000,mod_width_max/1000),
                           c("pMin",hist_pMin,mod_pMin),
                           c("pMax",hist_pMax,mod_pMax)))

names(results) <- c("Parameters","Historical","Modern")

results
print("All distances in km")

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
#  [1] hzar_0.2-5         foreach_1.5.2      MCMCpack_1.6-3    
#[4] MASS_7.3-57        coda_0.19-4        geosphere_1.5-14  
#[7] RColorBrewer_1.1-3 fields_14.0        viridis_0.6.2     
#[10] viridisLite_0.4.0  spam_2.9-0         ggrepel_0.9.1     
#[13] tess3r_1.1.0       forcats_0.5.1      stringr_1.4.0     
#[16] dplyr_1.0.9        purrr_0.3.4        readr_2.1.2       
#[19] tidyr_1.2.0        tibble_3.1.7       ggplot2_3.3.6     
#[22] tidyverse_1.3.1   
#
#loaded via a namespace (and not attached):
#  [1] httr_1.4.3          maps_3.4.0          bit64_4.0.5        
#[4] vroom_1.5.7         jsonlite_1.8.0      splines_4.2.1      
#[7] dotCall64_1.0-1     modelr_0.1.8        assertthat_0.2.1   
#[10] sp_1.5-0            cellranger_1.1.0    pillar_1.7.0       
#[13] backports_1.4.1     lattice_0.20-45     quantreg_5.93      
#[16] glue_1.6.2          digest_0.6.29       RcppEigen_0.3.3.9.2
#[19] rvest_1.0.2         colorspace_2.0-3    Matrix_1.4-1       
#[22] pkgconfig_2.0.3     broom_0.8.0         SparseM_1.81       
#[25] haven_2.5.0         scales_1.2.0        tzdb_0.3.0         
#[28] MatrixModels_0.5-0  farver_2.1.0        generics_0.1.2     
#[31] ellipsis_0.3.2      withr_2.5.0         cli_3.3.0          
#[34] survival_3.3-1      magrittr_2.0.3      crayon_1.5.1       
#[37] readxl_1.4.0        mcmc_0.9-7          fs_1.5.2           
#[40] fansi_1.0.3         xml2_1.3.3          tools_4.2.1        
#[43] hms_1.1.1           lifecycle_1.0.1     munsell_0.5.0      
#[46] reprex_2.0.1        compiler_4.2.1      rlang_1.0.2        
#[49] grid_4.2.1          iterators_1.0.14    rstudioapi_0.13    
#[52] labeling_0.4.2      gtable_0.3.0        codetools_0.2-18   
#[55] DBI_1.1.3           R6_2.5.1            gridExtra_2.3      
#[58] lubridate_1.8.0     bit_4.0.4           utf8_1.2.2         
#[61] stringi_1.7.6       parallel_4.2.1      Rcpp_1.0.9         
#[64] vctrs_0.4.1         dbplyr_2.2.0        tidyselect_1.1.2   
