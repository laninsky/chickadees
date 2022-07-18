# This code creates the geographic cline component of Fig. 2 
# in the main manuscript of Alexander et al.

# 1. Loading required libraries and scripts
library(tidyverse)
library(geosphere)
# Download package from https://cran.r-project.org/src/contrib/Archive/hzar/hzar_0.2-5.tar.gz
# install.packages("MCMCpack")
#install.packages("../../../../../Downloads/hzar_0.2-5.tar.gz", repos = NULL, type="source")
library(hzar)

# 2. Setwd
setwd("chickadee/output/")

# 3. Reading in data (tab delimited), dropping row
# corresponding to sample 99788 with low coverage
temp <- read_tsv("../data/Table_S1.txt")
temp <- temp %>% filter(Catalog_number!="99788")

# 4. Creating variables with our data of interest
mod_names <- as.matrix(temp %>% filter(Sampling_period=="MODERN") %>% filter(Included_in_tess3r=="YES") %>% dplyr::select(Catalog_number))

hist_names <- as.matrix(temp %>% filter(Sampling_period=="SMITHSONIAN") %>% filter(Included_in_tess3r=="YES") %>% dplyr::select(Catalog_number))

mod_coords <- as.matrix(temp %>% filter(Sampling_period=="MODERN") %>% filter(Included_in_tess3r=="YES") %>% dplyr::select(DecimalLongitude,DecimalLatitude))

hist_coords <- as.matrix(temp %>% filter(Sampling_period=="SMITHSONIAN") %>% filter(Included_in_tess3r=="YES") %>% dplyr::select(DecimalLongitude,DecimalLatitude))

mod_q.matrix <- as.matrix(temp %>% filter(Sampling_period=="MODERN") %>% filter(Included_in_tess3r=="YES") %>% dplyr::select(BC_genetic_cluster_assignment,CC_genetic_cluster_assignment))

hist_q.matrix <-  as.matrix(temp %>% filter(Sampling_period=="SMITHSONIAN") %>% filter(Included_in_tess3r=="YES") %>% dplyr::select(BC_genetic_cluster_assignment,CC_genetic_cluster_assignment))

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
write.table(mod_conf_min_max_mean,"../data/hzar_mod_structure.txt",quote=FALSE,col.names=TRUE,row.names=FALSE)
write.table(hist_conf_min_max_mean,"../data/hzar_hist_structure.txt",quote=FALSE,col.names=TRUE,row.names=FALSE)

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
  scale_y_continuous(name="Assignment to BC genetic cluster") +
  scale_x_continuous(name="Distance along transect (km)") +
  theme(legend.position = "none")
  
# Saved manually as a plot 4000 pixels wide * 2000 pixels wall
# Fig2_modern_transect.png

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
  scale_y_continuous(name="Assignment to BC genetic cluster") +
  scale_x_continuous(name="Distance along transect (km)") +
  theme(legend.position = "none")

# Saved manually as a plot 4000 pixels wide * 2000 pixels wall
# Fig2_historical_transect.png

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
#  [1] hzar_0.2-5       foreach_1.5.2    MCMCpack_1.6-3   MASS_7.3-57     
#[5] coda_0.19-4      geosphere_1.5-14 forcats_0.5.1    stringr_1.4.0   
#[9] dplyr_1.0.9      purrr_0.3.4      readr_2.1.2      tidyr_1.2.0     
#[13] tibble_3.1.7     ggplot2_3.3.6    tidyverse_1.3.1 
#
#loaded via a namespace (and not attached):
#  [1] lubridate_1.8.0    lattice_0.20-45    digest_0.6.29     
#[4] assertthat_0.2.1   utf8_1.2.2         R6_2.5.1          
#[7] cellranger_1.1.0   backports_1.4.1    MatrixModels_0.5-0
#[10] reprex_2.0.1       httr_1.4.3         pillar_1.7.0      
#[13] rlang_1.0.2        readxl_1.4.0       rstudioapi_0.13   
#[16] SparseM_1.81       Matrix_1.4-1       labeling_0.4.2    
#[19] splines_4.2.1      bit_4.0.4          munsell_0.5.0     
#[22] broom_0.8.0        compiler_4.2.1     modelr_0.1.8      
#[25] pkgconfig_2.0.3    mcmc_0.9-7         tidyselect_1.1.2  
#[28] codetools_0.2-18   fansi_1.0.3        crayon_1.5.1      
#[31] tzdb_0.3.0         dbplyr_2.2.0       withr_2.5.0       
#[34] grid_4.2.1         jsonlite_1.8.0     gtable_0.3.0      
#[37] lifecycle_1.0.1    DBI_1.1.3          magrittr_2.0.3    
#[40] scales_1.2.0       cli_3.3.0          stringi_1.7.6     
#[43] vroom_1.5.7        farver_2.1.0       fs_1.5.2          
#[46] sp_1.5-0           xml2_1.3.3         ellipsis_0.3.2    
#[49] generics_0.1.2     vctrs_0.4.1        iterators_1.0.14  
#[52] tools_4.2.1        bit64_4.0.5        glue_1.6.2        
#[55] hms_1.1.1          parallel_4.2.1     survival_3.3-1    
#[58] colorspace_2.0-3   rvest_1.0.2        haven_2.5.0       
#[61] quantreg_5.93  