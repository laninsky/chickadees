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
setwd("output/")

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

# 9. Printing to screen some parameters of interest to report in manuscript
# Because of a computer issue, had to write in these values from screen shot!
#hist_center <- 9055.73352448063
#mod_center <- 14718.3842749623
#hist_center_min <- 6634.54649994606
#mod_center_min <- 13310.1512531608
#hist_center_max <- 11231.7392024842
#mod_center_max <- 17306.6635855197
#hist_width <- 15671.4033837147
#mod_width <- 8686.17077826922
#hist_width_min <- 9867.47291908409
#mod_width_min <- 5009.68496049775
#hist_width_max <- 17307.3346304032
#mod_width_max <- 17303.7480806549
#hist_pMin <- 0
#mod_pMin <- 0
#hist_pMax <- 0.9998
#mod_pMax <- 0.999

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
## A tibble: 8 × 3
#Parameters        Historical       Modern          
#<chr>             <chr>            <chr>           
#1 Center            9.05573352448063 14.7183842749623
#2 min 95% CI Center 6.63454649994606 13.3101512531608
#3 max 95% CI Center 11.2317392024842 17.3066635855197
#4 Width             15.6714033837147 8.68617077826922
#5 min 95% CI Width  9.86747291908409 5.00968496049775
#6 max 95% CI Width  17.3073346304032 17.3037480806549
#7 pMin              0                0               
#8 pMax              0.9998           0.999   

print("All distances in km")

# Writing out the dataframes in case of needing to reload
write.table(mod_conf_min_max_mean,"../data/hzar_mod_structure.txt",quote=FALSE,col.names=TRUE,row.names=FALSE)
write.table(hist_conf_min_max_mean,"../data/hzar_hist_structure.txt",quote=FALSE,col.names=TRUE,row.names=FALSE)

# 8. Plotting
# Getting the data frames together
# Reading in from file for mod_conf_min_max_mean/hist_conf_min_max_mean because R crashed)
#mod_conf_min_max_mean <- read_table("../data/hzar_mod_structure.txt")
#hist_conf_min_max_mean <- read_table("../data/hzar_hist_structure.txt")

mod_dist <- as_tibble(cbind((temp %>%
                     filter(Sampling_period=="MODERN", Included_in_tess3r=="YES") %>% dplyr::select(Catalog_number)),
                  mod_dist))
names(mod_dist)[2] <- "distance"

hist_dist <- as_tibble(cbind((temp %>%
                               filter(Sampling_period=="SMITHSONIAN", Included_in_tess3r=="YES") %>% dplyr::select(Catalog_number)),
                            hist_dist))
names(hist_dist)[2] <- "distance"

distance_all <- rbind(mod_dist, hist_dist)

temp <- full_join(temp,distance_all,by="Catalog_number")

ggplot() +
  geom_ribbon(data=mod_conf_min_max_mean,mapping=aes(x=dist/1000,ymax = max, ymin = min),fill="grey50",color="black") +
  geom_line(data=mod_conf_min_max_mean,mapping=aes(x=dist/1000,y=mean),color="black", size=2) +
  geom_vline(xintercept = mod_center/1000, linetype="longdash", color = "black", size=3) +
  geom_vline(xintercept = mod_center_min/1000, color = "black", size=2) +
  geom_vline(xintercept = mod_center_max/1000, color = "black", size=2) +
  geom_jitter(data=(temp %>% filter(Sampling_period=="MODERN", Included_in_tess3r=="YES")),aes(x = distance/1000, y = Het, fill=BC_genetic_cluster_assignment),shape=21, color="black", stroke=2,size=15) +
  scale_fill_gradient2(
    low = "#15326C", mid = "#9437FF", high="#CE1B26", midpoint = 0.5, limits=c(-0.01,1.01)
  ) +
  theme_bw(base_size=36) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  theme(axis.title=element_text(size=40,face="bold")) +
  scale_y_continuous(name="Assignment to the BC genetic cluster/\nHeterozygosity (points)", limits=c(-0.01,1.01)) +
  scale_x_continuous(name="Distance along transect (km)", limits=c(-0.01,max(temp$distance,na.rm=TRUE)/1000+0.01)) +
  theme(legend.position = "none")

ggsave(filename="Fig2_modern_het_vs_dist.pdf",plot = last_plot(),width=24,height=12,units="in")

ggplot() +
  geom_ribbon(data=hist_conf_min_max_mean,mapping=aes(x=dist/1000,ymax = max, ymin = min),fill="grey50",color="black") +
  geom_line(data=hist_conf_min_max_mean,mapping=aes(x=dist/1000,y=mean),color="black", size=2) +
  geom_vline(xintercept = hist_center/1000, linetype="longdash", color = "black", size=3) +
  geom_vline(xintercept = hist_center_min/1000, color = "black", size=2) +
  geom_vline(xintercept = hist_center_max/1000, color = "black", size=2) +
  geom_jitter(data=(temp %>% filter(Sampling_period=="SMITHSONIAN", Included_in_tess3r=="YES")),aes(x = distance/1000, y = Het, fill=BC_genetic_cluster_assignment),shape=21, color="black", stroke=2,size=15) +
  scale_fill_gradient2(
    low = "#15326C", mid = "#9437FF", high="#CE1B26", midpoint = 0.5, limits=c(-0.01,1.01)
  ) +
  theme_bw(base_size=36) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  theme(axis.title=element_text(size=40,face="bold")) +
  scale_y_continuous(name="Assignment to the BC genetic cluster/\nHeterozygosity (points)", limits=c(-0.01,1.01)) +
  scale_x_continuous(name="Distance along transect (km)", limits=c(-0.01,max(temp$distance,na.rm=TRUE)/1000+0.01)) +
  theme(legend.position = "none")

ggsave(filename="Fig2_smithsonian_het_vs_dist.pdf",plot = last_plot(),width=24,height=12,units="in")

ggplot() +
  geom_point(data=(temp %>% filter(Sampling_period=="MODERN", Location_code %in% c(1,2,3,4))),aes(x = as.factor(Location_code), y = Het, fill=BC_genetic_cluster_assignment),shape=21, color="black", stroke=2,size=15) +
  scale_fill_gradient2(
    low = "#15326C", mid = "#9437FF", high="#CE1B26", midpoint = 0.5, limits=c(-0.01,1.01)
  ) +
  theme_bw(base_size=36) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  theme(axis.title.y=element_text(size=40,face="bold"),axis.title.x = element_blank()) +
  scale_y_continuous(name="Assignment to the BC genetic cluster/\nHeterozygosity (points)",
                     limits=c(-0.01,1.01), position = "right") +
  theme(legend.position = "none")

ggsave(filename="Fig2_modern_het_vs_dist_BC.pdf",plot = last_plot(),width=6,height=12,units="in")

ggplot() +
  geom_jitter(data=(temp %>% filter(Sampling_period=="MODERN", Location_code==50)),aes(x = as.factor(Location_code), y = Het, fill=BC_genetic_cluster_assignment),shape=21, color="black", stroke=2,size=15) +
  scale_fill_gradient2(
    low = "#15326C", mid = "#9437FF", high="#CE1B26", midpoint = 0.5, limits=c(-0.01,1.01)
  ) +
  theme_bw(base_size=36) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  theme(axis.title.y=element_text(size=40,face="bold"),axis.title.x = element_blank()) +
  scale_y_continuous(name="Assignment to the BC genetic cluster/\nHeterozygosity (points)",
                     limits=c(-0.01,1.01), position = "left") +
  theme(legend.position = "none")

ggsave(filename="Fig2_modern_het_vs_dist_CC.pdf",plot = last_plot(),width=4,height=12,units="in")


ggplot() +
  geom_point(data=(temp %>% filter(Sampling_period=="SMITHSONIAN", Location_code %in% c(1,2))),aes(x = as.factor(Location_code), y = Het, fill=BC_genetic_cluster_assignment),shape=21, color="black", stroke=2,size=15) +
  scale_fill_gradient2(
    low = "#15326C", mid = "#9437FF", high="#CE1B26", midpoint = 0.5, limits=c(-0.01,1.01)
  ) +
  theme_bw(base_size=36) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  theme(axis.title.y=element_text(size=40,face="bold"),axis.title.x = element_blank()) +
  scale_y_continuous(name="Assignment to the BC genetic cluster/\nHeterozygosity (points)",
                     limits=c(-0.01,1.01), position = "right") +
  theme(legend.position = "none")

ggsave(filename="Fig2_hist_het_vs_dist_BC.pdf",plot = last_plot(),width=4.5,height=12,units="in")

ggplot() +
  geom_point(data=(temp %>% filter(Sampling_period=="SMITHSONIAN", Location_code %in% c(20, 21, 22))),aes(x = as.factor(Location_code), y = Het, fill=BC_genetic_cluster_assignment),shape=21, color="black", stroke=2,size=15) +
  scale_fill_gradient2(
    low = "#15326C", mid = "#9437FF", high="#CE1B26", midpoint = 0.5, limits=c(-0.01,1.01)
  ) +
  theme_bw(base_size=36) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  theme(axis.title.y=element_text(size=40,face="bold"),axis.title.x = element_blank()) +
  scale_y_continuous(name="Assignment to the BC genetic cluster/\nHeterozygosity (points)",
                     limits=c(-0.01,1.01), position = "left") +
  theme(legend.position = "none")

ggsave(filename="Fig2_hist_het_vs_dist_CC.pdf",plot = last_plot(),width=5,height=12,units="in")

ggplot() +
  geom_point(data=(temp %>% filter(Sampling_period=="SMITHSONIAN", Location_code %in% c(11,12,15,16,17))),aes(x = as.factor(Location_code), y = Het, fill=BC_genetic_cluster_assignment),shape=21, color="black", stroke=2,size=15) +
  scale_fill_gradient2(
    low = "#15326C", mid = "#9437FF", high="#CE1B26", midpoint = 0.5, limits=c(-0.01,1.01)
  ) +
  theme_bw(base_size=36) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  theme(axis.title.y=element_text(size=40,face="bold"),axis.title.x = element_blank()) +
  scale_y_continuous(name="Assignment to the BC genetic cluster/\nHeterozygosity (points)",
                     limits=c(-0.01,1.01), position = "left") +
  theme(legend.position = "none")

ggsave(filename="Fig2_hist_het_vs_dist_mixed.pdf",plot = last_plot(),width=7,height=12,units="in")

ggplot() +
  geom_point(data=(temp %>% filter(Sampling_period=="MODERN", Included_in_tess3r=="NO", !Location_code %in% c(1,2,3,4,50))),aes(x = as.factor(Location_code), y = Het, fill=BC_genetic_cluster_assignment),shape=21, color="black", stroke=2,size=15) +
  scale_fill_gradient2(
    low = "#15326C", mid = "#9437FF", high="#CE1B26", midpoint = 0.5, limits=c(-0.01,1.01)
  ) +
  theme_bw(base_size=36) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  theme(axis.title.y=element_text(size=40,face="bold"),axis.title.x = element_blank()) +
  scale_y_continuous(name="Assignment to the BC genetic cluster/\nHeterozygosity (points)",
                     limits=c(-0.01,1.01), position = "left") +
  theme(legend.position = "none")

ggsave(filename="Fig2_modern_het_vs_dist_mixed.pdf",plot = last_plot(),width=11,height=12,units="in")


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
#  [1] hzar_0.2-5       foreach_1.5.2    MCMCpack_1.6-3   MASS_7.3-57      coda_0.19-4     
#[6] geosphere_1.5-14 forcats_0.5.1    stringr_1.4.0    dplyr_1.0.9      purrr_0.3.4     
#[11] readr_2.1.2      tidyr_1.2.0      tibble_3.1.7     ggplot2_3.3.6    tidyverse_1.3.1 
#
#loaded via a namespace (and not attached):
#  [1] lubridate_1.8.0    lattice_0.20-45    digest_0.6.29      assertthat_0.2.1   utf8_1.2.2        
#[6] R6_2.5.1           cellranger_1.1.0   backports_1.4.1    MatrixModels_0.5-0 reprex_2.0.1      
#[11] httr_1.4.3         pillar_1.7.0       rlang_1.0.2        readxl_1.4.0       rstudioapi_0.13   
#[16] SparseM_1.81       Matrix_1.4-1       labeling_0.4.2     splines_4.2.1      bit_4.0.4         
#[21] munsell_0.5.0      broom_0.8.0        compiler_4.2.1     modelr_0.1.8       pkgconfig_2.0.3   
#[26] mcmc_0.9-7         tidyselect_1.1.2   codetools_0.2-18   fansi_1.0.3        crayon_1.5.1      
#[31] tzdb_0.3.0         dbplyr_2.2.0       withr_2.5.0        grid_4.2.1         jsonlite_1.8.0    
#[36] gtable_0.3.0       lifecycle_1.0.1    DBI_1.1.3          magrittr_2.0.3     scales_1.2.0      
#[41] cli_3.3.0          stringi_1.7.6      vroom_1.5.7        farver_2.1.0       fs_1.5.2          
#[46] sp_1.5-0           xml2_1.3.3         ellipsis_0.3.2     generics_0.1.2     vctrs_0.4.1       
#[51] iterators_1.0.14   tools_4.2.1        bit64_4.0.5        glue_1.6.2         hms_1.1.1         
#[56] parallel_4.2.1     survival_3.3-1     colorspace_2.0-3   rvest_1.0.2        haven_2.5.0       
#[61] quantreg_5.93  