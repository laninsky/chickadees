# This code creates Fig 1 (individual bird structure assignments and coloring 
# site labels by “hybrid status” of Alexander et al.

# 1. Loading necessary libraries
library(tidyverse)
library(ggmap)
library(ggrepel)
library(ggsn)

# 2. Setwd
setwd("chickadee/output/")

# 3. Reading in data (tab delimited), dropping last blank row
temp <- read_tsv("../data/Table_S1.txt")
temp <- temp[1:165,]

# 4. Creating a new table with the variables we are interested in
reduced_table <- temp %>% filter(Sampling_period=="MODERN" | Sampling_period=="SMITHSONIAN") %>% select(Catalog_number, BC_genetic_cluster_assignment, CC_genetic_cluster_assignment, DecimalLongitude, DecimalLatitude, Sampling_period, Location_code,Included_in_tess3r,PC1)

reduced_table # should print out something similar to following
## A tibble: 160 × 9
#Catalog_number BC_genetic_cluster_a… CC_genetic_clus… DecimalLongitude DecimalLatitude
#<dbl>                 <dbl>            <dbl>            <dbl>           <dbl>
#  1         132041                 0.001            0.999            -94.1            38.0
#2         132042                 0.001            0.999            -94.1            38.0
#3         132043                 0.001            0.999            -94.1            38.0
#4         132044                 0.001            0.999            -94.1            38.0
#5         132045                 0.001            0.999            -94.1            38.0
#6         132046                 0.001            0.999            -94.1            38.0
#7         132047                 0.001            0.999            -94.1            38.0
#8         132048                 0.001            0.999            -93.7            37.9
#9         132049                 0.001            0.999            -93.7            37.9
#10         132050                 0.001            0.999            -93.7            37.9
## … with 150 more rows, and 4 more variables: Sampling_period <chr>,
##   Location_code <dbl>, Included_in_tess3r <chr>, PC1 <dbl>

# 5. Filtering data for each sampling period
modern <- reduced_table %>% filter(Sampling_period=="MODERN") %>% group_by(DecimalLongitude,DecimalLatitude,Location_code,Included_in_tess3r) %>% 
  summarise(BCsum=sum(BC_genetic_cluster_assignment),CCsum=sum(CC_genetic_cluster_assignment)) %>%
  mutate(r=BCsum+CCsum,status=ifelse((BCsum/r)>=0.95,"BC",ifelse((CCsum/r)>=0.95,"CC","Hybrid"))) 

historical <- reduced_table %>% filter(Sampling_period=="SMITHSONIAN") %>% group_by(DecimalLongitude,DecimalLatitude,Location_code,Included_in_tess3r) %>% 
  summarise(BCsum=sum(BC_genetic_cluster_assignment,na.rm=TRUE),CCsum=sum(CC_genetic_cluster_assignment,na.rm=TRUE)) %>%
  mutate(r=BCsum+CCsum,status=ifelse((BCsum/r)>=0.95,"BC",ifelse((CCsum/r)>=0.95,"CC","Hybrid"))) 

# 6. Making a box containing all sites in both modern and historical samples
# along with some buffer
latbound <- max(reduced_table$DecimalLatitude)-min(reduced_table$DecimalLatitude)+0.2
longbound <- max(reduced_table$DecimalLongitude)-min(reduced_table$DecimalLongitude)+0.2
if (latbound > longbound) {
  midlong <- (max(reduced_table$DecimalLongitude)+min(reduced_table$DecimalLongitude))/2
  minlong <- midlong-longbound/2
  maxlong <- midlong+longbound/2
  minlat <- min(reduced_table$DecimalLatitude)-0.1
  maxlat <- max(reduced_table$DecimalLatitude)+0.1
} else {
  midlat <- (max(reduced_table$DecimalLatitude)+min(reduced_table$DecimalLatitude))/2
  minlat <- midlat-longbound/2
  maxlat <- midlat+longbound/2
  minlong <- min(reduced_table$DecimalLongitude)-0.1
  maxlong <- max(reduced_table$DecimalLongitude)+0.1
}

# Taking the extent of the plot out to the Kansas border if it doesn't
# already overlap with it.
if (minlong > -94.6333333) {
  original_extent <- maxlong - minlong
  minlong <- -94.6333333
  new_extent <- maxlong - minlong
  minlat <- minlat - (new_extent-original_extent)/2
  maxlat <- maxlat + (new_extent-original_extent)/2
  sbbox <- make_bbox(lon=c(minlong,maxlong), lat=c(minlat,maxlat),f=0)
} else {
  sbbox <- make_bbox(lon=c(minlong,maxlong), lat=c(minlat,maxlat),f=0)
}

# 7. Using the sbbox object to retrieve a map covering the sample sites
sq_map <- get_map(location = sbbox, maptype = "terrain-background", source = "stamen", crop=TRUE)

# 8. Creating maps of sample locations by "hybrid status"
modernsamplelocations <- ggmap(sq_map) + 
  geom_vline(xintercept=-94.63333) +
  geom_point(data = modern, mapping = aes(x = DecimalLongitude, y = DecimalLatitude,fill = status), shape=21,color = "black",size=12)+
  scale_x_continuous(expand=c(0,0),position = "top") +
  scale_fill_manual(values=c("#CE1B26", "#15326C","#9437FF")) +
  coord_fixed(ratio=1) + xlab("Longitude") + ylab("Latitude") + 
  theme_bw(base_size = 20) +
  theme(legend.position="none",panel.border=element_rect(fill = NA)) + 
  theme(axis.title=element_text(size=28,face="bold"))

modernsamplelocations + geom_label_repel(modern,mapping=aes(fontface="bold",x=DecimalLongitude,y=DecimalLatitude,label=Location_code,fill=Included_in_tess3r,color=status),size=8, force=20, max.iter = 40000, segment.color="black") + scale_fill_manual(values=c("#CE1B26", "#15326C","#9437FF","#FFFFFF", "#F2B01F")) +
  scale_color_manual(values=c("#CE1B26", "#15326C","#9437FF")) +
  theme(plot.margin=unit(c(1,1,1,1),"cm"))

ggsave(filename="Fig1_modern_aggregated_by_site.pdf",plot = last_plot(),width=12.804,height=10.351,units="in")

historicalsamplelocations <- ggmap(sq_map)  + geom_vline(xintercept=-94.63333) + 
  geom_point(data = historical, mapping = aes(x = DecimalLongitude, y = DecimalLatitude,fill = status), shape=21,color = "black",size=12) +
  scale_x_continuous(expand=c(0,0),position = "bottom") +
  scale_fill_manual(values=c("#CE1B26", "#15326C","#9437FF")) +
  coord_fixed(ratio=1) + xlab("Longitude") + ylab("Latitude") + 
  theme_bw(base_size = 20) +
          theme(legend.position="none",panel.border=element_rect(fill = NA)) + 
  theme(axis.title=element_text(size=28,face="bold")) 

historicalsamplelocations + geom_label_repel(historical,mapping=aes(fontface="bold",x=DecimalLongitude,y=DecimalLatitude,label=Location_code,fill=Included_in_tess3r,color=status),size=8, force=10, max.iter = 40000, segment.color="black") + scale_fill_manual(values=c("#CE1B26", "#15326C","#9437FF","#FFFFFF", "#F2B01F")) +
scale_color_manual(values=c("#CE1B26", "#15326C","#9437FF")) +
  theme(plot.margin=unit(c(1,1,1,1),"cm"))

ggsave(filename="Fig1_historical_aggregated_by_site.pdf",plot = last_plot(),width=12.804,height=10.351,units="in")

# 9. Creating the structure plots ordered by location
# Modern by longitude, export as 1000 pixels wide: Fig1_structure_bars_long_modern.png
modernindividualnames <- temp %>% arrange(desc(BC_genetic_cluster_assignment)) %>% arrange(Location_code) %>% arrange(DecimalLongitude) %>% filter(Sampling_period=="MODERN")

modernindividual <- temp %>% filter(Sampling_period=="MODERN") %>% gather(cluster_assignment,assignment_value,c(BC_genetic_cluster_assignment,CC_genetic_cluster_assignment)) %>% arrange(desc(cluster_assignment)) %>% arrange(Location_code) %>% arrange(DecimalLongitude)

ggplot(modernindividual, aes(fill=cluster_assignment,y=assignment_value,x=as.factor(Catalog_number))) +
  geom_bar(stat="identity",color="black",width=1) + theme (legend.position="none", axis.text = element_blank(),axis.title=element_blank(), axis.ticks = element_blank()) + scale_y_continuous(limits=c(0,1),expand = c(0, 0)) +
  theme(aspect.ratio = 1/4) +
  scale_x_discrete(limits=modernindividualnames$Catalog_number) +
  scale_fill_manual(values = (c("#CE1B26","#15326C")))

# Historical by longitude, export as 1000 pixels wide: Fig1_structure_bars_long_historical.png
historicalindividualnames <- temp %>% filter(!is.na(BC_genetic_cluster_assignment)) %>% arrange(desc(BC_genetic_cluster_assignment)) %>% arrange(Location_code) %>% arrange(DecimalLongitude) %>% filter(Sampling_period=="SMITHSONIAN")

historicalindividual <- temp %>% filter(Sampling_period=="SMITHSONIAN") %>% filter(!is.na(BC_genetic_cluster_assignment)) %>% gather(cluster_assignment,assignment_value,c(BC_genetic_cluster_assignment,CC_genetic_cluster_assignment)) %>% arrange(desc(cluster_assignment)) %>% arrange(Location_code) %>% arrange(DecimalLongitude)

ggplot(historicalindividual, aes(fill=cluster_assignment,y=assignment_value,x=as.factor(Catalog_number))) +
  geom_bar(stat="identity",color="black",width=1) + theme (legend.position="none", axis.text = element_blank(),axis.title=element_blank(), axis.ticks = element_blank()) + scale_y_continuous(limits=c(0,1),expand = c(0, 0)) +
  theme(aspect.ratio = 1/4) +
  scale_x_discrete(limits=historicalindividualnames$Catalog_number) +
  scale_fill_manual(values = (c("#CE1B26","#15326C")))

# Because latitude largely recapitulates the longitudinal pattern (given the orientation)
# of the plot, only exporting based on longitude
# Modern by latitude, export as 1000 pixels height: Fig_S1_modern_lat.png
#modernindividualnames <- temp %>% arrange(desc(BC_genetic_cluster_assignment)) %>% arrange(Location_code) %>% arrange(DecimalLatitude) %>% filter(Sampling_period=="MODERN")

#modernindividual <- temp %>% filter(Sampling_period=="MODERN") %>% gather(cluster_assignment,assignment_value,c(BC_genetic_cluster_assignment,CC_genetic_cluster_assignment)) %>% arrange(desc(cluster_assignment)) %>% arrange(Location_code) %>% arrange(DecimalLatitude)

#ggplot(modernindividual, aes(fill=cluster_assignment,y=assignment_value,x=as.factor(Catalog_number))) +
#  geom_bar(stat="identity",color="black",width=1) + theme (legend.position="none", axis.text = element_blank(),axis.title=element_blank(), axis.ticks = element_blank()) + scale_y_continuous(limits=c(0,1),expand = c(0, 0)) +
#  theme(aspect.ratio = 4/1) +
#  scale_x_discrete(limits=modernindividualnames$Catalog_number) +
#  scale_fill_manual(values = (c("#CE1B26","#15326C"))) +
#  coord_flip()

# Modern by latitude, export as 1000 pixels height: Fig_S1_historical_lat.png
#historicalindividualnames <- temp %>% arrange(desc(BC_genetic_cluster_assignment)) %>% arrange(Location_code) %>% arrange(DecimalLatitude) %>% filter(Sampling_period=="SMITHSONIAN")

#historicalindividual <- temp %>% filter(Sampling_period=="SMITHSONIAN") %>% gather(cluster_assignment,assignment_value,c(BC_genetic_cluster_assignment,CC_genetic_cluster_assignment)) %>% arrange(desc(cluster_assignment)) %>% arrange(Location_code) %>% arrange(DecimalLatitude)

#ggplot(historicalindividual, aes(fill=cluster_assignment,y=assignment_value,x=as.factor(Catalog_number))) +
#  geom_bar(stat="identity",color="black",width=1) + theme (legend.position="none", axis.text = element_blank(),axis.title=element_blank(), axis.ticks = element_blank()) + scale_y_continuous(limits=c(0,1),expand = c(0, 0)) +
#  theme(aspect.ratio = 4/1) +
#  scale_x_discrete(limits=historicalindividualnames$Catalog_number) +
#  scale_fill_manual(values = (c("#CE1B26","#15326C"))) +
#  coord_flip()

# 10. Generating site labels, first for modern by longitude (the x-axis for the plot)
modern <- modern %>% arrange(Location_code) %>% arrange(DecimalLongitude)

# replacing "status" with the appropriate color
fill_codes <- as.matrix(modern$status)[,1]
fill_codes <- gsub("Hybrid","#9437FF",gsub("CC","#15326C",gsub("BC","#CE1B26",fill_codes)))

# creating a vector for the position of our sampling site text labels
labelx <- modern$r[length(modern$r)]/2
for (i in (dim(modern)[1]-1):1) {
  labelx <- c(labelx, (labelx[length(labelx)]+(modern$r[i+1]/2)+(modern$r[i]/2)))
}

# creating a vector for the y-position of the labels to stagger them
# and avoid overlap
labely <- rep(4,length(modern$r))
labely[seq(2,length(modern$r),4)] <- 1
labely[seq(3,length(modern$r),4)] <- 3
labely[seq(4,length(modern$r),4)] <- 0

# building the base colored plot
modernlong <- ggplot(modern,aes(y=r,x=2)) + geom_bar(stat="identity",color="black",aes(fill=factor(modern$Location_code,levels=as.numeric(modern$Location_code)))) +
  theme_classic() +
  theme(panel.border = element_blank(),axis.line=element_blank()) +
  scale_x_continuous(limits=c(-1,5)) +
  scale_fill_manual(values = c(fill_codes)) +
  theme (legend.position="none", axis.text = element_blank(),axis.title=element_blank(), axis.ticks = element_blank())  +
  theme(aspect.ratio = 1/5) +
  coord_flip() 

# Figuring out what the colors should be for fill and outline of the 
# first site label
if (modern$Included_in_tess3r[length(modern$Included_in_tess3r)]=="NO") {
  fillcolor <- "#FFFFFF"
} else {
  fillcolor <- "#F2B01F"
}
if (modern$status[length(modern$Included_in_tess3r)]=="Hybrid") {
  colorcolor <- "#9437FF"
} else {
  if (modern$status[length(modern$Included_in_tess3r)]=="CC") {
    colorcolor <- "#15326C"
  }  else {
    colorcolor <- "#CE1B26"
  }
}

# creating a vector to put in some vlines separating the labels
vlinepos <- c(0,modern$r[length(modern$r)])
for (i in (dim(modern)[1]-1):1) {
  vlinepos <- c(vlinepos, (vlinepos[length(vlinepos)]+(modern$r[i])))
}

# and adding this annotation on to the baseplot
modernlong <- modernlong + geom_hline(yintercept = vlinepos, linetype="dashed")

# and adding text label annotation on to the baseplot for left-most site
modernlong <- modernlong + annotate("label",y=labelx[1],x=labely[1],fontface="bold",label=modern$Location_code[length(modern$Included_in_tess3r)],fill=fillcolor,color=colorcolor,size=12)

# and now doing this for all the remaining labels
for (i in 2:length(labelx)) {
  if (modern$Included_in_tess3r[length(modern$Included_in_tess3r)-i+1]=="NO") {
    fillcolor <- "#FFFFFF"
  } else {
    fillcolor <- "#F2B01F"
  }
  if (modern$status[length(modern$Included_in_tess3r)-i+1]=="Hybrid") {
    colorcolor <- "#9437FF"
  } else {
    if (modern$status[length(modern$Included_in_tess3r)-i+1]=="CC") {
      colorcolor <- "#15326C"
    }  else {
    colorcolor <- "#CE1B26"
    }
  }
  modernlong <- modernlong + annotate("label",y=labelx[i],x=labely[i],fontface="bold",label=modern$Location_code[length(modern$Included_in_tess3r)-i+1],fill=fillcolor,color=colorcolor,size=12)
}  

# Exporting as 1500 pixels width, Fig1_modern_long_sites.png
modernlong + scale_y_reverse() 

# 11. Generating site labels, next historical by longitude (the x-axis for the plot)
historical <- historical %>% arrange(Location_code) %>% arrange(DecimalLongitude)

# replacing "status" with the appropriate color
fill_codes <- as.matrix(historical$status)[,1]
fill_codes <- gsub("Hybrid","#9437FF",gsub("CC","#15326C",gsub("BC","#CE1B26",fill_codes)))

# creating a vector for the position of our sampling site text labels
labelx <- historical$r[length(historical$r)]/2
for (i in (dim(historical)[1]-1):1) {
  labelx <- c(labelx, (labelx[length(labelx)]+(historical$r[i+1]/2)+(historical$r[i]/2)))
}

# creating a vector for the y-position of the labels to stagger them
# and avoid overlap
labely <- rep(4,length(historical$r))
labely[seq(2,length(historical$r),4)] <- 1
labely[seq(3,length(historical$r),4)] <- 3
labely[seq(4,length(historical$r),4)] <- 0

# building the base colored plot
historicallong <- ggplot(historical,aes(y=r,x=2)) + geom_bar(stat="identity",color="black",aes(fill=factor(historical$Location_code,levels=as.numeric(historical$Location_code)))) +
  theme_classic() +
  theme(panel.border = element_blank(),axis.line=element_blank()) +
  scale_x_continuous(limits=c(-1,5)) +
  scale_fill_manual(values = c(fill_codes)) +
  theme (legend.position="none", axis.text = element_blank(),axis.title=element_blank(), axis.ticks = element_blank())  +
  theme(aspect.ratio = 1/5) +
  coord_flip() 

# Figuring out what the colors should be for fill and outline of the 
# first site label
if (historical$Included_in_tess3r[length(historical$Included_in_tess3r)]=="NO") {
  fillcolor <- "#FFFFFF"
} else {
  fillcolor <- "#F2B01F"
}
if (historical$status[length(historical$Included_in_tess3r)]=="Hybrid") {
  colorcolor <- "#9437FF"
} else {
  if (historical$status[length(historical$Included_in_tess3r)]=="CC") {
    colorcolor <- "#15326C"
  }  else {
    colorcolor <- "#CE1B26"
  }
}

# creating a vector to put in some vlines separating the labels
vlinepos <- c(0,historical$r[length(historical$r)])
for (i in (dim(historical)[1]-1):1) {
  vlinepos <- c(vlinepos, (vlinepos[length(vlinepos)]+(historical$r[i])))
}

# and adding this annotation on to the baseplot
historicallong <- historicallong + geom_hline(yintercept = vlinepos, linetype="dashed")

# and adding text label annotation on to the baseplot for left-most site
historicallong <- historicallong + annotate("label",y=labelx[1],x=labely[1],fontface="bold",label=historical$Location_code[length(historical$Included_in_tess3r)],fill=fillcolor,color=colorcolor,size=12)

# and now doing this for all the remaining labels
for (i in 2:length(labelx)) {
  if (historical$Included_in_tess3r[length(historical$Included_in_tess3r)-i+1]=="NO") {
    fillcolor <- "#FFFFFF"
  } else {
    fillcolor <- "#F2B01F"
  }
  if (historical$status[length(historical$Included_in_tess3r)-i+1]=="Hybrid") {
    colorcolor <- "#9437FF"
  } else {
    if (historical$status[length(historical$Included_in_tess3r)-i+1]=="CC") {
      colorcolor <- "#15326C"
    }  else {
      colorcolor <- "#CE1B26"
    }
  }
  historicallong <- historicallong + annotate("label",y=labelx[i],x=labely[i],fontface="bold",label=historical$Location_code[length(historical$Included_in_tess3r)-i+1],fill=fillcolor,color=colorcolor,size=12)
}  

# Exporting as 1500 pixels width, Fig1_historical_long_sites.png
historicallong + scale_y_reverse() 

# Because latitude largely recapitulates the longitudinal pattern (given the orientation)
# of the plot, only exporting based on longitude
# 12. Generating site labels, next historical by latitude (the y-axis for the plot)
#historical <- historical %>% arrange(Location_code) %>% arrange(desc(DecimalLatitude))

# replacing "status" with the appropriate color
#fill_codes <- as.matrix(historical$status)[,1]
#fill_codes <- gsub("Hybrid","#9437FF",gsub("CC","#15326C",gsub("BC","#CE1B26",fill_codes)))

# creating a vector for the position of our sampling site text labels
#labelx <- historical$r[length(historical$r)]/2
#for (i in (dim(historical)[1]-1):1) {
#  labelx <- c(labelx, (labelx[length(labelx)]+(historical$r[i+1]/2)+(historical$r[i]/2)))
#}

# creating a vector for the y-position of the labels to stagger them
# and avoid overlap
#labely <- rep(4,length(historical$r))
#labely[seq(2,length(historical$r),4)] <- 1
#labely[seq(3,length(historical$r),4)] <- 3
#labely[seq(4,length(historical$r),4)] <- 0

# building the base colored plot
#historicallat <- ggplot(historical,aes(y=r,x=2)) + geom_bar(stat="identity",color="black",aes(fill=factor(historical$Location_code,levels=as.numeric(historical$Location_code)))) +
#  theme_classic() +
#  theme(panel.border = element_blank(),axis.line=element_blank()) +
#  scale_x_continuous(limits=c(-1,5)) +
#  scale_fill_manual(values = c(fill_codes)) +
#  theme (legend.position="none", axis.text = element_blank(),axis.title=element_blank(), axis.ticks = element_blank())  +
#  theme(aspect.ratio = 5/1)

# Figuring out what the colors should be for fill and outline of the 
# first site label
#if (historical$Included_in_tess3r[length(historical$Included_in_tess3r)]=="NO") {
#  fillcolor <- "#FFFFFF"
#} else {
#  fillcolor <- "#F2B01F"
#}
#if (historical$status[length(historical$Included_in_tess3r)]=="Hybrid") {
#  colorcolor <- "#9437FF"
#} else {
#  if (historical$status[length(historical$Included_in_tess3r)]=="CC") {
#    colorcolor <- "#15326C"
#  }  else {
#    colorcolor <- "#CE1B26"
#  }
#}

# creating a vector to put in some vlines separating the labels
#vlinepos <- c(0,historical$r[length(historical$r)])
#for (i in (dim(historical)[1]-1):1) {
#  vlinepos <- c(vlinepos, (vlinepos[length(vlinepos)]+(historical$r[i])))
#}

# and adding this annotation on to the baseplot
#historicallat <- historicallat + geom_hline(yintercept = vlinepos, linetype="dashed")

# and adding text label annotation on to the baseplot for left-most site
#historicallat <- historicallat + annotate("label",y=labelx[1],x=labely[1],fontface="bold",label=historical$Location_code[length(historical$Included_in_tess3r)],fill=fillcolor,color=colorcolor,size=12)

# and now doing this for all the remaining labels
#for (i in 2:length(labelx)) {
#  if (historical$Included_in_tess3r[length(historical$Included_in_tess3r)-i+1]=="NO") {
#    fillcolor <- "#FFFFFF"
#  } else {
#    fillcolor <- "#F2B01F"
#  }
#  if (historical$status[length(historical$Included_in_tess3r)-i+1]=="Hybrid") {
#    colorcolor <- "#9437FF"
#  } else {
#    if (historical$status[length(historical$Included_in_tess3r)-i+1]=="CC") {
#      colorcolor <- "#15326C"
#    }  else {
#      colorcolor <- "#CE1B26"
#    }
#  }
#  historicallat <- historicallat + annotate("label",y=labelx[i],x=labely[i],fontface="bold",label=historical$Location_code[length(historical$Included_in_tess3r)-i+1],fill=fillcolor,color=colorcolor,size=12)
#}  

# Exporting as 1500 pixels height, Fig_S1_historical_lat_sites.png
#historicallat

# 13. Generating site labels, finally modern by latitude (the y-axis for the plot)
#modern <- modern %>% arrange(Location_code) %>% arrange(desc(DecimalLatitude))

# replacing "status" with the appropriate color
#fill_codes <- as.matrix(modern$status)[,1]
#fill_codes <- gsub("Hybrid","#9437FF",gsub("CC","#15326C",gsub("BC","#CE1B26",fill_codes)))

# creating a vector for the position of our sampling site text labels
#labelx <- modern$r[length(modern$r)]/2
#for (i in (dim(modern)[1]-1):1) {
#  labelx <- c(labelx, (labelx[length(labelx)]+(modern$r[i+1]/2)+(modern$r[i]/2)))
#}

# creating a vector for the y-position of the labels to stagger them
# and avoid overlap
#labely <- rep(4,length(modern$r))
#labely[seq(2,length(modern$r),4)] <- 1
#labely[seq(3,length(modern$r),4)] <- 3
#labely[seq(4,length(modern$r),4)] <- 0

# building the base colored plot
#modernlat <- ggplot(modern,aes(y=r,x=2)) + geom_bar(stat="identity",color="black",aes(fill=factor(modern$Location_code,levels=as.numeric(modern$Location_code)))) +
#  theme_classic() +
#  theme(panel.border = element_blank(),axis.line=element_blank()) +
#  scale_x_continuous(limits=c(-1,5)) +
#  scale_fill_manual(values = c(fill_codes)) +
#  theme (legend.position="none", axis.text = element_blank(),axis.title=element_blank(), axis.ticks = element_blank())  +
#  theme(aspect.ratio = 5/1)

# Figuring out what the colors should be for fill and outline of the 
# first site label
#if (modern$Included_in_tess3r[length(modern$Included_in_tess3r)]=="NO") {
#  fillcolor <- "#FFFFFF"
#} else {
#  fillcolor <- "#F2B01F"
}
#if (modern$status[length(modern$Included_in_tess3r)]=="Hybrid") {
#  colorcolor <- "#9437FF"
#} else {
#  if (modern$status[length(modern$Included_in_tess3r)]=="CC") {
#    colorcolor <- "#15326C"
#  }  else {
#    colorcolor <- "#CE1B26"
#  }
#}

# creating a vector to put in some vlines separating the labels
#vlinepos <- c(0,modern$r[length(modern$r)])
#for (i in (dim(modern)[1]-1):1) {
#  vlinepos <- c(vlinepos, (vlinepos[length(vlinepos)]+(modern$r[i])))
#}

# and adding this annotation on to the baseplot
#modernlat <- modernlat + geom_hline(yintercept = vlinepos, linetype="dashed")

# and adding text label annotation on to the baseplot for left-most site
#modernlat <- modernlat + annotate("label",y=labelx[1],x=labely[1],fontface="bold",label=modern$Location_code[length(modern$Included_in_tess3r)],fill=fillcolor,color=colorcolor,size=11)

# and now doing this for all the remaining labels
#for (i in 2:length(labelx)) {
#  if (modern$Included_in_tess3r[length(modern$Included_in_tess3r)-i+1]=="NO") {
#    fillcolor <- "#FFFFFF"
#  } else {
#    fillcolor <- "#F2B01F"
#  }
#  if (modern$status[length(modern$Included_in_tess3r)-i+1]=="Hybrid") {
#    colorcolor <- "#9437FF"
#  } else {
#    if (modern$status[length(modern$Included_in_tess3r)-i+1]=="CC") {
#      colorcolor <- "#15326C"
#    }  else {
#      colorcolor <- "#CE1B26"
#    }
#  }
#  modernlat <- modernlat + annotate("label",y=labelx[i],x=labely[i],fontface="bold",label=modern$Location_code[length(modern$Included_in_tess3r)-i+1],fill=fillcolor,color=colorcolor,size=11)
#}  

# Exporting as 1500 pixels height, Fig_S1_modern_lat_sites.png
#modernlat

# Now making the PC1 vs BC_cluster assignment plots
temp <- temp %>% filter(!is.na(BC_genetic_cluster_assignment)) %>% filter(Sampling_period %in% c("MODERN","SMITHSONIAN"))

# Modern first
ggplot() +
  geom_point(temp,
             mapping=aes(x=BC_genetic_cluster_assignment,y=PC1,alpha=Sampling_period,fill=BC_genetic_cluster_assignment),shape=21,color="black",size=15) +
  scale_x_reverse(expand=c(0,0),name = "STRUCTURE assignment to black-capped cluster") + 
  scale_y_reverse(expand=c(0,0), name = "PC1 (81.3% variance explained)") +
  scale_alpha_discrete(range=c(1,0.3)) +
  scale_fill_gradient2(
    low = "#15326C", mid = "#9437FF", high="#CE1B26", midpoint = 0.5
  ) +
  theme_bw(base_size = 28) +
  theme(legend.position = "none") +
  theme(axis.title=element_text(size=36,face="bold")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggsave(filename="Fig1_modern_PC1_vs_structure.pdf",plot = last_plot(),width=38,height=35,units="cm")

# Historical
ggplot() +
  geom_point(temp,
             mapping=aes(x=BC_genetic_cluster_assignment,y=PC1,alpha=Sampling_period,fill=BC_genetic_cluster_assignment),shape=21,color="black",size=15) +
  scale_x_reverse(expand=c(0,0),name = "STRUCTURE assignment to black-capped cluster") + 
  scale_y_reverse(expand=c(0,0), name = "PC1 (81.3% variance explained)") +
  scale_alpha_discrete(range=c(0.3,1)) +
  scale_fill_gradient2(
    low = "#15326C", mid = "#9437FF", high="#CE1B26", midpoint = 0.5
  ) +
  theme_bw(base_size = 28) +
  theme(legend.position = "none") +
  theme(axis.title=element_text(size=36,face="bold")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggsave(filename="Fig1_historical_PC1_vs_structure.pdf",plot = last_plot(),width=38,height=35,units="cm")

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
#  [1] grid      stats     graphics  grDevices utils     datasets 
#[7] methods   base     

#other attached packages:
#  [1] scales_1.2.0    ggsn_0.5.0      ggrepel_0.9.1   ggmap_3.0.0    
#[5] forcats_0.5.1   stringr_1.4.0   dplyr_1.0.9     purrr_0.3.4    
#[9] readr_2.1.2     tidyr_1.2.0     tibble_3.1.7    ggplot2_3.3.6  
#[13] tidyverse_1.3.1

#loaded via a namespace (and not attached):
#  [1] Rcpp_1.0.8.3        lubridate_1.8.0     lattice_0.20-45    
#[4] png_0.1-7           class_7.3-20        digest_0.6.29      
#[7] assertthat_0.2.1    utf8_1.2.2          R6_2.5.1           
#[10] cellranger_1.1.0    plyr_1.8.7          backports_1.4.1    
#[13] reprex_2.0.1        e1071_1.7-11        httr_1.4.3         
#[16] pillar_1.7.0        RgoogleMaps_1.4.5.3 rlang_1.0.2        
#[19] curl_4.3.2          readxl_1.4.0        rstudioapi_0.13    
#[22] labeling_0.4.2      foreign_0.8-82      bit_4.0.4          
#[25] munsell_0.5.0       proxy_0.4-27        broom_0.8.0        
#[28] compiler_4.2.1      modelr_0.1.8        pkgconfig_2.0.3    
#[31] tidyselect_1.1.2    fansi_1.0.3         crayon_1.5.1       
#[34] tzdb_0.3.0          dbplyr_2.2.0        withr_2.5.0        
#[37] sf_1.0-7            bitops_1.0-7        jsonlite_1.8.0     
#[40] gtable_0.3.0        lifecycle_1.0.1     DBI_1.1.3          
#[43] magrittr_2.0.3      units_0.8-0         KernSmooth_2.23-20 
#[46] vroom_1.5.7         cli_3.3.0           stringi_1.7.6      
#[49] farver_2.1.0        fs_1.5.2            sp_1.5-0           
#[52] xml2_1.3.3          ellipsis_0.3.2      generics_0.1.2     
#[55] vctrs_0.4.1         rjson_0.2.21        tools_4.2.1        
#[58] bit64_4.0.5         glue_1.6.2          hms_1.1.1          
#[61] jpeg_0.1-9          parallel_4.2.1      colorspace_2.0-3   
#[64] maptools_1.1-4      classInt_0.4-7      rvest_1.0.2        
#[67] haven_2.5.0 
