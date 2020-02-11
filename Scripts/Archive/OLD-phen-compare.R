######################
# 5-phen-compare - correlations among phenology products
#
# Correlations between veg phenology products and bird phenology, plotted on hex grid
#
# formerly phen_compare_hex6.R
######################

# This script examines the interannual correlations between different phenological products on a moderately
# coarse hexagonal grid. In addition to the MODIS, AVHRR, and SI veg phenology products, this script
# brings in the all bird index from 2009 onward. In the future, this script will also include
# some nestwatch results.

library(ggplot2)
library(dggridR)
setwd("/Users/TingleyLab/Dropbox/Work/Phenomismatch/Veg_phenology/Hex_gridded_products/")

hexgrid7 <- dgconstruct(res=7)
hexgrid6 <- dgconstruct(res=6)
all_states <- map_data("state")
w2hr <- map_data("world")

usamap <- data.frame(map("world", "USA", plot = FALSE)[c("x", "y")])
canadamap <- data.frame(map("world", "Canada", plot = FALSE)[c("x", "y")])
mexicomap <- data.frame(map("world", "Mexico", plot = FALSE)[c("x", "y")])

load("hex6_avhrr.Rdata")
load("hex6_eMODIS.Rdata")
load("hex6_naphen_si.Rdata")
load("hex6_VIPPHEN.Rdata")
load("hex6_mcd12Q2.Rdata")
load("hex6_mcd12Q2_max.Rdata")
load("hex6_phenocam.Rdata")
load("/Users/TingleyLab/Dropbox/Work/Phenomismatch/NA_birdPhen/hex6_birdAll2009.Rdata")
load("/Users/TingleyLab/Dropbox/Work/Phenomismatch/NA_birdPhen/hex6_nestwatch.Rdata")


names(hex6_mcd12Q2)[1] <- "cell"
names(hex6_mcd12Q2_max)[1] <- "cell"

hex6list <- list(avhrr=hex6_avhrr, emodis=hex6_eMODIS, si=hex6_naphen_si, vipphen=hex6_VIPPHEN,
                 mcd_increase=hex6_mcd12Q2, mcd_maximum=hex6_mcd12Q2_max, phenocam_10=hex6_phenocam$hex6_phenocam10,
                 phenocam_25=hex6_phenocam$hex6_phenocam25, phenocam_50=hex6_phenocam$hex6_phenocam50,
                 birdAll2009=hex6_birdAll2009)

for(i in 1:8){
  for(j in (i+1):9){
    d1 <- hex6list[[i]]
    d2 <- hex6list[[j]]
    sharedcells <- unique(d1$cell[which(d1$cell %in% d2$cell)])
    sharedyears <- unique(d1$year[which(d1$year %in% d2$year)])
    d1s <- d1[which(d1$cell %in% sharedcells & d1$year %in% sharedyears),]
    d2s <- d2[which(d2$cell %in% sharedcells & d2$year %in% sharedyears),]
    D <- dplyr::inner_join(d1s, d2s, by = c('cell', 'year'))
    dr <- dp <- rep(NA, length(sharedcells))
    for(k in 1:length(sharedcells)){
      DCs <- D[which(D$cell == sharedcells[k]),]
      if(dim(DCs)[1] > 2){
        fit <- tryCatch(lm(DCs[,3]~DCs[,4]), error=function(e){return(NA)})
      }else{
        dr[k] <- NA
        dp[k] <- NA
      }
      
      if(!is.na(fit)){
        dr[k] <- sqrt(summary(fit)$r.squared)*(-1+2*as.numeric(summary(fit)$coefficients[2,1] > 0))
        dp[k] <- summary(fit)$coefficients[2,4]
      }else{
        dr[k] <- NA
        dp[k] <- NA
      }
      
    }
    Dframe <- data.frame(cell=sharedcells, correlation=dr, p_value=dp, cp=dr*(dp<.05))
    grid <- dgcellstogrid(hexgrid6, Dframe$cell, frame=TRUE,wrapcells=TRUE)
    grid <- merge(grid, Dframe, by = "cell")
    
    p<- ggplot() + 
      geom_polygon(data=grid,      aes(x=long, y=lat, group=group, fill=correlation), alpha=1)  +
      #      geom_path   (data=grid,      aes(x=long, y=lat, group=group), alpha=1, color="white", size=.2) +
      scale_fill_gradient2(low = "red", high = "blue", limits = c(-1,1)) +
      #      geom_path(data=all_states, aes(x=long, y=lat, group = group), color="black", alpha=.4, size = .2) +
      geom_path(data = usamap, aes(x, y), colour = "black", alpha = .5, size = .2) +
      geom_path(data = canadamap, aes(x, y), colour = "black", alpha = .5, size = .2) +
      geom_path(data = mexicomap, aes(x, y), colour = "black", alpha = .5, size = .2) +
      coord_map("ortho", orientation = c(35, -95, 0))+
      xlab('')+ylab('')+
      theme(axis.ticks.x=element_blank())+
      theme(axis.ticks.y=element_blank())+
      theme(axis.text.x=element_blank())+
      theme(axis.text.y=element_blank())+
      ggtitle(paste(names(hex6list)[[i]], "vs.", names(hex6list)[[j]]))
    print(p)
    
  }
}



d2 <- hex6_birdAll2009

for(i in 1:9){
    d1 <- hex6list[[i]]
    sharedcells <- unique(d1$cell[which(d1$cell %in% d2$cell)])
    sharedyears <- unique(d1$year[which(d1$year %in% d2$year)])
    d1s <- d1[which(d1$cell %in% sharedcells & d1$year %in% sharedyears),]
    d2s <- d2[which(d2$cell %in% sharedcells & d2$year %in% sharedyears),]
    D <- dplyr::inner_join(d1s, d2s, by = c('cell', 'year'))
    dr <- dp <- rep(NA, length(sharedcells))
    for(k in 1:length(sharedcells)){
      DCs <- D[which(D$cell == sharedcells[k]),]
      if(dim(DCs)[1] > 2){
        fit <- tryCatch(lm(DCs[,4]~DCs[,3], weights = 1/(DCs[,5]^2)), error=function(e){return(NA)})
      }else{
        dr[k] <- NA
        dp[k] <- NA
      }
      
      if(!is.na(fit)){
        dr[k] <- sqrt(summary(fit)$r.squared)*(-1+2*as.numeric(summary(fit)$coefficients[2,1] > 0))
        dp[k] <- summary(fit)$coefficients[2,4]
      }else{
        dr[k] <- NA
        dp[k] <- NA
      }
      
    }
    Dframe <- data.frame(cell=sharedcells, correlation=dr, p_value=dp, cp=dr*(dp<.05))
    grid <- dgcellstogrid(hexgrid6, Dframe$cell, frame=TRUE,wrapcells=TRUE)
    grid <- merge(grid, Dframe, by = "cell")
    
    p<- ggplot() + 
      geom_polygon(data=grid,      aes(x=long, y=lat, group=group, fill=correlation), alpha=1)  +
      #      geom_path   (data=grid,      aes(x=long, y=lat, group=group), alpha=1, color="white", size=.2) +
      scale_fill_gradient2(low = "red", high = "blue", limits = c(-1,1)) +
      #      geom_path(data=all_states, aes(x=long, y=lat, group = group), color="black", alpha=.4, size = .2) +
      geom_path(data = usamap, aes(x, y), colour = "black", alpha = .5, size = .2) +
      geom_path(data = canadamap, aes(x, y), colour = "black", alpha = .5, size = .2) +
      geom_path(data = mexicomap, aes(x, y), colour = "black", alpha = .5, size = .2) +
      coord_map("ortho", orientation = c(35, -95, 0))+
      xlab('')+ylab('')+
      theme(axis.ticks.x=element_blank())+
      theme(axis.ticks.y=element_blank())+
      theme(axis.text.x=element_blank())+
      theme(axis.text.y=element_blank())+
      ggtitle(paste("birds 2009-2016 vs.", names(hex6list)[[i]]))
    print(p)
}

d2 <- hex6_nestwatch
d2 <- hex6_nestwatch[which(hex6_nestwatch$n >= 4), ]

for(i in 1:10){
  d1 <- hex6list[[i]]
  sharedcells <- unique(d1$cell[which(d1$cell %in% d2$cell)])
  sharedyears <- unique(d1$year[which(d1$year %in% d2$year)])
  d1s <- d1[which(d1$cell %in% sharedcells & d1$year %in% sharedyears),]
  d2s <- d2[which(d2$cell %in% sharedcells & d2$year %in% sharedyears),]
  D <- dplyr::inner_join(d1s, d2s, by = c('cell', 'year'))
  dr <- dp <- rep(NA, length(sharedcells))
  for(k in 1:length(sharedcells)){
    DCs <- D[which(D$cell == sharedcells[k]),]
    if(dim(DCs)[1] > 3){
      fit <- tryCatch(lm(DCs[,4]~DCs[,3]), error=function(e){return(NA)})
    }else{
      dr[k] <- NA
      dp[k] <- NA
    }
    
    if(!is.na(fit)){
      dr[k] <- sqrt(summary(fit)$r.squared)*(-1+2*as.numeric(summary(fit)$coefficients[2,1] > 0))
      dp[k] <- summary(fit)$coefficients[2,4]
    }else{
      dr[k] <- NA
      dp[k] <- NA
    }
    
  }
  Dframe <- data.frame(cell=sharedcells, correlation=dr, p_value=dp, cp=dr*(dp<.05))
  grid <- dgcellstogrid(hexgrid6, Dframe$cell, frame=TRUE,wrapcells=TRUE)
  grid <- merge(grid, Dframe, by = "cell")
  
  p<- ggplot() + 
    geom_polygon(data=grid,      aes(x=long, y=lat, group=group, fill=correlation), alpha=1)  +
    #      geom_path   (data=grid,      aes(x=long, y=lat, group=group), alpha=1, color="white", size=.2) +
    scale_fill_gradient2(low = "red", high = "blue", limits = c(-1,1)) +
    #      geom_path(data=all_states, aes(x=long, y=lat, group = group), color="black", alpha=.4, size = .2) +
    geom_path(data = usamap, aes(x, y), colour = "black", alpha = .5, size = .2) +
    geom_path(data = canadamap, aes(x, y), colour = "black", alpha = .5, size = .2) +
    geom_path(data = mexicomap, aes(x, y), colour = "black", alpha = .5, size = .2) +
    coord_map("ortho", orientation = c(35, -95, 0))+
    xlab('')+ylab('')+
    theme(axis.ticks.x=element_blank())+
    theme(axis.ticks.y=element_blank())+
    theme(axis.text.x=element_blank())+
    theme(axis.text.y=element_blank())+
    ggtitle(paste("nestwatch TRES vs.", names(hex6list)[[i]]))
  print(p)
  
}
