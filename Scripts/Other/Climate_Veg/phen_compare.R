library(ggplot2)
library(dggridR)
setwd("/Users/TingleyLab/Dropbox/Work/Phenomismatch/Veg_phenology/Hex_gridded_products/")
setwd("/Users/samuelsocolar/Dropbox/Work/Phenomismatch/Veg_phenology/Hex_gridded_products/")


hexgrid7 <- dgconstruct(res=7)
hexgrid6 <- dgconstruct(res=6)
all_states <- map_data("state")
w2hr <- map_data("world")

usamap <- data.frame(map("world", "USA", plot = FALSE)[c("x", "y")])
canadamap <- data.frame(map("world", "Canada", plot = FALSE)[c("x", "y")])
mexicomap <- data.frame(map("world", "Mexico", plot = FALSE)[c("x", "y")])

load("hex7_avhrr.Rdata")
load("hex7_eMODIS.Rdata")
load("hex7_naphen_si.Rdata")
load("hex7_VIPPHEN.Rdata")
load("hex7_mcd12Q2.Rdata")
load("hex7_mcd12Q2_max.Rdata")
load("hex7_phenocam.Rdata")


hex7list <- list(avhrr=hex7_avhrr, emodis=hex7_eMODIS, si=hex7_naphen_si, vipphen=hex7_VIPPHEN,
                 mcd_increase=hex7_mcd12Q2, mcd_maximum=hex7_mcd12Q2_max, phenocam10=hex7_phenocam$hex7_phenocam10,
                 phenocam25=hex7_phenocam$hex7_phenocam25, phenocam50=hex7_phenocam$hex7_phenocam50)

for(i in 1:8){
  for(j in (i+1):9){
    d1 <- hex7list[[i]]
    d2 <- hex7list[[j]]
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
        if(dim(summary(fit)$coefficients)[1] == 2){
          dr[k] <- sqrt(summary(fit)$r.squared)*(-1+2*as.numeric(summary(fit)$coefficients[2,1] > 0))
          dp[k] <- summary(fit)$coefficients[2,4]
        }else{
          dr[k] <- NA
          dp[k] <- NA
        }
      }else{
        dr[k] <- NA
        dp[k] <- NA
      }

    }
    Dframe <- data.frame(cell=sharedcells, correlation=dr, p_value=dp)
    grid <- dgcellstogrid(hexgrid7, Dframe$cell, frame=TRUE,wrapcells=TRUE)
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
      ggtitle(paste(names(hex7list)[[i]], "vs.", names(hex7list)[[j]]))
    print(p)
  }
}
