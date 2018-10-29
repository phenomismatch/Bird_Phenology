######################
# 3 - proces cubic output
#
# Aggeregate posterior info and diagnostic info from 2-logit-cubic.R to be used in IAR model
#
# Parts formerly in ICAR_parallel.R
######################


# Top-level dir -----------------------------------------------------------

#desktop/laptop
#dir <- '~/Google_Drive/R/'

#Xanadu
dir <- '/home/CAM/cyoungflesh/phenomismatch/'



# db/hm query dir ------------------------------------------------------------

db_dir <- 'db_query_2018-10-15'
hm_dir <- 'halfmax_species_2018-10-16'
IAR_dir <- 'IAR_2018-10-24'


# runtime -----------------------------------------------------------------

tt <- proc.time()



# Load packages -----------------------------------------------------------

library(dplyr)


# Set wd ------------------------------------------------------------------

setwd(paste0(dir, 'Bird_Phenology/Data/'))



# import eBird species list -----------------------------------------------------

species_list_i <- read.table('eBird_species_list.txt', stringsAsFactors = FALSE)
species_list <- species_list_i[,1]
nsp <- length(species_list)


#DATA ONLY VALID THROUGH 2017 (2018 data only goes to ~ jday 60 as of 2018-10-15 query)
years <- 2002:2017
nyr <- length(years)



# combine logit cubic results and diagnostic info -----------------------------------------------------------------

counter <- 0
for (i in 1:nsp)
{
  #i <- 1
  
  #import presence absence ebird data for each specices
  setwd(paste0(dir, 'Bird_Phenology/Data/Processed/', db_dir))
  spdata <- readRDS(paste0('ebird_NA_phen_proc_', species_list[i], '.rds'))
  
  #import halfmax estimates and diagnostics from logit cubic model
  setwd(paste0(dir, 'Bird_Phenology/Data/Processed/', hm_dir))
  temp_halfmax <- readRDS(paste0('halfmax_matrix_list_', species_list[i], '.rds'))
  temp_diag <- readRDS(paste0('halfmax_fit_diag_', species_list[i], '.rds'))
  
  if (i == 1)
  {
    #get number of unique cells
    cells <- unique(spdata$cell6)
    ncel <- length(cells)
    
    #create data.frame to fill
    diagnostics_frame <- as.data.frame(matrix(data = NA, nrow = nsp*ncel*nyr, ncol = 19))
    names(diagnostics_frame) <- c("species", "cell", "year", "n1", "n1W", "n0", "n0i", "njd1", "njd0", "njd0i",
                                  "nphen_bad", "min_n.eff", "max_Rhat", "HM_n.eff", "HM_Rhat",
                                  "HM_mean", "HM_sd", "HM_LCI", "HM_UCI")
  }
  
  #loop through years
  for (j in 1:nyr)
  {
    #j <- 7
    print(paste(i,j))
    ysdata <- dplyr::filter(spdata, year == years[j])
    
    for (k in 1:ncel)
    {
      #k <- 1
      counter <- counter + 1
      diagnostics_frame$species[counter] <- species_list[i]
      diagnostics_frame$year[counter] <- years[j]
      diagnostics_frame$cell[counter] <- cells[k]
      
      cysdata <- dplyr::filter(ysdata, cell6 == cells[k])
      
      #number of surveys where species was detected
      diagnostics_frame$n1[counter] <- sum(cysdata$detect)
      #number of surveys where species was not detected
      diagnostics_frame$n0[counter] <- sum(cysdata$detect == 0)
      #number of detections that came before jday 60
      diagnostics_frame$n1W[counter] <- sum(cysdata$detect*as.numeric(cysdata$day < 60))
      
      if (diagnostics_frame$n1[counter] > 0)
      {
        #number of non-detections before first detection
        diagnostics_frame$n0i[counter] <- length(which(cysdata$detect == 0 & 
                                                         cysdata$day < min(cysdata$day[which(cysdata$detect == 1)])))
        #number of unique days with detections
        diagnostics_frame$njd1[counter] <- length(unique(cysdata$day[which(cysdata$detect == 1)]))
        #number of unique days of non-detections before first detection
        diagnostics_frame$njd0i[counter] <- length(unique(cysdata$day[which(cysdata$detect == 0 & 
                                                                              cysdata$day < min(cysdata$day[which(cysdata$detect == 1)]))]))
      }
      
      #number of unique days with non-detection
      diagnostics_frame$njd0[counter] <- length(unique(cysdata$day[which(cysdata$detect == 0)]))
      
      
      if (diagnostics_frame$n1[counter] > 29 & 
          diagnostics_frame$n1W[counter] < (diagnostics_frame$n1[counter] / 50) &
          diagnostics_frame$n0[counter] > 29 &
          diagnostics_frame$njd0i[counter] > 29 &
          diagnostics_frame$njd1[counter] > 19)
      {
        diagnostics_frame$min_n.eff[counter] <- temp_diag[[j]][[k]]$mineffsSize
        diagnostics_frame$max_Rhat[counter] <- temp_diag[[j]][[k]]$maxRhat
        
        halfmax_posterior <- as.vector(temp_halfmax[[j]][k,])
        
        #convert to mcmc.list and calc n_eff and Rhat using coda (DIFFERENT THAN STAN ESTIMATES)
        halfmax_mcmcList <- coda::mcmc.list(coda::as.mcmc(halfmax_posterior[1:500]), 
                                            coda::as.mcmc(halfmax_posterior[501:1000]),
                                            coda::as.mcmc(halfmax_posterior[1001:1500]), 
                                            coda::as.mcmc(halfmax_posterior[1501:2000]))
        
        diagnostics_frame$HM_n.eff[counter] <- round(coda::effectiveSize(halfmax_mcmcList), digits = 0)
        diagnostics_frame$HM_Rhat[counter] <- round(coda::gelman.diag(halfmax_mcmcList)$psrf[1], digits = 2)
        
        #determine how many estimates are 1 and not 1 (estimates of 1 are bogus)
        diagnostics_frame$nphen_bad[counter] <- sum(halfmax_posterior == 1)
        #halfmax_posterior2 <- halfmax_posterior[which(halfmax_posterior != 1)]
        
        #calculate posterior mean and sd
        diagnostics_frame$HM_mean[counter] <- mean(halfmax_posterior)
        diagnostics_frame$HM_sd[counter] <- sd(halfmax_posterior)
        
        diagnostics_frame$HM_LCI[counter] <- quantile(halfmax_posterior, probs = 0.025)
        diagnostics_frame$HM_UCI[counter] <- quantile(halfmax_posterior, probs = 0.975)
        
        # #fit normal and cauchy distributions to data
        # normfit <- fitdistrplus::fitdist(halfmax_posterior2, "norm")
        # diagnostics_frame$Nloglik[counter] <- normfit$loglik
        # 
        # cfit <- NA
        # #cfit <- tryCatch(fitdistrplus::fitdist(halfmax_posterior2,"cauchy"), error=function(e){return(NA)})
        # if(!is.na(cfit))
        # {
        #   diagnostics_frame$HM_c_loc[counter] <- cfit$estimate[1]
        #   diagnostics_frame$HM_c_scale[counter] <- cfit$estimate[2]
        #   diagnostics_frame$Cloglik[counter] <- cfit$loglik
        # }
      }
    } # k -cell
  } # j - year
} # i - species


# #how many have crazy estimates for halfmax (1 for some iter)
# sum(diagnostics_frame$nphen_bad > 0, na.rm = TRUE)

#add species_cell column
#diagnostics_frame$spCel <- paste(diagnostics_frame$species, diagnostics_frame$cell, sep="_")

#add 'meets criteria' column
diagnostics_frame$m_crit <- NA



# Filter data based on criteria -----------------------------------------------------------

# Which species-years are worth modeling? At a bare minimum:
#     Species with at least 'NC' cells in all three years from 2015-2017
#     Species-years with at least 'NC' cells for those species

NC <- 3

#which species/years meet criteria for model
for (i in 1:length(species_list))
{
  #i <- 1
  #filter by species
  t_sp <- dplyr::filter(diagnostics_frame, species == species_list[i])
  
  #number of cells with good data in each year frmo 2015-2017
  nobs_yr <- c()
  for (j in 2015:2017)
  {
    #j <- 2017
    ty_sp <- dplyr::filter(t_sp, year == j)
    ind <- which(!is.na(ty_sp$HM_mean))
    nobs_yr <- c(nobs_yr, length(ind))
    #ty_sp[ind,]
  }
  
  
  #if all three years have greater than or = to 'NC' cells of data, figure out which years have at least 'NC' cells
  if (sum(nobs_yr >= NC) == 3)
  {
    #see which years have more than 3 cells of data
    nobs_yr2 <- c()
    for (j in min(years):max(years))
    {
      #j <- 2012
      ty_sp2 <- dplyr::filter(t_sp, year == j)
      ind2 <- which(!is.na(ty_sp2$HM_mean))
      nobs_yr2 <- c(nobs_yr2, length(ind2))
    }
    
    #years to keep (more than three cells of data)
    yrs_kp <- years[which(nobs_yr2 >= NC)]
    
    #of this species, which years to keep
    t_sp_kp <- t_sp[which(t_sp$year %in% yrs_kp),]
    
    #figure out which cells can be classified as non-winter cells
    for (k in 1:ncel)
    {
      #k <- 1
      t_sp_kp_c <- dplyr::filter(t_sp_kp, cell == cells[k])
      
      t_sp_kp_c$n1[counter] > 29 
      t_sp_kp_c$n1W[counter] < (diagnostics_frame$n1[counter] / 50) &
        t_sp_kp_c$n0[counter] > 29 &
      
      
      
      
      #of the cells for this particular species/year, does this trigger the winter threshold?
      t_sp_kp_c_f <- t_sp_kp_c[which(!is.na(t_sp_kp_c$HM_mean)),]
      
      #make sure 
      if (sum(t_sp_kp_c_f$n1W < (t_sp_kp_c_f$n1 / 50)) > 0)
    }
    
    
    
    #add m_crit == TRUE if that species/year meets criteria
    diagnostics_frame[which(diagnostics_frame$species == species_list[i] & 
                              diagnostics_frame$year %in% yrs_kp),]$m_crit <- TRUE
  }
}


# write diagnostics data.frame to RDS

IAR_dir_path <- paste0(dir, 'Bird_phenology/Data/Processed/IAR_', Sys.Date())

dir.create(IAR_dir_path)
setwd(IAR_dir_path)

saveRDS(diagnostics_frame, paste0('diagnostics_frame-', Sys.Date(), '.rds'))



# create dfs that show # cells with data in each year/species, and # years with data in each cell/species -----------------

#create vector of year names
yrs_vec <- c()
for (i in min(years):max(years))
{
  yrs_vec <- c(yrs_vec, paste0('yr_', i))
}

#create vector of cell names
cell_vec <- c()
for (i in 1:ncel)
{
  cell_vec <- c(cell_vec, paste0('cell_', cells[i]))
}


#create df with species/cells/n_yrs per sp_cell
cells_frame <- data.frame(species = rep(species_list, each = length(cells)), 
                          cell = rep(cells, length(species_list)), 
                          n_yrs = NA)

yrs_frame <- data.frame(species = rep(species_list, each = length(years)), 
                        year = rep(years, length(species_list)), 
                        n_cells = NA)


#add years and cell cols
cells_frame[yrs_vec] <- NA
yrs_frame[cell_vec] <- NA

#fill cells_frame and yrs_frame
for (i in 1:nsp)
{
  #i <- 80
  
  for (k in 1:ncel)
  {
    #k <- 30
    t_cell <- dplyr::filter(diagnostics_frame, species == species_list[i], cell == cells[k])
    
    #get index for species/cell row
    idx <- which(cells_frame$species == species_list[i] & 
                   cells_frame$cell == cells[k])
    
    #insert number of yrs with data
    yrs_d <- t_cell$year[which(!is.na(t_cell$HM_mean))]
    cells_frame[idx,'n_yrs'] <- length(yrs_d)
    
    #find out which years have data and insert TRUE into df - leave NA otherwise
    temp_yrs <- paste0('yr_', yrs_d)
    cells_frame[idx, which(colnames(cells_frame) %in% temp_yrs)] <- TRUE
  }
  
  for (j in min(years):max(years))
  {
    #j <- 2015
    t_yr <- dplyr::filter(diagnostics_frame, species == species_list[i], year == j)
    
    #get index for species/year row
    idx_yr <- which(yrs_frame$species == species_list[i] & 
                      yrs_frame$year == j)
    
    #insert number of cells with data
    cells_d <- t_yr$cell[which(!is.na(t_yr$HM_mean))]
    yrs_frame[idx_yr,'n_cells'] <- length(cells_d)
    
    #find out which years have data and insert TRUE into df - leave NA otherwise
    temp_cells <- paste0('cell_', cells_d)
    yrs_frame[idx_yr, which(colnames(yrs_frame) %in% temp_cells)] <- TRUE
  }
}

#write to RDS

saveRDS(cells_frame, paste0('cells_frame-', Sys.Date(), '.rds'))
saveRDS(yrs_frame, paste0('yrs_frame-', Sys.Date(), '.rds'))

