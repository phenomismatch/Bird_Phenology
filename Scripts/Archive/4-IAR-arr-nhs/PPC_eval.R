############################
# PPC evaluation for model fits
#
############################


# create y_PPC function ------------------------------------------------------------


#outputs y_PPC
y_PPC_fun <- function(species, IAR_in_dir)
{
  #Load packages
  library(geosphere)
  library(dplyr)
  library(dggridR)
  
  #Set wd
  setwd(IAR_in_dir)
  s1 <- strsplit(IAR_in_dir, split = '/')[[1]][7]
  IAR_in_date <- strsplit(s1, split = '_')[[1]][3]
  
  #species arg
  args <- as.character(species)
  
  #Filter data by species/years
  #read in master df
  df_master <- readRDS(paste0('IAR_input-', IAR_in_date, '.rds'))
  
  #filter by species and year to be modeled
  f_out <- dplyr::filter(df_master, species == args & MODEL == TRUE)
  
  #define cells and years to be modeled
  cells <- unique(f_out$cell)
  years <- unique(f_out$year)
  nyr <- length(years)
  ncell <- length(cells)
  
  #create adjacency matrix
  #make hexgrid
  hexgrid6 <- dggridR::dgconstruct(res = 6)
  
  #get hexgrid cell centers
  cellcenters <- dggridR::dgSEQNUM_to_GEO(hexgrid6, cells)
  
  #add lat col to df
  f_out$lat <- cellcenters$lat_deg
  
  #create adjacency matrix - 1 if adjacent to cell, 0 if not
  adjacency_matrix <- matrix(data = NA, nrow = ncell, ncol = ncell)
  
  for (i in 1:ncell)
  {
    #i <- 1
    for (j in i:ncell)
    {
      #j <- 69
      dists <- geosphere::distm(c(cellcenters$lon_deg[i], cellcenters$lat_deg[i]),
                                c(cellcenters$lon_deg[j], cellcenters$lat_deg[j]))
      adjacency_matrix[i,j] <- as.numeric((dists/1000) > 0 & (dists/1000) < 311)
    }
  }
  
  #indices for 1s
  ninds <- which(adjacency_matrix == 1, arr.ind = TRUE)
  
  #if a cell doesn't border any other cells, drop it and redefine objects
  DROP <- FALSE
  
  s_cols <- apply(adjacency_matrix, 2, function(x) sum(x, na.rm = TRUE))
  s_rows <- apply(adjacency_matrix, 1, function(x) sum(x, na.rm = TRUE))
  to.rm.ind <- which((s_cols + s_rows) == 0)
  
  if (length(to.rm.ind) > 0)
  {
    DROP <- cells[to.rm.ind]
    
    cells <- cells[-to.rm.ind]
    ncell <- length(cells)
    f_out <- dplyr::filter(f_out, cell %in% cells)
    cellcenters <- dggridR::dgSEQNUM_to_GEO(hexgrid6, cells)
    
    #create adjacency matrix - 1 if adjacent to cell, 0 if not
    adjacency_matrix <- matrix(data = NA, nrow = ncell, ncol = ncell)
    
    for (i in 1:ncell)
    {
      #i <- 1
      for (j in i:ncell)
      {
        #j <- 4
        dists <- geosphere::distm(c(cellcenters$lon_deg[i], cellcenters$lat_deg[i]),
                                  c(cellcenters$lon_deg[j], cellcenters$lat_deg[j]))
        adjacency_matrix[i,j] <- as.numeric((dists/1000) > 0 & (dists/1000) < 311)
      }
    }
    #indices for 1s
    ninds <- which(adjacency_matrix == 1, arr.ind = TRUE)
  }
  
  #fill y_PPC
  y_PPC <- rep(NA, ncell * nyr)
  #counter to fill y_PPC
  counter <- 1
  for (j in 1:nyr)
  {
    #j <- 16
    temp_yr <- dplyr::filter(f_out, year == years[j])
    for (n in 1:ncell)
    {
      #n <- 1
      #matrix with observed values with NAs
      y_PPC[counter] <- temp_yr$HM_mean[n]
      counter <- counter + 1
    }
  }
  return(y_PPC)
}

#outputs y_rep
y_rep_fun <- function(obj, na.vec)
{
  temp_ch <- MCMCvis::MCMCpstr(obj, params = 'y_rep', type = 'chains')[[1]]
  t_temp_ch <- t(temp_ch)
  out <- t_temp_ch[,-na.vec]
  return(out)
}

#posterior predictive check function
PPC_p_fun <- function(y_PPC, y_rep, func = mean)
{
  tvec <- c()
  for (i in 1:NROW(y_rep))
  {
    #i <- 1
    #temp <- sum(t_y_rep[i,] > y_PPC, na.rm = TRUE)
    temp <- func(y_rep[i,], na.rm = TRUE) >
                  func(y_PPC, na.rm = TRUE)
    tvec <- c(tvec, temp)
  }
  
  #PPC p-value
  PPC_p <- mean(tvec)
  return(PPC_p)
}

#mean resid (pred - actual) for each datapoint
resid_ind_fun <- function(y_PPC, y_rep)
{
  temp <- rep(NA, length(y_PPC))
  for (i in 1:length(y_PPC))
  {
    #i <- 1
    temp[i] <- mean(y_rep[,i] - y_PPC[i])
  }
  return(temp)
}



# master function -------------------------------------------------------------------

#species (with underscore), dir with IAR input, dir with IAR output, and dir for fig out 
PPC_fig_fun <- function(SPECIES, IAR_IN, IAR_OUT, OUT_DATE, FIG_OUT)
{
  # SPECIES <- 'Vireo_olivaceus'
  # IAR_IN <- '~/Google_Drive/R/Bird_Phenology/Data/Processed/IAR_input_2019-05-03'
  # IAR_OUT <- '~/Desktop/Bird_Phenology_Offline/Data/Processed/IAR_output_2019-05-26'
  # FIG_OUT <- '~/Desktop/Bird_Phenology_Offline/Figures/PPC_figs'
  
  #get y_PPC
  y_PPC <- y_PPC_fun(species = SPECIES, IAR_in_dir = IAR_IN) 

  #find NAs in y_PPC (missing data)
  na.y.rm <- which(is.na(y_PPC))
  n_y_PPC <- y_PPC[-na.y.rm]

  
  #read in RDS obj
  RDS_OBJ <- readRDS(paste0(IAR_OUT, '/', SPECIES, '-', OUT_DATE, '-iar-stan_output.rds'))
  
  #get y_rep
  y_rep <- y_rep_fun(obj = RDS_OBJ, na.vec = na.y.rm)
  
  #get avg residuals for each pnt
  ind_resid <- resid_ind_fun(n_y_PPC, y_rep)
  
  
  #PPC stats
  #mean
  PPC_mn <- round(PPC_p_fun(n_y_PPC, y_rep, mean), 2)
  mn_bias <- round(mean(ind_resid), 2)
  #sd
  # PPC_sd <- round(PPC_p_fun(n_y_PPC, y_rep, sd), 2)
  # sd_bias <- round(sd(ind_resid), 2)
  
  #additional diagnostics
  num_diverge <- rstan::get_num_divergent(RDS_OBJ)
  model_summary <- MCMCvis::MCMCsummary(RDS_OBJ, Rhat = TRUE, n.eff = TRUE, round = 2)
  rhat_output <- as.vector(model_summary[, grep('Rhat', colnames(model_summary))])
  neff_output <- as.vector(model_summary[, grep('n.eff', colnames(model_summary))])
  
  #PPC plots
  #density overlay plot
  #modified bayesplot::ppc_dens_overlay function
  tdata <- bayesplot::ppc_data(n_y_PPC, y_rep[1:100,])
  
  annotations <- data.frame(xpos = c(-Inf, -Inf, -Inf, -Inf, -Inf),
                            ypos = c(Inf, Inf, Inf, Inf, Inf),
                            annotateText = c(paste0('PPC mean: ', PPC_mn),
                                             paste0('Mean bias: ', mn_bias),
                                             paste0('# divergences: ', num_diverge),
                                             paste0('Max Rhat: ', max(rhat_output)),
                                             paste0('Min n.eff: ', min(neff_output))),
                            hjustvar = c(0, 0, 0, 0, 0),
                            vjustvar = c(4, 6, 8, 10, 12))
  
  p <- ggplot(tdata) + 
    aes_(x = ~value) + 
    stat_density(aes_(group = ~rep_id, color = "yrep"), 
                 data = function(x) dplyr::filter(x, !tdata$is_y), 
                 geom = "line", position = "identity", size = 0.25, 
                 alpha = 0.3, trim = FALSE, bw = 'nrd0', adjust = 1, 
                 kernel = 'gaussian', n = 1024) + 
    stat_density(aes_(color = "y"), data = function(x) dplyr::filter(x, tdata$is_y), 
                 geom = "line", position = "identity", lineend = "round", size = 1, trim = FALSE, 
                 bw = 'nrd0', adjust = 1, kernel = 'gaussian', n = 1024) + 
    theme_classic() + 
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          plot.title = element_text(size = 24)) +
    labs(colour = '') +
    scale_color_manual(values = c('red', 'black')) +
    ggtitle(paste0(SPECIES)) +
    geom_text(data = annotations, aes(x = xpos, y = ypos,
                                      hjust = hjustvar, vjust = vjustvar,
                                      label = annotateText),
              size = 3, col = 'black')
  
  setwd(FIG_OUT)
  ggsave(paste0(SPECIES, '_PPC.pdf'), p)
  
  #bayesplot::ppc_intervals(n_y_PPC, y_rep[1:200,])
  #bayesplot::ppc_ribbon(n_y_PPC, y_rep[1:200,])
  
  #average yrep for each pnt
  yrm <- apply(y_rep, 2, mean)
  pdf(paste0(SPECIES, '_pred_true.pdf'))
  plot(n_y_PPC, yrm, pch = 19, col = rgb(0,0,0,0.4), 
       xlim = range(n_y_PPC, yrm), ylim = range(n_y_PPC, yrm),
       xlab = 'y', ylab = 'y_rep', main = paste0(SPECIES))
  abline(a = 0, b = 1, lty = 2, lwd = 2, col = 'red')
  dev.off()
}




# run function ------------------------------------------------------------

#test
PPC_fig_fun(SPECIES = 'Vireo_bellii',
            IAR_IN = '~/Google_Drive/R/Bird_Phenology/Data/Processed/IAR_input_2019-05-03',
            IAR_OUT = '~/Desktop/',
            OUT_DATE = '2019-06-07',
            FIG_OUT = '~/Desktop/')


setwd('~/Google_Drive/R/Bird_Phenology/Data/')

sp_list <- read.table('IAR_species_list.txt', stringsAsFactors = FALSE)[,1]

#loop through all species
for (i in 1:length(sp_list))
{
  #i <- 2
  if (length(grep(paste0('^', sp_list[i], '.*output.rds$'), 
                  list.files('~/Desktop/Bird_Phenology_Offline/Data/Processed/IAR_output_2019-05-26'))) > 0)
  {
    PPC_fig_fun(SPECIES = sp_list[i],
                IAR_IN = '~/Google_Drive/R/Bird_Phenology/Data/Processed/IAR_input_2019-05-03',
                IAR_OUT = '~/Desktop/Bird_Phenology_Offline/Data/Processed/IAR_output_2019-05-26',
                FIG_OUT = '~/Desktop/Bird_Phenology_Offline/Figures/PPC_figs')
  }
}



