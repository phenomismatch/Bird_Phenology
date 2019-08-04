species.list <- as.character(read.csv('/Users/jacobsocolar/Dropbox/Work/Code/Phenomismatch/Bird_Phenology/Data/IAR_species_list.txt', header = F)[,1])

setwd('/Users/jacobsocolar/Desktop/useful_datasets/Processed/IAR_output_2019-05-26')
f <- list.files()

# This script assembles a list of P dataframes, each corresponding to an iteration 
# of a posterior chain. Each dataframe has a row for each species-cell-year, giving the predicted arrival, the year
# offset, the combined spatial and nonspatial random effect, and whether the species-cell-year had input data are was estimated exclusively from the IAR model.  


counter <- 0
for(s in 1:length(species.list)){
  species <- species.list[s]
  input <- readRDS(paste0(species, '-2019-05-26-iar-stan_input.rds'))
  counter <- counter + input$J*input$N
}

P <- 500
24000 %% P == 0  # P must evenly divide 24000 or the line that begins "output <-" below will fail.  Can be fixed in that
# if some other number of posterior iterations is really desired.

frame.list <- list()
for(i in 1:P){
  frame.list[[i]] <- data.frame(species = as.character(rep(NA, counter)), cell = as.character(rep(NA, counter)), 
                                year = as.character(rep(NA, counter)), arrival = as.numeric(rep(NA, counter)), 
                                yr_effect = as.numeric(rep(NA, counter)), nu = as.numeric(rep(NA, counter)), 
                                data.in = as.numeric(rep(NA, counter)), stringsAsFactors = F)
}

# There must be a faster way to do this...
counter <- 0
for(s in 1:length(species.list)){
  species <- species.list[s]
  input <- readRDS(paste0(species, '-2019-05-26-iar-stan_input.rds'))
  output <- as.data.frame(readRDS(paste0(species, '-2019-05-26-iar-stan_output.rds')))[(24000/P)*c(1:P), ]
  cell_output <- output[,grep("y_true", names(output))]
  cell_output2 <- output[,grep("nu", names(output))]
  cell_output3 <- output[,grep("beta0\\[", names(output))]
  
  for(y in 1:input$J){
    year <- 2018 - input$J + y
    cols_y <- grep(paste0(',', y, '\\]'), names(cell_output))
    cols_y2 <- grep(paste0(',', y, '\\]'), names(cell_output2))
    for(c in 1:input$N){
      print(paste0(s, " ", y, " ", c))
      cell <- input$cells[c]
      data.in <- c %in% input$ii_obs[, y]
      counter <- counter + 1
      for(p in 1:P){
        arr.day <- cell_output[p, cols_y[c]]
        nu <- cell_output2[p, cols_y2[c]]
        yr.eff <- cell_output3[p, y]
        
        frame.list[[p]]$species[[counter]] <- species
        frame.list[[p]]$cell[[counter]] <- cell
        frame.list[[p]]$year[[counter]] <- year
        frame.list[[p]]$arrival[[counter]] <- arr.day
        frame.list[[p]]$nu[[counter]] <- nu
        frame.list[[p]]$yr_effect[[counter]] <- yr.eff
        frame.list[[p]]$data.in[[counter]] <- data.in
      }
    }
  }
}

save(frame.list, file = "framelist.Rdata")

library(rstanarm)

migD <- read.csv('/Users/JacobSocolar/Dropbox/Work/Phenomismatch/IAR_output/migD.csv')  # gives whether long or short distance migrant, and a summary habitat code for each species
migD[is.na(migD)] <- 0

hex6 <- dggridR::dgconstruct(res = 6)

fl <- frame.list[[1]]
fl$sc <- paste0(fl$species, fl$cell)

fl$lat <- dggridR::dgSEQNUM_to_GEO(hex6,as.numeric(fl$cell))$lat_deg

fl <- fl[fl$data.in == 1, ]
fl$migD <- fl$species %in% migD$species[migD$LDM==1]



model <- lmerTest::lmer(nu ~ yr_effect + lat + yr_effect*lat + (1|sc), data = fl)
summary(model)
model <- lmerTest::lmer(nu ~ yr_effect*lat*migD + (1|sc), data = fl)
summary(model)
model <- lmerTest::lmer(nu ~ yr_effect + lat + yr_effect*lat + (1|sc), data = fl.s)
summary(model)


model <- rstanarm::stan_lmer(nu ~ yr_effect + lat + yr_effect*lat + (1|sc), data = fl)


