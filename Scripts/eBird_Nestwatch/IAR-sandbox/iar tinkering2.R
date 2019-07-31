species.list <- as.character(read.csv('/Users/jacobsocolar/Dropbox/Work/Code/Phenomismatch/Bird_Phenology/Data/IAR_species_list.txt', header = F)[,1])

setwd('/Users/jacobsocolar/Desktop/useful_datasets/Processed/IAR_output_2019-05-26')
f <- list.files()
# This script assembles a list of P (currently 500) dataframes, each corresponding to the 24*nth posterior iteration 
# of a posterior chain. Each dataframe has a row for each species-cell-year, giving the predicted arrival and whether 
# the species-cell-year had input data are was estimated exclusively from the IAR model.  


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
                                data.in = as.numeric(rep(NA, counter)), stringsAsFactors = F)
}

# This is horrifically slow, but for now memory pressure on my computer prevents parallelizing.
# I imagine there's a faster way to do this...
counter <- 0
for(s in 1:length(species.list)){
  species <- species.list[s]
  input <- readRDS(paste0(species, '-2019-05-26-iar-stan_input.rds'))
  output <- as.data.frame(readRDS(paste0(species, '-2019-05-26-iar-stan_output.rds')))[(24000/P)*c(1:P), ]
  cell_output <- output[,grep("y_true", names(output))]
  for(y in 1:input$J){
    year <- 2018 - input$J + y
    cols_y <- grep(paste0(',', y, '\\]'), names(cell_output))
    for(c in 1:input$N){
      print(paste0(s, " ", y, " ", c))
      cell <- input$cells[c]
      data.in <- c %in% input$ii_obs[, y]
      counter <- counter + 1
      for(p in 1:P){
        value <- cell_output[p, cols_y[c]]
        
        frame.list[[p]]$species[[counter]] <- species
        frame.list[[p]]$cell[[counter]] <- cell
        frame.list[[p]]$year[[counter]] <- year
        frame.list[[p]]$arrival[[counter]] <- value
        frame.list[[p]]$data.in[[counter]] <- data.in
      }
    }
  }
}
