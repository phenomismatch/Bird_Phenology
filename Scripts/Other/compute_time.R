##################
# Calculate total compute time for phenology analysis pipeline
#
##################



# 2-halfmax-arr -----------------------------------------------------------

#transfer files to local machine with shell
#scp cyoungflesh@transfer.cam.uchc.edu:/home/CAM/cyoungflesh/phenomismatch/Bird_Phenology/Data/Processed/halfmax_species_2019-03-29/*.o* /Users/caseyyoungflesh/Desktop/Bird_Phenology_Offline/Data/Processed/halfmax_species_2019-03-29


setwd('~/Google_Drive/R/Bird_Phenology/Data/Processed/arrival_GAM_2020-07-10/')

#list files
fls <- list.files()
fls2 <- fls[grep('.out', fls)]

out_2 <- c()
for (i in 1:length(fls2))
{
  #i <- 1
  #get text of total runtime from STDOUT
  rt_temp <- system(paste0('grep Runtime ', fls2[i]), intern = TRUE)
  
  if (length(rt_temp) > 0)
  {
    #strip out text
    min1 <- strsplit(rt_temp, split = ' ')[[1]][4]
    #remove slash and convert to numeric
    min2 <- as.numeric(substring(min1, 1, (nchar(min1) - 1)))
    out_2 <- c(out_2, min2)  
  }
}


#compute time in days
GAM_time <- (sum(out_2)/60)/24



# 4-IAR -------------------------------------------------------------------

#final species list
setwd('~/Google_Drive/R/pheno_trends/Results/arr-gr-SVC-sens-2020-08-07/')
DATA <- readRDS('arr-gr-SVC-sens-data-2020-08-07.rds')
usp <- unique(DATA$mrg_f6$species)

setwd('~/Google_Drive/R/Bird_Phenology/Data/Processed/arrival_IAR_hm_2020-07-21')
fls_4 <- list.files()

fls_5 <- unique(grep(paste(usp, collapse="|"), 
                        fls_4, value = TRUE))

fls_6 <- fls_5[grep('results', fls_5)]

out_6 <- c()
for (i in 1:length(fls_6))
{
  #i <- 1
  print(i)
  rt_temp <- system(paste0('grep minutes ', fls_6[i]), intern = TRUE)
  
  if (length(rt_temp) > 0)
  {
    #strip out text
    min1 <- as.numeric(strsplit(rt_temp, split = ' ')[[1]][3])
    out_6 <- c(out_6, min1)  
  }
}

#compute time in days
IAR_time <- (sum(out_6)/60)/24


#sens model
setwd('~/Google_Drive/R/pheno_trends/Results/arr-gr-SVC-sens-2020-08-07/')
list.files()
rt_temp <- system('grep minutes arr-gr-SVC-sens-stan-results-2020-08-07.txt', intern = TRUE)

#strip out text
min1 <- as.numeric(strsplit(rt_temp, split = ' ')[[1]][3])
sens_time <- min1 / 60 / 24
  
  

#TOTAL DAYS - GAM + IAR + SENS
TDAYS <- GAM_time + IAR_time + sens_time

#TOTAL PROCESSOR DAYS
TDAYS * 4

