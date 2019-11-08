##################
# Calculate total compute time for phenology analysis pipeline
#
##################



# 2-halfmax-arr -----------------------------------------------------------

#transfer files to local machine with shell
#scp cyoungflesh@transfer.cam.uchc.edu:/home/CAM/cyoungflesh/phenomismatch/Bird_Phenology/Data/Processed/halfmax_species_2019-03-29/*.o* /Users/caseyyoungflesh/Desktop/Bird_Phenology_Offline/Data/Processed/halfmax_species_2019-03-29


setwd('~/Desktop/Bird_Phenology_Offline/Data/Processed/halfmax_species_2019-03-29/')

#list files
fls <- list.files()

#filter out .out files that were cancelled or had errors (those were rerun but had different STDOUT names because of xanadu wierdness)
fls2 <- fls[-grep('CANCEL', fls)]
fls3 <- fls2[-grep('ERROR', fls2)]

out_2 <- c()
for (i in 1:length(fls3))
{
  #i <- 22
  #get text of total runtime from STDOUT
  rt_temp <- system(paste0('grep Runtime ', fls3[i]), intern = TRUE)
  
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
(sum(out_2)/60)/24



# 4-IAR -------------------------------------------------------------------

setwd('~/Google_Drive/R/Bird_Phenology/Data/Processed/IAR_output_2019-05-26/')

fls_4 <- list.files()
fls_42 <- fls_4[grep('results', fls_4)]

out_4 <- c()
for (i in 1:length(fls_42))
{
  #i <- 2
  rt_temp <- system(paste0('grep minutes ', fls_42[i]), intern = TRUE)
  
  if (length(rt_temp) > 0)
  {
    #strip out text
    min1 <- as.numeric(strsplit(rt_temp, split = ' ')[[1]][3])
    out_4 <- c(out_4, min1)  
  }
}

#compute time in days
(sum(out_4)/60)/24


#TOTAL MINUTES
t_minutes <- sum(c(out_2, out_4))

#TOTAL HOURS
t_hours <- t_minutes/60

#TOTAL DAYS
t_days <- t_hours/24

#TOTAL YEARS
t_days/365
