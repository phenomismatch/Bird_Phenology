#########################
#Check skipped years in arrival data
#########################

tt <- readRDS('arrival_master_2019-05-26.rds')
usp <- unique(tt$species)
for (i in 1:length(usp))
{
  #i <- 1
  temp <- dplyr::filter(tt, species == usp[i])
  uy <- unique(temp$year)
  ly <- length(uy)
  dd <- max(uy) - min(uy) + 1
  print(usp[i])
  print(paste0('Missing years: ', (dd - ly)))
  print(uy)
  print('')
}