
#terminal code for imagemagick:
# convert -density 300 -delay 100 -resize 1000x1000 *.pdf ~/Desktop/IAR.gif

IAR_species <- as.character(read.table('~/Google_Drive/R/Bird_Phenology/Data/IAR_species_list.txt')[,1])

setwd('~/Google_Drive/R/Bird_Phenology/Scripts/eBird_Nestwatch/AOS/')
for (i in 1:length(IAR_species))
{
  #i <- 1
  print(IAR_species[i])
  system(paste0('Rscript gif_plotting.R ', IAR_species[i]))
}