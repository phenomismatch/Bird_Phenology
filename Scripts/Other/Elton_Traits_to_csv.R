#read in Elton Traits and convert to .csv

setwd('~/Google_Drive/R//Bird_Phenology/Data/Traits/Elton_traits/')

tt <- read.table('BirdFuncDat.txt', 
                 header = TRUE, 
                 sep = '\t',
                 quote = '',
                 fill = TRUE,
                 stringsAsFactors = FALSE)
head(tt)


setwd('~/Google_Drive/R/Bird_Phenology/Data/Traits/Databases/')
write.csv(tt, 'Elton_traits_birds.csv', row.names = FALSE)


#IAR species names to csv
setwd('~/Google_Drive/R/Bird_Phenology/Data/')

tt2 <- read.table('IAR_species_list.txt', 
                 header = FALSE, 
                 stringsAsFactors = FALSE)
head(tt2)

species_list_i2 <- as.vector(apply(tt2, 2, function(x) gsub("_", " ", x)))

setwd('~/Google_Drive/R/Bird_Phenology/Data/Traits/Databases/')
write.csv(species_list_i2, 'IAR_species_list.csv', row.names = FALSE)
