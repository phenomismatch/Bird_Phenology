nlcd2016 <- raster::raster('/Users/jacobsocolar/Desktop/useful_datasets/NLCD/NLCD_Land_Cover_L48_2019424_full_zip/NLCD_2016_Land_Cover_L48_20190424.img')
nlcd2013 <- raster::raster('/Users/jacobsocolar/Desktop/useful_datasets/NLCD/NLCD_Land_Cover_L48_2019424_full_zip/NLCD_2013_Land_Cover_L48_20190424.img')
nlcd2011 <- raster::raster('/Users/jacobsocolar/Desktop/useful_datasets/NLCD/NLCD_Land_Cover_L48_2019424_full_zip/NLCD_2011_Land_Cover_L48_20190424.img')
nlcd2008 <- raster::raster('/Users/jacobsocolar/Desktop/useful_datasets/NLCD/NLCD_Land_Cover_L48_2019424_full_zip/NLCD_2006_Land_Cover_L48_20190424.img')
nlcd2006 <- raster::raster('/Users/jacobsocolar/Desktop/useful_datasets/NLCD/NLCD_Land_Cover_L48_2019424_full_zip/NLCD_2006_Land_Cover_L48_20190424.img')
nlcd2004 <- raster::raster('/Users/jacobsocolar/Desktop/useful_datasets/NLCD/NLCD_Land_Cover_L48_2019424_full_zip/NLCD_2004_Land_Cover_L48_20190424.img')
nlcd2001 <- raster::raster('/Users/jacobsocolar/Desktop/useful_datasets/NLCD/NLCD_Land_Cover_L48_2019424_full_zip/NLCD_2001_Land_Cover_L48_20190424.img')

begin <- proc.time()
decid <- (nlcd2001 %in% c(41,43)) & (nlcd2004 %in% c(41,43)) & (nlcd2006 %in% c(41,43)) &
  (nlcd2008 %in% c(41,43)) & (nlcd2011 %in% c(41,43)) & (nlcd2013 %in% c(41,43)) &
  (nlcd2016 %in% c(41,43))
proc.time() - begin






phenextract <- function(raster1, raster2, aggregate = 'n',  nlcd_years = 'all'){
  if(!is.na(nlcd_years)){

  }
}

nlcd2016 <- raster::raster('/Users/jacobsocolar/Desktop/useful_datasets/NLCD/NLCD_Land_Cover_L48_2019424_full_zip/NLCD_2016_Land_Cover_L48_20190424.img')
nlcd2001 <- raster::raster('/Users/jacobsocolar/Desktop/useful_datasets/NLCD/NLCD_Land_Cover_L48_2019424_full_zip/NLCD_2001_Land_Cover_L48_20190424.img')
s <- (ncld2016 %in% c(41, 43)) & (ncld2001 %in% c(41, 43))
extent(nlcd2016) == extent(nlcd2001)
