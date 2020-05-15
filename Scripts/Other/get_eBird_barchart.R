# From Allen Hurlbert
#
# Need to know eBird’s species and county codes, e.g. as in a Histogram data page URL like this: https://ebird.org/barchart?byr=2019&eyr=2019&bmo=1&emo=12&r=CA-ON-NP


getEbirdBarchartData = function(countyCode, speciesCode = 'reevir1', year) {
  require(data.table)
  url = paste0('https://ebird.org/barchartData?r=', countyCode,
               '&bmo1&emo=12&byr=', year, '&eyr=', year, '&spp=', speciesCode, '&fmt=tsv')
  fileIn = fread(url, skip = 16, header = F)
  
  fileOut = data.frame(date = paste0(rep(year, 48), '-', rep(1:12, each = 4), '-', rep(c(1,8,15,22), times = 12))) %>%
    mutate(speciesCode = speciesCode,
           county = countyCode,
           Year = year,
           julianday = yday(date),
           freq = unlist(fileIn[1, 2:49])) %>%
    select(-date)
  return(fileOut)
}