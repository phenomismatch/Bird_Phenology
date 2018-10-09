# This script imports data from the eBird reference dataset (ERD) for a vector of species, defined by 
# species_list and a vector of years defined by years.  The script is designed to interact with the ERD
# in the format downloaded directly from eBird, with a separate folder for each year from 2002 onwards.

years <- 2002:2016
nyr <- length(years)

species_list <- c("Empidonax_virescens", "Myiarchus_crinitus", "Contopus_virens", "Vireo_olivaceus",
                  "Vireo_solitarius", "Vireo_gilvus", "Vireo_flavifrons", "Catharus_fuscescens",
                  "Dumetella_carolinensis", "Setophaga_dominica", "Limnothlypis_swainsonii",
                  "Setophaga_citrina", "Geothlypis_formosa", "Parkesia_motacilla", "Parkesia_noveboracensis",
                  "Mniotilta_varia", "Setophaga_americana", "Setophaga_ruticilla", "Setophaga_virens",
                  "Setophaga_virens", "Setophaga_caerulescens", "Protonotaria_citrea", "Setophaga_cerulea",
                  "Seiurus_aurocapilla", "Cardellina_canadensis", "Piranga_olivacea", "Piranga_rubra",
                  "Pheucticus_ludovicianus", "Icterus_galbula", "Empidonax_traillii", "Empidonax_alnorum",
                  "Empidonax_minimus", "Tyrannus_tyrannus", "Vireo_bellii", "Vireo_griseus", "Tachycineta_bicolor",
                  "Stelgidopteryx_serripennis","Hirundo_rustica","Riparia_riparia","Petrochelidon_pyrrhonota",
                  "Progne_subis", "Vermivora_cyanoptera","Vermivora_chrysoptera","Oreothlypis_ruficapilla",
                  "Setophaga_pensylvanica", "Setophaga_petechia", "Setophaga_discolor", "Geothlypis_philadelphia",
                  "Pooecetes_gramineus", "Ammodramus_nelsoni", "Passerina_cyanea",
                  "Passerina_caerulea","Spiza_americana","Icterus_spurius","Dolichonyx_oryzivorus","Contopus_cooperi",
                  "Empidonax_flaviventris","Regulus_satrapa","Regulus_calendula","Vireo_philadelphicus",
                  "Troglodytes_hiemalis","Catharus_guttatus","Catharus_ustulatus","Catharus_bicknelli",
                  "Catharus_minimus","Setophaga_fusca","Setophaga_striata","Setophaga_tigrina","Oreothlypis_peregrina",
                  "Setophaga_castanea","Setophaga_palmarum","Oreothlypis_celata","Cardellina_pusilla",
                  "Oporornis_agilis","Setophaga_magnolia","Setophaga_coronata","Passerella_iliaca",
                  "Melospiza_lincolnii","Spizelloides_arborea","Junco_hyemalis","Zonotrichia_leucophrys",
                  "Zonotrichia_albicollis","Euphagus_carolinus","Ictinia_mississippiensis",
                  "Elanoides_forficatus","Archilochus_colubris","Coccyzus_erythropthalmus","Coccyzus_americanus",
                  "Antrostomus_vociferus","Antrostomus_carolinensis","Sphyrapicus_varius","Sayornis_phoebe",
                  "Bombycilla_cedrorum","Cyanocitta_cristata","Turdus_migratorius","Polioptila_caerulea",
                  "Setophaga_pinus","Peucaea_aestivalis","Spinus_tristis",
                  "Geothlypis_trichas","Icteria_virens","Mimus_polyglottos","Lanius_ludovicianus",
                  "Melospiza_georgiana","Quiscalus_quiscula","Passerculus_sandwichensis","Ammodramus_caudacutus",
                  "Spizella_passerina","Pipilo_erythrophthalmus","Troglodytes_aedon","Sturnella_magna",
                  "Sialia_sialis","Cistothorus_palustris","Corvus_ossifragus","Agelaius_phoeniceus")
nsp <- length(species_list)

# # run awk script directly from R to get the column indicies of the species of interest, make sure
# # the indicies are the same across all year files:
# species_indices <- as.data.frame(matrix(data=NA,nrow=length(species_list), ncol=length(years)))
# for(j in 1:length(years)){
#   setwd(paste0("/Users/TingleyLab/Desktop/useful_datasets/eBird/ERD2016SS/",years[j],"/"))
#   for(i in 1:length(species_list)){
#     unix_command <- paste0("head -1 checklists.csv | awk -v RS=\",\" '/",species_list[i],
#                            "/{print NR;}'")
#     species_indices[i,j] <- system(unix_command, intern = T)
#   }
# }
# equalTest <- vector()
# for(j in 2:length(years)){
#   equalTest[j-1] <- all.equal(species_indices[,1], species_indices[,j])
# }
# sum(equalTest)
# spInd <- as.vector(species_indices[,1])
# 
# # For each year, run awk script to extract just the columns of interest
# unix_command <- paste0("cut -d \",\" -f1-8,12-19,",paste(spInd,collapse=","), " checklists.csv > checklists_NA_birdPhen.csv")
# for(j in 1:length(years)){
#   setwd(paste0("/Users/TingleyLab/Desktop/useful_datasets/eBird/ERD2016SS/",years[j],"/"))
#   system(unix_command)
# }
# 
# # For each year, run awk script to extract only observations that occurred on or before jday 200:
# setwd("/Users/TingleyLab/Desktop/useful_datasets/eBird/ERD2016SS/")
# for(j in 1:length(years)){
#   unix_command <- paste0("awk -f extract_awk_time.txt ",years[j],"/checklists_NA_birdPhen.csv > ",
#                          years[j],"/checklists_NA_birdPhen_time.csv")
#   system(unix_command)
# }



# For each year, import the data, and store in a list:
dataList <- as.list(rep(NA, nyr))
for(i in 1:nyr){
  assign(paste0("data", years[i]), read.csv(paste0("/Users/TingleyLab/Desktop/useful_datasets/eBird/ERD2016SS/",years[i], "/checklists_NA_birdPhen_time.csv")))
  dataList[[i]] <- eval(parse(text = paste("data", years[i], sep="")))
  rm(list = paste("data", years[i], sep=""))
  print(years[i])
}

# rbind dataframes to a single dataframe and save:
data_NA_birdPhen <- data.table::rbindlist(dataList)
save(data_NA_birdPhen, file = '/Users/Tingleylab/Dropbox/Work/Phenomismatch/data_NA_birdPhen.Rdata')
