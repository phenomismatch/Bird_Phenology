library(dggridR)
hexgrid6 <- dgconstruct(res=6) # Construct geospatial hexagonal grid
dg_closest_res_to_spacing(hexgrid6, 200) # In terms of spacing, a resolution of 7 would be closer to Mayor et al's 200km, but in terms of cell area Mayor's cells fall between resolution 6 and 7. To keep computation time low, we'll use the coarser resolution (6) for now.

load('/Users/Tingleylab/Dropbox/Work/Phenomismatch/dataAll.Rdata')

dataAll$cell <- dgGEO_to_SEQNUM(hexgrid6, dataAll$LONGITUDE, dataAll$LATITUDE)$seqnum
dataAll2 <- dataAll[which(dataAll$MONTH < 8), ]
dataAll2 <- dataAll2[which(dataAll2$LONGITUDE > -100 & dataAll2$LATITUDE > 32), ]
dataAll2 <- as.data.frame(dataAll2)
cells <- unique(dataAll2$cell)
ncell <- length(cells)
years <- c(2002:2014)
nyr <- length(years)
species_list <- c("Empidonax_virescens", "Myiarchus_crinitus", "Contopus_virens", "Vireo_olivaceus",
                  "Vireo_solitarius", "Vireo_flavifrons", "Vireo_gilvus", "Catharus_fuscescens",
                  "Dumetella_carolinensis", "Setophaga_dominica", "Limnothlypis_swainsonii", "Setophaga_citrina",
                  "Geothlypis_formosa", "Parkesia_motacilla", "Parkesia_noveboracensis", "Mniotilta_varia",
                  "Setophaga_americana", "Setophaga_ruticilla", "Helmitheros_vermivorum", "Setophaga_virens",
                  "Setophaga_caerulescens", "Protonotaria_citrea", "Setophaga_cerulea", "Seiurus_aurocapilla",
                  "Cardellina_canadensis", "Piranga_rubra", "Piranga_olivacea", "Pheucticus_ludovicianus",
                  "Icterus_galbula")
nsp <- length(species_list)

model_results <- as.data.frame(matrix(data=NA, nrow=nyr*nsp*ncell, ncol=17))
colnames(model_results) <- c("year", "cell", "species", "n0", "n1", "n1J", "GAM_peak", "GAM_half", "pol2_peak", "pol2_half", "pol3_peak", "pol3_half",
                          "pol4_peak", "pol4_half", "pol5_peak", "pol5_half", "sLogi_half")
counter <- 0
for(i in 1:nyr){
  for(j in 1:ncell){
    data_ij <- dataAll2[which(dataAll2$YEAR==years[i] & dataAll2$cell==cells[j]),]
    for(k in 1:nsp){
      counter <- counter+1
      model_results$year[counter] <- years[i]
      model_results$cell[counter] <- cells[j]
      model_results$species[counter] <- species_list[k]
      data_ijk <- data_ij[ , as.numeric(which(colnames(data_ij) %in% c('DAY', species_list[k])))]
      if('X' %in% data_ijk[,2]){
        data_ijk[which(data_ijk[,2] == 'X'),2] <- 1
      }
      data_ijk[,2] <- as.numeric(as.numeric(as.character(data_ijk[,2])) > 0)
      n1 <- sum(data_ijk[,2])
      n0 <- dim(data_ijk)[1] - sum(1-data_ijk[,2])
      n1J <- sum(data_ijk[,2]*as.numeric(data_ijk$DAY > 160))
      model_results$n0[counter] <- n0
      model_results$n1[counter] <- n1
      model_results$n1J[counter] <- n1J
      if(n1 > 500 & n0 > 500 & n1J > 100){
        colnames(data_ijk)[2] <- "PrAb"
        gam_ijk <- gam(PrAb ~ s(DAY), family='binomial', data = data_ijk)
        predictions <- predict(gam_ijk, newdata=newdata)
        transf_pred <- 1/(1+exp(-predictions))
        model_results$GAM_peak[counter] <- min(which(transf_pred==max(transf_pred)))
        model_results$GAM_half[counter] <- min(which(transf_pred > max(transf_pred)/2))
        data_ijk$Day2 <- data_ijk$DAY^2
        data_ijk$Day3 <- data_ijk$DAY^3
        data_ijk$Day4 <- data_ijk$DAY^4
        data_ijk$Day5 <- data_ijk$DAY^5
        
        poly_ijk <- glm(PrAb ~ DAY + Day2, family='binomial', data = data_ijk)
        predictions <- predict(poly_ijk, newdata=newdata)
        transf_pred <- 1/(1+exp(-predictions))
        model_results$pol2_peak[counter] <- min(which(transf_pred==max(transf_pred)))
        model_results$pol2_half[counter] <- min(which(transf_pred > max(transf_pred)/2))
        
        poly_ijk <- glm(PrAb ~ DAY + Day2 + Day3, family='binomial', data = data_ijk)
        predictions <- predict(poly_ijk, newdata=newdata)
        transf_pred <- 1/(1+exp(-predictions))
        model_results$pol3_peak[counter] <- min(which(transf_pred==max(transf_pred)))
        model_results$pol3_half[counter] <- min(which(transf_pred > max(transf_pred)/2))
        
        poly_ijk <- glm(PrAb ~ DAY + Day2 + Day3 + Day4, family='binomial', data = data_ijk)
        predictions <- predict(poly_ijk, newdata=newdata)
        transf_pred <- 1/(1+exp(-predictions))
        model_results$pol4_peak[counter] <- min(which(transf_pred==max(transf_pred)))
        model_results$pol4_half[counter] <- min(which(transf_pred > max(transf_pred)/2))
        
        poly_ijk <- glm(PrAb ~ DAY + Day2 + Day3 + Day4 + Day5, family='binomial', data = data_ijk)
        predictions <- predict(poly_ijk, newdata=newdata)
        transf_pred <- 1/(1+exp(-predictions))
        model_results$pol5_peak[counter] <- min(which(transf_pred==max(transf_pred)))
        model_results$pol5_half[counter] <- min(which(transf_pred > max(transf_pred)/2))
        
        print(paste(i,j,k))
      }
        
    }
  }
}



Logistic <- function(x){
  list(predictors = list(Asym = 1, xmid = 1, scal = 1),
       variables = list(substitute(x)),
       term = function(predLabels, varLabels){
         paste(predLabels[1], "/(1 + exp((", predLabels[2], "-",
               varLabels[1], ")/", predLabels[3], "))")
       },
       start = function(theta){
         theta[3] <- 1
         theta
       }
  )
}
class(Logistic) <- "nonlin"


counter <- 0
for(i in 1:nyr){
  for(j in 1:ncell){
    data_ij <- dataAll2[which(dataAll2$YEAR==years[i] & dataAll2$cell==cells[j]),]
    for(k in 1:nsp){
      counter <- counter+1
      data_ijk <- data_ij[ , as.numeric(which(colnames(data_ij) %in% c('DAY', species_list[k])))]
      if('X' %in% data_ijk[,2]){
        data_ijk[which(data_ijk[,2] == 'X'),2] <- 1
      }
      data_ijk[,2] <- as.numeric(as.numeric(as.character(data_ijk[,2])) > 0)
      n1 <- sum(data_ijk[,2])
      n0 <- dim(data_ijk)[1] - sum(1-data_ijk[,2])
      n1J <- sum(data_ijk[,2]*as.numeric(data_ijk$DAY > 160))
      if(n1 > 500 & n0 > 500 & n1J > 100){
        colnames(data_ijk)[2] <- "PrAb"
        data_ijk$scaleDay <- scale(data_ijk$DAY)
        scaledLog <- tryCatch(gnm::gnm(PrAb ~ Logistic(scaleDay), data=data_ijk, family='binomial'), error=function(e){})
        if(!is.null(scaledLog)){
          model_results$sLogi_half[counter] <- scaledLog$coefficients[3]*sd(data_ijk$DAY) + mean(data_ijk$DAY)
          print(paste(i,j,k))
        }
      }
      
    }
  }
}

save(model_results, file = '/Users/TingleyLab/Dropbox/Work/Phenomismatch/model_results_nJ100.Rdata')

model_results_clipped <- model_results[which(model_results$n0 > 499 & model_results$n1 > 499 & model_results$n1J > 99),]
peaked_results <- model_results_clipped[which(model_results_clipped$GAM_peak > 50 & model_results_clipped$GAM_peak < 180),]
dim(peaked_results)
# There were 1495 species-year-cell combinations with at least 500 zeros, 500 ones, and 100 ones after 10 June.
# Of these, 1491 showed GAM peaks after jDay 50. 1491 (not a typo, coincidence) showed GAM peaks before the final day of the modeled time window (jDay 213).
# 1454 showed peaks before jDay 180 (end of June).

sensible_sLog <- model_results_clipped[which(model_results_clipped$sLogi_half < 180 & model_results_clipped$sLogi_half > 30),]
dim(sensible_sLog)

sensible_quint_half <- model_results_clipped[which(model_results_clipped$pol5_half < 180 & model_results_clipped$pol5_half > 30),]
dim(sensible_quint_half)
sensible_quint_peak <- model_results_clipped[which(model_results_clipped$pol5_peak < 180 & model_results_clipped$pol5_peak > 30),]
dim(sensible_quint_peak)

sensible_quart_half <- model_results_clipped[which(model_results_clipped$pol4_half < 180 & model_results_clipped$pol4_half > 30),]
dim(sensible_quart_half)
sensible_quart_peak <- model_results_clipped[which(model_results_clipped$pol4_peak < 180 & model_results_clipped$pol4_peak > 30),]
dim(sensible_quart_peak)

sensible_cub_half <- model_results_clipped[which(model_results_clipped$pol3_half < 180 & model_results_clipped$pol3_half > 30),]
dim(sensible_cub_half)
sensible_cub_peak <- model_results_clipped[which(model_results_clipped$pol3_peak < 180 & model_results_clipped$pol3_peak > 30),]
dim(sensible_cub_peak)

sensible_quad_half <- model_results_clipped[which(model_results_clipped$pol2_half < 180 & model_results_clipped$pol2_half > 30),]
dim(sensible_quad_half)
sensible_quad_peak <- model_results_clipped[which(model_results_clipped$pol2_peak < 180 & model_results_clipped$pol2_peak > 30),]
dim(sensible_quad_peak)

# The polynomials consistently yield sensible results.

# Now we want to check which parametric models yield phenological metrics that reliably correlate with the GAM-derived metrics.
plot(model_results_clipped$GAM_half[which(model_results_clipped$sLogi_half>0)], model_results_clipped$sLogi_half[which(model_results_clipped$sLogi_half>0)],
     xlab='GAM prediction', ylab='Scaled Logistic Prediction', main='Half Max: scaled logistic')
abline(1,1)
summary(lm(model_results_clipped$sLogi_half[which(model_results_clipped$sLogi_half>0)] ~ model_results_clipped$GAM_half[which(model_results_clipped$sLogi_half>0)]))
# R-squared = .62, n = 1347

plot(model_results_clipped$GAM_half, model_results_clipped$pol5_half,
     xlab='GAM prediction', ylab='Quintic prediction', main='Half Max: quintic')
abline(1,1)
summary(lm(model_results_clipped$pol5_half[which(model_results_clipped$pol5_half<10^6)] ~ model_results_clipped$GAM_half[which(model_results_clipped$pol5_half<10^6)]))
# R-squared = .98, n = 1487

plot(model_results_clipped$GAM_half, model_results_clipped$pol4_half,
     xlab='GAM prediction', ylab='Quartic prediction', main='Half Max: quartic')
abline(1,1)
summary(lm(model_results_clipped$pol4_half[which(model_results_clipped$pol4_half<10^6)] ~ model_results_clipped$GAM_half[which(model_results_clipped$pol4_half<10^6)]))
# R-squared = .96, n = 1491

plot(model_results_clipped$GAM_half, model_results_clipped$pol3_half,
     xlab='GAM prediction', ylab='Cubic prediction', main='Half Max: cubic')
abline(1,1)
summary(lm(model_results_clipped$pol3_half[which(model_results_clipped$pol3_half<10^6)] ~ model_results_clipped$GAM_half[which(model_results_clipped$pol3_half<10^6)]))
# R-squared = .93, n = 1489

plot(model_results_clipped$GAM_half, model_results_clipped$pol2_half,
     xlab='GAM prediction', ylab='Quadratic prediction', main='Half Max: quadratic')
abline(1,1)
summary(lm(model_results_clipped$pol2_half[which(model_results_clipped$pol2_half<10^6)] ~ model_results_clipped$GAM_half[which(model_results_clipped$pol2_half<10^6)]))
# R-squared = .86, n = 1491


# Peaks
plot(peaked_results$GAM_peak[which(peaked_results$pol5_peak<214 & peaked_results$pol5_peak>50)], peaked_results$pol5_peak[which(peaked_results$pol5_peak<214 & peaked_results$pol5_peak>30)],
     xlab='GAM prediction', ylab='Quintic prediction', main='Peak: quintic')
abline(1,1)
summary(lm(peaked_results$pol5_peak[which(peaked_results$pol5_peak<214 & peaked_results$pol5_peak>50)] ~ peaked_results$GAM_peak[which(peaked_results$pol5_peak<214 & peaked_results$pol5_peak>50)]))
# R-squared = .62, n = 1450
summary(lm(peaked_results$pol4_peak[which(peaked_results$pol4_peak<214 & peaked_results$pol4_peak>50)] ~ peaked_results$GAM_peak[which(peaked_results$pol4_peak<214 & peaked_results$pol4_peak>50)]))
# R-squared = .48
summary(lm(peaked_results$pol3_peak[which(peaked_results$pol3_peak<214 & peaked_results$pol3_peak>50)] ~ peaked_results$GAM_peak[which(peaked_results$pol3_peak<214 & peaked_results$pol3_peak>50)]))
# R-squared = .41
summary(lm(peaked_results$pol2_peak[which(peaked_results$pol2_peak<214 & peaked_results$pol2_peak>50)] ~ peaked_results$GAM_peak[which(peaked_results$pol2_peak<214 & peaked_results$pol2_peak>50)]))
# R-squared = .36

# So far we've looked at correlations between GAM-derived phenophases and parametrically derived phenophases across species-cell-years
# We especially care about our ability to detect an early or late year at a given location.  So we can estimate a model that controls for cell*species, and isolates the ability to predict interannual differences.
rsq::rsq.partial(lm(model_results_clipped$sLogi_half[which(model_results_clipped$sLogi_half>0)] ~ model_results_clipped$GAM_half[which(model_results_clipped$sLogi_half>0)] + model_results_clipped$species[which(model_results_clipped$sLogi_half>0)]*model_results_clipped$cell[which(model_results_clipped$sLogi_half>0)]))$partial.rsq[1]
# partial R-squared: .16
rsq::rsq.partial(lm(model_results_clipped$pol5_half[which(model_results_clipped$pol5_half<10^6)] ~ model_results_clipped$GAM_half[which(model_results_clipped$pol5_half<10^6)] + model_results_clipped$species[which(model_results_clipped$pol5_half<10^6)]*model_results_clipped$cell[which(model_results_clipped$pol5_half<10^6)]))$partial.rsq[1]
# partial R-squared: .87
rsq::rsq.partial(lm(model_results_clipped$pol4_half[which(model_results_clipped$pol4_half<10^6)] ~ model_results_clipped$GAM_half[which(model_results_clipped$pol4_half<10^6)] + model_results_clipped$species[which(model_results_clipped$pol4_half<10^6)]*model_results_clipped$cell[which(model_results_clipped$pol4_half<10^6)]))$partial.rsq[1]
# partial R-squared: .80
rsq::rsq.partial(lm(model_results_clipped$pol3_half[which(model_results_clipped$pol3_half<10^6)] ~ model_results_clipped$GAM_half[which(model_results_clipped$pol3_half<10^6)] + model_results_clipped$species[which(model_results_clipped$pol3_half<10^6)]*model_results_clipped$cell[which(model_results_clipped$pol3_half<10^6)]))$partial.rsq[1]
# partial R-squared: .67
rsq::rsq.partial(lm(model_results_clipped$pol2_half[which(model_results_clipped$pol2_half<10^6)] ~ model_results_clipped$GAM_half[which(model_results_clipped$pol2_half<10^6)] + model_results_clipped$species[which(model_results_clipped$pol2_half<10^6)]*model_results_clipped$cell[which(model_results_clipped$pol2_half<10^6)]))$partial.rsq[1]
# partial R-squared: .48

# We see that the scaled logistic performs poorly; and the polynomials range from ok to quite good, with higher-order polynomials performing better.

# For completeness, we can look at the data for the peak-derived phenophases:
rsq::rsq.partial(lm(peaked_results$pol5_peak[which(peaked_results$pol5_peak<214 & peaked_results$pol5_peak>30)] ~ peaked_results$GAM_peak[which(peaked_results$pol5_peak<214 & peaked_results$pol5_peak>30)] + peaked_results$species[which(peaked_results$pol5_peak<214 & peaked_results$pol5_peak>30)]*peaked_results$cell[which(peaked_results$pol5_peak<214 & peaked_results$pol5_peak>30)]))$partial.rsq[1]
# partial R-squared: .47
rsq::rsq.partial(lm(peaked_results$pol4_peak[which(peaked_results$pol4_peak<214 & peaked_results$pol4_peak>50)] ~ peaked_results$GAM_peak[which(peaked_results$pol4_peak<214 & peaked_results$pol4_peak>30)] + peaked_results$species[which(peaked_results$pol4_peak<214 & peaked_results$pol4_peak>30)]*peaked_results$cell[which(peaked_results$pol4_peak<214 & peaked_results$pol4_peak>30)]))$partial.rsq[1]
# partial R-squared: .33
rsq::rsq.partial(lm(peaked_results$pol3_peak[which(peaked_results$pol3_peak<214 & peaked_results$pol3_peak>50)] ~ peaked_results$GAM_peak[which(peaked_results$pol3_peak<214 & peaked_results$pol3_peak>30)] + peaked_results$species[which(peaked_results$pol3_peak<214 & peaked_results$pol3_peak>30)]*peaked_results$cell[which(peaked_results$pol3_peak<214 & peaked_results$pol3_peak>30)]))$partial.rsq[1]
# partial R-squared: .23
rsq::rsq.partial(lm(peaked_results$pol2_peak[which(peaked_results$pol2_peak<214 & peaked_results$pol2_peak>50)] ~ peaked_results$GAM_peak[which(peaked_results$pol2_peak<214 & peaked_results$pol2_peak>30)] + peaked_results$species[which(peaked_results$pol2_peak<214 & peaked_results$pol2_peak>30)]*peaked_results$cell[which(peaked_results$pol2_peak<214 & peaked_results$pol2_peak>30)]))$partial.rsq[1]
# partial R-squared: .16

# These are generally worse, but essentially all of them are better than the scaled logistic***
# *** the big assumption here is that the GAM is producing good indices for the biological phenophases of interest. That's not necessarily the case, but I think it's safe to assume that where it and the scaled logistic differ, the GAM is more likely provide a useful phenophase.

summary(lm(peaked_results$pol5_half[which(peaked_results$pol5_peak<214 & peaked_results$pol5_peak>50)] ~ peaked_results$GAM_peak[which(peaked_results$pol5_peak<214 & peaked_results$pol5_peak>50)]))
