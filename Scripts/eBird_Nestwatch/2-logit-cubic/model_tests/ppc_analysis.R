#### I'm not positive this next part works, since I haven't had enough computer power/time
# to fit the model for a bunch of species-cell-years and check. The basic idea here is we want
# to histogram the posterior predictive check values. We'll look at all values, the values near
# the estimated halfmax.

halfmaxppc <- list()

for(f in 1:length(formulas)){
  ppcsf <- unlist(ppcs[[f]], recursive = T)  # If we get the posterior predictive data across multiple species into a list, then we can run this code over that object to histogram the ppc values for every bin for every species, which could be useful for comparing polynomials overall.
  hist(ppcsf, main = paste("formula", f))
  
  # Now pick out the bin that contains the (posterior mean) halfmax. This is an inefficient way 
  # to do it, but it's reliable and I'm confident it works (as long as I didn't screw up)!
  halfmaxppc[[f]] <- vector()
  counter <- 0
  for (j in 1:nyr){
    for (k in 1:ncel){
      counter <- counter + 1
      halfmaxbin <- 1 + mean(halfmax_matrix_list[[f]][[j]][[k,]]) %/% 20
      halfmaxppc[[f]][counter] <- ppcs[[f]][[j]][[k]][[halfmaxbin]]
    }
  }
}