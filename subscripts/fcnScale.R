# 'standardize' each var by subtract mean and divide by standard deviation
fcnScale = function (origvar, logvar){
  origvar[origvar == 0] = 0.00001 # if zero add tiny amount
  if(logvar == T){
    origvar = log(origvar)
  } 
  COVmax <- max(origvar, na.rm = T)  # get min max first
  COVmin <- min(origvar, na.rm = T)
  COVmeans <- mean(origvar, na.rm = T)  # get mean, sd
  COVsds <- sd(origvar, na.rm = T) 
  # then scale standardising the covs (needed for min max (on standardised scale!))
  scaledvar <- (origvar - COVmeans)/COVsds  
  return(list(scaledvar = scaledvar, 
              COVmax, COVmin, COVmeans, COVsds))
}

# if sd 1 ( = variance one too)
# standard deviation is the square root of the variance
# vals = rnorm(0, sd = .5, n = 1000)
# sd(vals)
# var(vals)
# (0.5^2) # sd^2 = var
# hist(vals)
# so be careful if want to change to 2 sds priors! Gelman