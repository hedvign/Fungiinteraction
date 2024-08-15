# Models for: "Combining observational and experimental data to estimate 
# environmental and species drivers of fungal metacommunity dynamics" 
# Hedvig Nenzen, SLU Artdatabanken, Sweden, 2023.
# works with static envrionment covarariates and dynamic decay covariate
# based on Kery & Royle 2020_chp 4, p 219
# hierarchical model based on Kery & Royle 2020_chp 3, p 180.
rm(list = ls())
library(jagsUI)
library(rjags) 
library(ggplot2)
library(bayesplot)
library(MCMCvis)
library(wiqid)

setwd("~/Documents/Polypores/MODELfungiinteractions")


useInformativePriors = F # T or F

# use data with only 'live' fruit bodies, or dead and live fruitbodies: 'dead'
mainfoldername = "livefruitbodies_InfPriors/"
#mainfoldername = "livefruitbodies_noInfPriors/"

#mainfoldername = "deadfruitbodies_InfPriors/"
#mainfoldername = "deadfruitbodies_noInfPriors/"

dir.create(mainfoldername, showWarnings = F)

########################################
# create matrix to indicate which species and variables should be tested
spnames = c( "fomipini","heteparv", "tricabie", # early successional species
             "phlecent","phelferr", "fomirose", "postcaes", "antrseri", # mid sp
             "phelviti", "gloesepi","skelbrev", "phelnigr") # late

varnames = c("beta", "beta2", "betadecay", "betadecay2" # col 
             , "gamma", "gamma2", "gammadecay", "gammadecay2" # ext
)
envvarss = matrix(1, ncol = length(varnames), nrow = length(spnames), dimnames = list(spnames, varnames))
colvarss = matrix(1, ncol = length(spnames), nrow = length(spnames), dimnames = list(spnames, paste0("col.", spnames)))
extvarss = matrix(1, ncol = length(spnames), nrow = length(spnames), dimnames = list(spnames, paste0("ext.", spnames)))
spvarss = matrix(1, ncol = 2, nrow = length(spnames), dimnames = list(spnames, c("speciesCol",	"speciesExt")))

diag(colvarss) = diag(extvarss) = 0 # not test itself
whichvarsModels = as.data.frame(cbind(focalsp = spnames, envvarss, spvarss, colvarss, extvarss)) 
rownames(whichvarsModels) = 1:nrow(whichvarsModels) 
whichvarsModels
rm(envvarss, colvarss, extvarss, spvarss)
write.table(whichvarsModels, file = paste0(mainfoldername,"whichvarsModels.csv"), sep = ";", row.names = F)


#####################################################################m#########################
# RUN MODELS ##################################################################################
modelx = 12 ####################################################################################
modelx
foldername = mainfoldername
#for(modelx in 1:nrow(whichvarsModels)){ #nrow(whichvarsModels)################################
fcnRunmodel = function(modelx, foldername){ 
  
  envnames = colparms = extparms = c("DIAM", "DIAM2") 
  setwd("~/Documents/Polypores/MODELfungiinteractions")
  
  whichvarsModels = read.table(paste0(foldername, "whichvarsModels.csv"), sep = ";", header = T)
  spnames = unique(whichvarsModels$focalsp)
  
  ni <- 2000
  nc <- 3 # n chains
  nb = ifelse(ni==10000, 2000, floor(ni/2)) # n burn-i
  nt = 5 # thin
  
  whichvars = as.data.frame(whichvarsModels[modelx,])
  simname = paste0(foldername, "data_sp", modelx)
  simname
  
  sp = as.character(whichvars["focalsp"]) 
  othersp = spnames[!spnames %in% sp]
  othersp
  # open data #############################################
  if(F){
    if(grepl(foldername, pattern = "live")){
      load('DATA/ColExtdata_arrays_livefbs.Rdata')
    }
    if(grepl(foldername, pattern = "dead")){
      load('DATA/ColExtdata_arrays_deadfbs.Rdata')
    }
    envdata = envtable
    decaydat = as.matrix(decayarray[, "decay", ]) # remove intercept, so only matrix
    rm(decayarray)
    # make new variable, DIAM2.
    envdata = cbind(envdata, DIAM2 = (envdata$DIAM))
    
    # remove sp not modelled
    sparray = sparray[,spnames,]
    dim(sparray)
    apply(sparray, MARGIN = c(2, 3), sum)
    table(rowSums(sparray[,,1], na.rm=T))
    
    #############################
    # subselect to some patches ######################################
    # select two study sites
    studysites = envinfo$ss # vector of study sites ID
    studysitestoanalyse = list("Rorstrand", "Gardfjallet")
    
    selectsites = studysites = rep(F, dim(sparray)[1])
    for(sxxx in studysitestoanalyse){
      selectsites[envinfo[,"ss"] == sxxx] = T
    }
    selectsites
    table(selectsites)
    studysiteschara = envinfo[selectsites, "ss"]
    studysitestoanalyse
    factor(studysiteschara, levels = unlist(studysitestoanalyse))
    studysites = as.numeric(factor(studysiteschara, levels = unlist(studysitestoanalyse)))
    
    
    sparray = sparray[selectsites, ,,drop=F]
    envdata = envdata[selectsites,,drop=F]
    decaydat = decaydat[selectsites,, drop = F]
    length(studysites)
    unique(studysites)
    
    # remove years that have only NAs
    nrow(sparray[,1,])
    colSums(apply(is.na(sparray[,,]), FUN = sum, MARGIN = c(2,3)) != 0)
    yearskeep = colSums(is.na(sparray[,1,])) != nrow(sparray[,1,]) # for the first sp...
    yearskeep
    sparray = sparray[, ,yearskeep]
    decaydat = decaydat[,yearskeep]
    
    apply(sparray, MARGIN = c(2,3), FUN = sum)
    dimnames(sparray)
    save(sparray, envdata, decaydat, file = 'DATA/sparray_livefbs_RorstrandGardfjället.Rdata') # envinfo, studysiteschara, 
    
    write.table(cbind(rownames(envdata), envdata[,c("DIAM", "DIAM2")]), file = "DATA/envdata.csv", sep = ",", row.names = F)
    write.table(cbind(rownames(envdata), decaydat), file = "DATA/decaydata.csv", sep = ",", row.names = F)
    
    if(grepl(foldername, pattern = "live")){
      write.table(cbind(rownames(sparray[,,1]), sparray[,,1]), file = "DATA/sparray_t1_livefbs.csv", sep = ",", row.names = F)
      write.table(cbind(rownames(sparray[,,2]), sparray[,,2]), file = "DATA/sparray_t2_livefbs.csv", sep = ",", row.names = F)
      
    }  
    if(grepl(foldername, pattern = "dead")){
      write.table(cbind(rownames(sparray[,,1]), sparray[,,1]), file = "DATA/sparray_t1_deadfbs.csv", sep = ",", row.names = F)
      write.table(cbind(rownames(sparray[,,2]), sparray[,,2]), file = "DATA/sparray_t2_deadfbs.csv", sep = ",", row.names = F)
    }
  } # end subselectt data
  #############################
  
  if(grepl(foldername, pattern = "live")){
    sparray_t1 = read.table("DATA/sparray_t1_livefbs.csv", sep = ",", row.names = 1, header = T)
    sparray_t2 = read.table("DATA/sparray_t2_livefbs.csv", sep = ",", row.names = 1, header = T)
  }
  
  if(grepl(foldername, pattern = "dead")){
    sparray_t1 = read.table("DATA/sparray_t1_deadfbs.csv", sep = ",", row.names = 1, header = T)
    sparray_t2 = read.table("DATA/sparray_t2_deadfbs.csv", sep = ",", row.names = 1, header = T)
  }
  
  
  sparray = array(NA, dim = list(nrow(sparray_t1), ncol(sparray_t1), 2))
  dimnames(sparray) = list(rownames(sparray_t1), colnames(sparray_t1), c("t1", "t7"))
  sparray[,,1] = as.matrix(sparray_t1)
  sparray[,,2] = as.matrix(sparray_t2)
  
  envdata = read.table("DATA/envdata.csv", sep = ",", row.names = 1, header = T)
  decaydat = read.table("DATA/decaydata.csv", sep = ",", row.names = 1, header = T)
  decaydat2 = decaydat
  
  
  OUTDIR = paste0(simname, "/")
  dir.create(OUTDIR, showWarnings = F)
  
  
  # make species-specific data ##########################################
  sparray[, sp, ]
  spdat = sparray[, sp, , drop = T]
  spdat = matrix(spdat, ncol = dim(sparray)[[3]], dimnames = list(NULL, dimnames(sparray)[[3]]))
  spdat
  
  otherspdat = sparray[, othersp, , drop = T]# array 
  otherspdat
  
  # If random or hierarchical effects
  nstudysites = 1
  studysites = rep(1, dim(sparray)[1]) # no hierarchical effects, all sites equal
  
  spnames
  
  # save orig values, and dynamics in long format
  source("subscripts/fcnMakelongdat.R")
  alldatallsporig = makelongdat(sparray = sparray, spnames = spnames, othersp = othersp,
                                colparms = colparms, extparms = extparms, 
                                envdata = envdata, decaydat = decaydat, decaydat2 = decaydat^2)
  
  ###############################################################
  # standardize each var by subtract mean and divide by standard deviation
  source("subscripts/fcnScale.R")
  
  # first standardize, then square
  envx="DIAM2"
  head(envdata)
  for (envx in envnames) {
    res = fcnScale(origvar = envdata[, envx], logvar = F) #  only scale
    envdata[, envx] = res[["scaledvar"]] # replace original  
    
  }
  head(envdata)
  #scale decaydata, which is separate 3D object
  decaydat[,1] = fcnScale(origvar = decaydat[,1], logvar = F)[["scaledvar"]]
  decaydat[,2] = fcnScale(origvar = decaydat[,1], logvar = F)[["scaledvar"]]
  
  # Square values ###########################
  envdata$DIAM2 = envdata$DIAM2^2
  decaydat2 = decaydat * decaydat
  
  # save transformed env values for all sp #####
  alldatallsptrans = makelongdat(sparray = sparray, spnames = spnames, othersp = othersp,
                                 colparms = colparms, extparms = extparms, 
                                 envdata = envdata, decaydat = decaydat, decaydat2 = decaydat2)
  save(alldatallsporig, alldatallsptrans, allspnamesorig = spnames, envdata, # save long data files for plotting later
       file = paste0(OUTDIR, "alldatallsp", sp,".RData")) 
  
  # make col and ext env, subset to selected input variables
  # constant climate (but decay is time-specific)
  envcol = as.matrix(envdata[, colparms])
  envext = as.matrix(envdata[, extparms])
  colnames(envcol) = colparms
  colnames(envext) = extparms
  
  decay = decaydat # decaydat not needed to subset for sp, same for allsp
  decay2 = decaydat2
  y = spdat
  
  
  ############################
  # detection data:
  det = c(0,1,1,1,1,0,1,1,0,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,1,1,1,1,0,1,1,1,1,1,1,1,0,0,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1)
  if(sp == "tricabie") det = rep(1,13)
  if(sp == "fomipini") det = rep(1,14)
  if(sp == "antrseri") det = c(1,1,0,1,1,1,1,1,0,1,1,1)
  ndet = length(det)
  
  #####################################
  # make initial values
  # All other inits zero. Then start second chain at sp prev + 0.001.
  cloglog <- function(theta){log(-log(1 - theta))}
  logit <- function(theta){log(theta / (1-theta))}
  invlogit <- function(theta){exp(theta) / (1 + exp(theta))}  # note: inbuilt plogis() does the same
  
  
  extprobb = 1-mean(alldatallsporig[[sp]]$ext[,"t2"])
  colprobb = mean(alldatallsporig[[sp]]$col[,"t2"])
  
  # default initial values:
  if(colprobb == 0 | colprobb == 1 | is.nan(colprobb)| is.na(colprobb)){
    colprobb = as.numeric(spinfoSites[["all"]][sp,"pCol"]) 
  }
  if(extprobb == 0 | extprobb == 1 | is.nan(extprobb) | is.na(extprobb)) {
    extprobb = as.numeric(spinfoSites[["all"]][sp,"pExt"]) #0.01
  }
  
  colprobb
  extprobb
  
  colinit = logit(colprobb) 
  extinit = logit(extprobb) 
  inits <- list(list(alphacol.prob=colinit, alphaext.prob=extinit), 
                list(alphacol.prob=colinit+0.0001, alphaext.prob=extinit+0.0001), 
                list(alphacol.prob=colinit-0.0001, alphaext.prob=extinit-0.0001) # 3 chains
  )
  inits
  length(inits)
  
  ################################################
  # get names of interacting species
  dim(otherspdat)
  
  colspnames = colnames(whichvars[,paste0("col.",othersp)] == 1)[whichvars[,paste0("col.",othersp)] == 1]
  colspnames
  colspnamesshort = substr(colspnames, 5, 12)
  
  extspnames = colnames(whichvars[,paste0("ext.",othersp)] == 1)[whichvars[,paste0("ext.",othersp)] == 1]
  extspnamesshort = substr(extspnames, 5, 12)
  extspnames
  extspnamesshort
  
  print(paste0(spnames[modelx], ", ", modelx, " of ", length(spnames), " species (nOcc = ", sum(colSums(spdat, na.rm=T), na.rm=T), ")",sep = " "))
  nspcol = dim(otherspdat)[[2]]
  nspext = dim(otherspdat)[[2]]
  nspcol
  
  nsites = dim(spdat)[1]
  nyears = dim(spdat)[2]
  
  # give the model also Z as data, with 0s as NAs:
  z <- y
  z[z == 0] <- NA 
  
  ################################################
  # get informative priors from experiments 
  priorsOtherspCol = priorsOtherspPrecisionCol = matrix(0, ncol = length(colspnamesshort), nrow = 1, dimnames = list(NULL, colspnamesshort))
  priorsOtherspExt = priorsOtherspPrecisionExt = matrix(0, ncol = length(extspnamesshort), nrow = 1, dimnames = list(NULL, extspnamesshort))
  
  priorsOtherspPrecisionCol[] = priorsOtherspPrecisionExt[] = 0.7 # default
  priorsOtherspCol
  
  if(useInformativePriors){
    load("DATA/finalInformativePriors.Rdata") 
    priorsOtherspCol = priorsOtherspCol[sp,colspnamesshort]
    priorsOtherspPrecisionCol = priorsOtherspPrecisionCol[sp,colspnamesshort]
    priorsOtherspExt = priorsOtherspExt[sp,extspnamesshort]
    priorsOtherspPrecisionExt = priorsOtherspPrecisionExt[sp,extspnamesshort]
  }
  
  ######################################################################
  bundledata <- list(y = y, z = z, nsites = nsites, nyears = nyears, nstudysites = nstudysites, studysites = studysites,
                     envcol = envcol, envext = envext, 
                     decay = decay, decay2 = decay2, 
                     priorsOtherspCol = priorsOtherspCol, priorsOtherspPrecisionCol = priorsOtherspPrecisionCol,
                     priorsOtherspExt = priorsOtherspExt, priorsOtherspPrecisionExt = priorsOtherspPrecisionExt,
                     det = det, ndet = ndet)
  
  # Parameters monitored:
  params <- c("alphacol", "alphaext", 
              "pcol", "pext" , "pdet", "ncols", "nexts", "noccs", colnames(whichvars)[2:9])
  params
  
  if(nspcol > 0){
    params = c(params, "delta")
    otherspdatcol = otherspdat
    bundledata2 <- list(otherspdatcol = otherspdatcol, nspcol = nspcol)
    bundledata = append(bundledata, bundledata2)         
    names(bundledata)
  }
  if(nspext > 0){ # sp in ext
    params = c(params, "epsilon")
    otherspdatext = otherspdat
    bundledata2 <- list(otherspdatext = otherspdatext, nspext = nspext)
    bundledata = append(bundledata, bundledata2)   
  } 
  names(bundledata)
  
  params = c(params, "Indbeta", "Indbeta2","IndbetaNspt1", "IndbetaQuad", "Indbetadecay", "Indbetadecay2","IndbetadecayQuad",
             "Indgamma", "Indgamma2", "IndgammaNspt1", "IndgammaQuad", "Indgammadecay", "Indgammadecay2","IndgammadecayQuad", 
             "Inddelta", "Indepsilon") 
  
  # write species-specific model (with its informative priors)
  source(paste0("subscripts/makeDynoocmodel.R"), local = T)
  
  # Run model ###################################################
  ###############################################################
  out <- jags(bundledata, parameters.to.save = params, inits = inits, 
              model.file = paste0(OUTDIR, "dynocc.txt"), parallel = T,
              n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, verbose = T, store.data = TRUE)
  out
  print("ran model")
  outmc = as.mcmc(out)
  post = outmc$samples
  
  # save model results in txt #############################################
  sink(file = paste(OUTDIR, "ModelSummary", sp, ".txt", sep = ""))
  print(out)
  print(whichvars)
  print(paste0(sp,sep = " dynamics: "))
  print(paste0("colprob (init) = ", round(colprobb, 4)))
  print(paste0("extprob (init) = ", round(extprobb, 4)))
  if(nspcol > 0) print(paste0(colspnames,sep = " othersp "))
  if(nspext > 0) print(paste0(extspnames,sep = " othersp "))
  source("subscripts/calc_summary.R", local = T) ##################################
  sink() ################################################
  # plots
  params = colnames(post[[1]]) # remove some params from plots
  params = params[!params %in% c("deviance", "pcol", "pext", "pdet")]
  params = params[!grepl(params, pattern = "ncols")]
  params = params[!grepl(params, pattern = "nexts")]
  params = params[!grepl(params, pattern = "noccs")]
  params = params[!grepl(pattern = "^Ind", x = params)]
  color_scheme_set(scheme = "brewer-Spectral")
  
  nparas = 1:length(params)
  parasets = nparas[(nparas %% 6) == 1]
  if(length(nparas)>6){
    for(paraset in parasets){
      someparas = paraset:(paraset+5) 
      paranamess  = params[someparas]
      paranamess = paranamess[!is.na(paranamess)]
      
      jpeg(paste0(OUTDIR, "traceplots_", paraset, ".jpeg"), width = 500, height = 400, units = "px")
      print(mcmc_trace(post, pars = paranamess))
      dev.off()
    }
  }
  jpeg(paste(OUTDIR, "histplots.jpeg",sep=""), width = 300, height = 500, units = "px")
  print(mcmc_areas(post, pars = params, prob = 0.7, # shaded region
                   prob_outer = 1, point_est = "median")
        + labs(title = sp, subtitle = "Medians, 90% and 100% intervals")) 
  
  dev.off()
  jpeg(paste(OUTDIR, "histplots_delta.jpeg",sep=""), width = 500, height = 800, units = "px")
  print(mcmc_areas(post, regex_pars = "^delta", prob = 0.7, # shaded region
                   prob_outer = 1, point_est = "median")
        + labs(title = sp, subtitle = "Medians, 90% and 100% intervals")) 
  dev.off()
  
  
  # save results and final data set for plotting:
  save(post, alldatallsporig, alldatallsptrans,
       file = paste(OUTDIR, sp, "_model.RData", sep = ""))
  
  print(OUTDIR)
  print("end sp")
  print("##################################")
  
} ########################################## end model
mainfoldername
nmodels <- nrow(whichvarsModels)
nmodels

# run one model:
#fcnRunmodel(modelx, mainfoldername)

# run model for each species in loop
for(modelx in 1:nmodels){
  fcnRunmodel(modelx, mainfoldername)
}

#######################################################
# to parallelize
if(F){ 
  #library("Rmpi")
  #library("parallel")
  #library("snow")
  nnodes = 19
  cl <- makeCluster(nnodes, type = "MPI")
  # run the parallel simulations
  clusterMap(cl, fcnRunmodel, modelx = 1:nmodels, foldername = foldername)
  stopCluster(cl)
}
print(mainfoldername)
