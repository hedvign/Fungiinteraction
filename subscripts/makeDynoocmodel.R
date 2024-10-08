# code to write col-ext model for each species

# informative or non-informative species

# JAGS


sink(paste0(OUTDIR, "dynocc.txt"))
cat("
# model generated by 'fcn_dynocc9ssvs.R'
model{
### PRIORS ###############################
psi1 ~ dunif(0, 1)            # initial prevalence prior

    ",fill=T)

cat("# Priors alpha intercept colonization ############## \n")
 cat("
    alphacol ~ dnorm(0, 0.7)   
        ",fill=T)
  

# alphaext
  cat("# Priors alpha intercept extinction ##############\n") # #alphaext.prob ~ dunif(0,1)  #alphaext <- logit(alphaext.prob)  
  cat("
    alphaext ~ dnorm(0, 0.7)   
        ",fill=T)

#####################################################
if(whichvars$beta){ #  colonization environment linear #######################
  cat("\n# Priors beta colonization environment linear ##############\n")
   cat("
    pbeta[2] <- 0.2 # Prob. of being in the model
    pbeta[1] <- 1-pbeta[2]
    Precbeta[1] <- 100 # Prior precision if not in the model
    Precbeta[2] <- 0.7 # Prior precision if in the model
    for(k in 1:nstudysites){
      # ssvs:                         # alt.3, ssvs
      Indbeta[k] ~ dcat(pbeta[]) # Indicator for if in the model 
      beta[k] ~ dnorm(0, Precbeta[Indbeta[k]])
    }
      ")
  
}
if(whichvars$beta2){ # colonization environment quadratic #######################################
  cat("# Priors beta colonization environment quadratic ##############\n")
    cat("
    pbeta2[2] <- 0.2 # Prob. of being in the model
    pbeta2[1] <- 1-pbeta2[2]
    Precbeta2[1] <- 100 # Prior precision if not in the model
    Precbeta2[2] <- 0.7 # Prior precision if in the model
    for(k in 1:nstudysites){
      # ssvs:                         # alt.3, ssvs
      Indbeta2[k] ~ dcat(pbeta2[]) # Indicator for if in the model 
      beta2[k] ~ dnorm(0, Precbeta2[Indbeta2[k]])
    }
      ")
  
}


################################################################
if(whichvars$gamma){ # extinction environment linear ###########################################
  cat("# Priors gamma extinction environment linear ##############\n")
    cat("
    pgamma[2] <- 0.2 # Prob. of being in the model
    pgamma[1] <- 1-pgamma[2]
    Precgamma[1] <- 100 # Prior precision if not in the model
    Precgamma[2] <- 0.7 # Prior precision if in the model
    for(k in 1:nstudysites){
      # ssvs:                         # alt.3, ssvs
      Indgamma[k] ~ dcat(pgamma[]) # Indicator for if in the model 
      gamma[k] ~ dnorm(0, Precgamma[Indgamma[k]])
    }
      ")
}
if(whichvars$gamma2){ # extinction environment quadratic ############################
  cat("# Priors gamma extinction environment quadratic ##############\n")
    cat("
    pgamma2[2] <- 0.2 # Prob. of being in the model
    pgamma2[1] <- 1-pgamma2[2]
    Precgamma2[1] <- 100 # Prior precision if not in the model
    Precgamma2[2] <- 0.7 # Prior precision if in the model
    for(k in 1:nstudysites){
      # ssvs:                         # alt.3, ssvs
      Indgamma2[k] ~ dcat(pgamma2[]) # Indicator for if in the model 
      gamma2[k] ~ dnorm(0, Precgamma2[Indgamma2[k]])
    }
      ")
}


################################################################
# The SSVS priors would be: 
# noninformative but narrow dnorm(0, 0.7) for the spike, and informative but wide dnorm(0.7, 100) for the slab.

cat("# Priors delta, colonization species interactions ##############\n")

    cat("
    # SSVS
  for (s in 1:nspcol) {
  pDelta[2,s] <- 0.2 # Prob. of being in the model
  pDelta[1,s] <- 1-pDelta[2,s]
  PrecDelta[1,s] <- 100 # Prior precision if not in the model
  PrecDelta[2,s] <- 0.5 # Prior precision if in the model") # XX
    cat("
  }
 ") 
    # ssvs:
    if(useInformativePriors){
      
   cat("# informative priors and SSVS, colonization
     for(k in 1:nstudysites){ ")
      for (s in 1:nspcol) {
        cat("
          Inddelta[",s,",k] ~ dcat(pDelta[,",s,"])
          MeanDelta[1,",s,"] <- ")
            cat(eval(priorsOtherspCol[s]), "# Prior informative mean if in the model (option 2)") # option 2
            #cat("0 # Prior mean if not in the model (option 3)\n") # option 3
          cat("
          MeanDelta[2,",s,"] <- ")
            cat(eval(priorsOtherspCol[s]), "# Prior informative mean if in the model ")
        cat("
          delta[",s,",k] ~ dnorm(MeanDelta[Inddelta[",s,",k],",s,"], PrecDelta[Inddelta[",s,",k],",s,"])") 
    
      } # end for R
      cat("\n\t}")
    }
    if(!useInformativePriors){
      cat("# no informative priors, SSVS, colonization
          for(k in 1:nstudysites){ ")
      for (s in 1:nspcol) {
        cat("
          Inddelta[",s,",k] ~ dcat(pDelta[,",s,"])    
          delta[",s,",k] ~ dnorm(0, PrecDelta[Inddelta[",s,",k],",s,"])") 
        } # end for R
      cat("\n\t}")
    } # end if not inf priors and ssvs

  
  
cat("\n
# missing values in other species
  for (i in 1:nsites) {
    for (t in 1:nyears) {
      for (s in 1:nspcol) {
        # variable is 0/1 sp pres/abs, so should be Bernouilli distribution
        otherspdatcol[i, s, t] ~ dbern(0.5) # xx
      }
    }
  }
  ")




if(nspext > 0 ){ # epsilon, extinction species interactions ##############################
  cat("# Priors epsilon, extinction species interactions ##############\n")

    cat("
    # SSVS
  for (s in 1:nspext) {
  pEpsilon[2,s] <- 0.2 # Prob. of being in the model
  pEpsilon[1,s] <- 1-pEpsilon[2,s]
  PrecEpsilon[1,s] <- 100 # Prior precision if not in the model
  PrecEpsilon[2,s] <- 0.7 # Prior precision if in the model")
    cat("
  }
 ") 
    # ssvs:
    if(useInformativePriors){
      
      cat("# informative priors and SSVS, extinction
     for(k in 1:nstudysites){ ")
      for (s in 1:nspext) {
        cat("
          Indepsilon[",s,",k] ~ dcat(pEpsilon[,",s,"])
          epsilon[",s,",k] ~ dnorm(")
        cat(eval(priorsOtherspExt[s]))
        cat(",")
        cat(paste0("PrecEpsilon[Indepsilon[",s,",k],",s,"]"))
        cat(")")
      } # end for R
      cat("\n\t}")
    }
    if(!useInformativePriors){
      cat("# no informative priors, SSVS, extinction
          for(k in 1:nstudysites){ ")
      for (s in 1:nspext) {
        cat("
          Indepsilon[",s,",k] ~ dcat(pEpsilon[,",s,"])    
          epsilon[",s,",k] ~ dnorm(0, PrecEpsilon[Indepsilon[",s,",k],",s,"])") 
      } # end for R
      cat("\n\t}")
    } # end if not inf priors and ssvs
  }
  
  
  cat("\n
# missing values in other species
  for (i in 1:nsites) {
    for (t in 1:nyears) {
      for (s in 1:nspext) {
        # variable is 0/1 sp pres/abs, so should be Bernouilli distribution
        otherspdatext[i, s, t] ~ dbern(0.5) # xx
      }
    }
  }
  ")




###################################################
### decay 
if(whichvars$betadecay){ # col
  cat("# Priors beta colonization decay linear ##############\n")
  cat("
#####################################################
  pbetadecay[2] <- 0.2 # Prob. of being in the model
  pbetadecay[1] <- 1-pbetadecay[2]
  Precbetadecay[1] <- 100 # Prior precision if not in the model
  Precbetadecay[2] <- 0.7 # Prior precision if in the model

  for(k in 1:nstudysites){ 
    # ssvs:
    Indbetadecay[k] ~ dcat(pbetadecay[]) # Indicator for if in the model
    betadecay[k] ~ dnorm(0, Precbetadecay[Indbetadecay[k]])
  }")
  
}

if(whichvars$betadecay2){ # col
  cat("# Priors beta colonization decay quadratic ############## \n")
    cat("
  pbetadecay2[2] <- 0.2 # Prob. of being in the model
  pbetadecay2[1] <- 1-pbetadecay2[2]
  Precbetadecay2[1] <- 100 # Prior precision if not in the model
  Precbetadecay2[2] <- 0.7 # Prior precision if in the model

  for(k in 1:nstudysites){ 
    # ssvs:
    Indbetadecay2[k] ~ dcat(pbetadecay2[]) # Indicator for if in the model
    betadecay2[k] ~ dnorm(0, Precbetadecay2[Indbetadecay2[k]])
  }")
  
}


### decay ext
if(whichvars$gammadecay){ # ext
  cat("# Priors gamma extinction decay linear ##############\n")
 
    cat("
  pgammadecay[2] <- 0.2 # Prob. of being in the model
  pgammadecay[1] <- 1-pgammadecay[2]
  Precgammadecay[1] <- 100 # Prior precision if not in the model
  Precgammadecay[2] <- 0.7 # Prior precision if in the model

  for(k in 1:nstudysites){ 
    # ssvs:
    Indgammadecay[k] ~ dcat(pgammadecay[]) # Indicator for if in the model
    gammadecay[k] ~ dnorm(0, Precgammadecay[Indgammadecay[k]])
  }")
  
}

if(whichvars$gammadecay2){ # col
  cat("# Priors gamma extinction decay quadratic ##############\n")
    cat("
  pgammadecay2[2] <- 0.2 # Prob. of being in the model
  pgammadecay2[1] <- 1-pgammadecay2[2]
  Precgammadecay2[1] <- 100 # Prior precision if not in the model
  Precgammadecay2[2] <- 0.7 # Prior precision if in the model

  for(k in 1:nstudysites){ 
    # ssvs:
    Indgammadecay2[k] ~ dcat(pgammadecay2[]) # Indicator for if in the model
    gammadecay2[k] ~ dnorm(0, Precgammadecay2[Indgammadecay2[k]])
  }")
  
}

###############################################
###############################################
cat("
###############################################
# Linear models
for (i in 1:nsites) {
  # Initial conditions of system
   z[i, 1] ~ dbern(psi1) # replace y with z here for detection XX staticmodel
  for (t in 1:(nyears)) {", fill = T) # XX


# ext ################################################
cat("        ep[i, t] <- ilogit(")

if(nstudysites == 1){
  cat("alphaext\n\t")
}
if(nstudysites > 1){
  cat("alphaext[studysites[i]]\n\t")
}

if(whichvars$gamma) {cat("+ (gamma[studysites[i]] * envext[i,1])\n\t") }
if(whichvars$gamma2) {cat("+ (gamma2[studysites[i]] * envext[i,2])\n\t") }
if(whichvars$gammadecay) {cat("+ (gammadecay[studysites[i]] * decay[i, t])\n\t") } 
if(whichvars$gammadecay2) {cat("+ (gammadecay2[studysites[i]] * decay2[i, t])\n\t") }   
if(ncol(otherspdatext) > 0 ){ # sp in col
  cat("+ inprod(epsilon[,studysites[i]], otherspdatext[i, , t])   \n\t") 
}
cat(")\n") # end ilogit 


# start cps ##########################
cat("        cp[i, t] <- ilogit(") 

if(nstudysites == 1){
  cat("alphacol\n\t")
}
if(nstudysites > 1){
  cat("alphacol[studysites[i]]\n\t")
}
if(whichvars$beta) {cat("\t+ (beta[studysites[i]] * envcol[i,1]) \n\t") }
if(whichvars$beta2) {cat("+ (beta2[studysites[i]] * envcol[i,2]) \n\t") }
if(whichvars$betadecay) {cat("+ (betadecay[studysites[i]] * decay[i, t])  \n\t") } 
if(whichvars$betadecay2) {cat("+ (betadecay2[studysites[i]] * decay2[i, t])  \n\t") }   
if(ncol(otherspdatcol) > 0 ){ # sp in col
  cat("+ inprod(delta[,studysites[i]], otherspdatcol[i, , t])\n") 
}
cat(")\n") # end ilogit 
#cat("col test end\n")

cat("}  # end all years
}  # end all sites
 ",fill=T)
cat("############################
# State transitions (replace y with z here for detection)
  for (t in 2:nyears) { # XX 2:nyears # static model
	  for (i in 1:nsites) {
	       z[i, t] ~ dbern(z[i,t - 1] * ep[i,t - 1]  + ((1 - z[i,t - 1]) * cp[i,t - 1]))
       # z[i, t] ~ dbern(cp[i,t]) # xx static model
		}
  }")
cat("
  detpar ~ dnorm(0, 1) # detection probability prior
  for (i in 1:ndet) { 
  	  det[i] ~ dbern(detmu[i]) 
  	  logit(detmu[i]) <- detpar   
    }
  pdet <- exp(detpar) / (1 + exp(detpar)) # pdet <- ilogit(detpar)
  # Observation model
  for (i in 1:nsites) {
    for (t in 1:nyears) {
      y[i, t] ~ dbern(zdet[i,t])
      zdet[i,t] <- z[i, t] * pdet
      dist[i, t] ~ dbern(zdet[i,t]) # saved for counting states
    }
  }
  for (t in 2:nyears){
    for (i in 1:nsites){
      sumcols[i,t] <- equals((dist[i,t-1] - dist[i,t]),  -1)
      sumexts[i,t] <- equals((dist[i,t-1] - dist[i,t]), 1)
      tmpncol[i,t] <- sum(sumcols[i,t]) # counts nCol
      tmpnext[i,t] <- sum(sumexts[i,t]) # counts nExt
      tmpnocc[i,t] <- sum(dist[i,t]) # counts Occurrences
      #counts coccurr of sp and all other sp in t2
    }
   noccs[t] <- sum(tmpnocc[,t])
   ncols[t] <- sum(tmpncol[,t])
   nexts[t] <- sum(tmpnext[,t])
    }
   pcol <- mean(cp[,1])
  pext <- 1-mean(ep[,1]) 
  
  ")

cat("
} # end model \n") # end model

##################

sink()

#file.copy(from = list.files(SIMDIR, full.names = T), to = OUTDIR, overwrite = T) 
