# make a model output file with all info 
# and an object?

##################################################################
# does the parameter 95% estimates overlap zero
# both intervals have to be true
library(MCMCvis)



post = post[,(!colnames(post[[1]]) %in% "deviance")]
     
quants = summary(post[,], quantiles = c(.025,0.05,0.95,0.975))[[2]]
quants
overlap095 = quants[,"2.5%"] <= 0 &  0 <= quants[,"97.5%"]
overlap090 = quants[,"5%"] <= 0 &  0 <= quants[,"95%"]
overlap090
overlap095
overlap095[overlap095 == 0] # these do not overlap zero and are important
overlap095[overlap095 == 1] # these overlap zero and are not important

# calculate Rhat with mcmcvis package
summ = MCMCsummary(post)
summ[,c("Rhat")]
summ


# calculate deviance
pD = round(out$pD, 3)
DIC = round(out$DIC, 3)
############################################
# MC error

mcmcerrsSDpc = getMCerror(post, n.chains=2, SDpc = T)

mcmcerrsSDpcbelow5 = mcmcerrsSDpcbelow1.5 = rep(F, times = length(mcmcerrsSDpc))
mcmcerrsSDpcbelow1.5
mcmcerrsSDpc < 5
mcmcerrsSDpcbelow5[mcmcerrsSDpc < 5] = T
mcmcerrsSDpcbelow1.5[mcmcerrsSDpc < 1.5] = T
mcmcerrsSDpcbelow1.5
mcmcerrsSDpcbelow5

# combine measures ##################################################
summary_dev = cbind(summ, overlap090, overlap095, 
                    MCerr = mcmcerrsSDpc, MCerrbelow5 = mcmcerrsSDpcbelow5, MCerrbelow1.5 = mcmcerrsSDpcbelow1.5)

# save as Robject
summary_devs = list(summary_dev = summary_dev, pD = pD, DIC = DIC) # 

####################################################
  #table(studysites)
  #studysiteschara
#  unlist(studysitestoanalyse)
 # studysitesorder = unlist(studysitestoanalyse) #unique(envinfo[,"ss"])[unique(envinfo[,"ss"]) %in% unlist(studysitestoanalyse)]
  #print(studysitesorder)
  # 
  # dynsPersite = data.frame(matrix(NA, ncol = 7, nrow = 0))
  # 
  # dynsPersite
  # for(sxxx in studysitesorder){
  #   print(sxxx)
  #   tabb = (t(as.matrix(spinfoSites[[sxxx]][sp,1:7])))
  #   (tabb)
  #   dynsPersite = rbind(dynsPersite, as.numeric(tabb))
  #   colnames(dynsPersite) = colnames(tabb)
  # }
  # dynsPersite
  # StudySite = paste(studysitesorder, 1:nrow(dynsPersite))
  # rownames(dynsPersite) = StudySite 
  # dynsPersite = cbind(StudySite, dynsPersite)
  # dynsPersite
  # 
  # print(dynsPersite)
  
print(summary_dev)
print(tail(summary_dev))

