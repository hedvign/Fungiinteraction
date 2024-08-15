# compares plots
rm(list = ls())
setwd("~/Documents/Polypores/MODELfungiinteractions")

library(mcmcplots)
library(bayesplot)
library(gtools)
library(corrplot)
library(RColorBrewer)
library(corrplot)

# compare informative and uninformative priors: 
OUTDIR = "livefruitbodies_InfPriors/" # Informative priors, OUTDIR
#OUTDIR = "livefruitbodies_noInfPriors/" # Uninformative priors, OUTDIR2

list.files(OUTDIR)


whichvarsModels = read.table(paste0(OUTDIR, "whichvarsModels.csv"), sep = ";", header =T)
whichvarsModels

# read spnames
spnames = c( "fomipini","heteparv",  "tricabie", # early successional species
             "phlecent","phelferr", "fomirose", "postcaes", "antrseri", # mid sp
             "phelviti", "gloesepi","skelbrev", "phelnigr") # late
# taxonomically updated species names: 
spnameslong = c("Fomitopsis pinicola", "Heterobasidion parviporum", "Trichaptum abietinum", 
                "Phlebia centrifuga", "Phellinidium ferrugineofuscum", "Fomitopsis rosea", "Postia caesia", "Neoantrodia serialis", 
                "Fuscoporia viticola", "Gloeophyllum sepiarium", "Skeletocutis brevispora", "Phellopilus nigrolimitatus")
# no genus
nams = strsplit(spnameslong, split = " ")
genus = lapply(nams, `[[`, 1)
speccs = lapply(nams, `[[`, 2)
genus = paste0(substr(unlist(genus), 1,1), ".")
spnamesmediumlong = paste(genus, speccs, sep = "~")
spnamesmediumlongspace = paste(genus, speccs, sep = " ")
spnamesmediumlong
spnamesmediumlong[spnamesmediumlong == "P.~ferrugineofuscum"] = "P.~ferrug."

processes = c("col", "ext")
processeslong = c("colonization", "extinction")
names(processeslong) = processes


modelnames = list.files(OUTDIR, pattern = paste0("^data")) 
modelnames = mixedsort(modelnames)
rownames(whichvarsModels) = modelnames

# make cols ######################

cols <- rev(brewer.pal(n = length(spnames), name = "Set3"))
names(cols) = spnames
cols
#sequence of hexadecimal characters giving (alternatively) the length of the segments and the blanks separating them.
ltys = rep(c("22", "4212","62", "8212"), 3)[1:length(spnames)]

VariancePartition = array(NA, dim = list(length(spnames), 5))
dimnames(VariancePartition)[[2]] = c("VarSp", "VarEnv", "TotVar", "NvarsSp", "NvarsEnv") 
dimnames(VariancePartition)[[1]] = spnames

# variable names 
envpars = c("DIAM", "DIAM2", "decay", "decay2")
envparscolsymb = c("beta", "beta2",  "betadecay", "betadecay2") 
envparsextsymb = c("gamma", "gamma2", "gammadecay", "gammadecay2") 

allspparscol = paste0("col.", spnames)
allspparsext = paste0("ext.", spnames)

## svss results
probinclusion = 0.2
whichvarsModelsSSVS = whichvarsModels
whichvarsModelsSSVS = whichvarsModelsSSVS[,!colnames(whichvarsModelsSSVS) %in% c("focalsp", "speciesCol", "speciesExt")]
whichvarsModelsSSVS[,] = 0

# make gradient of climate data for plotting
# load model data
load(paste(OUTDIR, modelnames[1], "/",whichvarsModels$focalsp[1], "_model.RData", sep="") )
nsamples = 50
rawEnvdata = scaledEnvdata = matrix(NA, nrow = nsamples, ncol = 4)
colnames(rawEnvdata) = colnames(scaledEnvdata) = envpars
for(parx in c("DIAM", "decay")){
  dat = c(alldatallsporig[[1]]$col[,parx],alldatallsporig[[1]]$ext[,parx]) 
  dat = seq(from=min(dat), to=max(dat), length.out=nsamples)
  rawEnvdata[,parx] = dat
  dattrans = (dat - mean(dat)) / sd(dat)
  scaledEnvdata[,parx] = dattrans
}
# for squared 
scaledEnvdata[,"DIAM2"] = scaledEnvdata[,"DIAM"]^2 
scaledEnvdata[,"decay2"] = scaledEnvdata[,"decay"]^2 

MxCol_Sub = MxExt_Sub = MxCol_Perc = MxExt_Perc = array(0, dim = list(length(spnames), length(spnames)), 
                                                        dimnames = list(spnames, spnames))

PPRallspecies = list()
mx=1 # big loop
for(mx in 1:nrow(whichvarsModels)){ #######################################################
  
  sp = whichvarsModels[mx,"focalsp"]
  sp
  splong = spnameslong[sp == spnames]
  
  print(paste0("species ", sp, ", ", mx," of ", nrow(whichvarsModels), " models", sep = " "))
  whichvarsModels[mx,]
  
  # load moalldatallsptransdel data
  load(paste(OUTDIR, modelnames[mx], "/",sp, "_model.RData", sep="") )
  ls()
  
  # save evaluation measures
  outt = MCMCsummary(post)
  outt

  # get parameter names for this species ############################
  # exclude itself
  otherspnames = spnames[spnames != sp]
  otherspnamescolsymb = rownames(outt)[grepl(rownames(outt), pattern = "^delta")]
  otherspnamesextsymb = rownames(outt)[grepl(rownames(outt), pattern = "^epsilon")]
  otherspnamesextsymb
  otherspnamescol = paste0("col.", otherspnames)
  otherspnamesext = paste0("ext.", otherspnames)
  
  onespvars = c(envpars, envpars, otherspnamescol, otherspnamesext)
  onespvarssymb = c(envparscolsymb, envparsextsymb, otherspnamescolsymb, otherspnamesextsymb)
  onespwhichvars = c(envparscolsymb, envparsextsymb, otherspnamescol, otherspnamesext)
  
  ###################
  # get svss results
  ssvsvarssymb = c(envparscolsymb, envparscolsymb, otherspnamescolsymb, otherspnamesextsymb)
  ssvsvarssymb = paste0("Ind", ssvsvarssymb)
  
  ssvsonesp = as.data.frame(whichvarsModelsSSVS[mx,onespwhichvars] )
  ssvsonesp = (outt[ssvsvarssymb,"mean"] - 1) / probinclusion
  names(ssvsonesp) = onespwhichvars
  whichvarsModelsSSVS
  # if BF above 3, the variable is important
  ssvsonesp = ifelse(ssvsonesp > 3, 1, 0)  
  
  # if quadratic is important, also make linear important and plot it
  if(ssvsonesp[["beta2"]] == 1)  ssvsonesp["beta"] = 1
  if(ssvsonesp[["betadecay2"]] == 1)  ssvsonesp["betadecay"] = 1
  if(ssvsonesp[["gamma2"]] == 1)  ssvsonesp["gamma"] = 1
  if(ssvsonesp[["gammadecay2"]] == 1)  ssvsonesp["gammadecay"] = 1
  
  whichvarsModelsSSVS[mx,onespwhichvars] = ssvsonesp
  print(whichvarsModelsSSVS[,1:3])
  
  # Variance partitioning, from Bob O'Hara  ###############  
  # logit(p[i]) <- beta0 + beta1*X1[i] + beta2*X2[i] 
  # then you can partition the variance on the link scale with 
  # TotVar <- beta1^2*Var(X1[]) + beta2^2*Var(X2[]) + beta1*beta2*Cov(X1[],X2[])
  # And if we ignore the final term, we can say that the proportion of variance due to X1 is = (beta1^2*Var(X1[]))/TotVar
  # We can extend this so the numerator is the sum of several covariate effects.
  
  # first which variables are important
  ssvsenvparscolsymb = envparscolsymb[whichvarsModelsSSVS[mx,envparscolsymb] == 1] 
  ssvsenvparsextsymb = envparsextsymb[whichvarsModelsSSVS[mx,envparsextsymb] == 1] 
  ssvsspparscol = otherspnamescol[whichvarsModelsSSVS[mx,otherspnamescol] == 1]
  ssvsspparsext = otherspnamesext[whichvarsModelsSSVS[mx,otherspnamesext] == 1]
  ssvsspparscolshort = substr(ssvsspparscol, 5, 12)
  ssvsspparsextshort = substr(ssvsspparsext, 5, 12)
  
  
  # get which parameter corresponds to which variable name
  onespvars 
  ssvsenvparscol = onespvars[onespvarssymb %in% ssvsenvparscolsymb]
  ssvsenvparsext = onespvars[onespvarssymb %in% ssvsenvparsextsymb]
  
  ssvsspparscolsymb = onespvarssymb[onespvars %in% ssvsspparscol]
  ssvsspparsextsymb = onespvarssymb[onespvars %in% ssvsspparsext]
  
  # get variance of environmentaland sp data
  envvariancecol = apply(as.matrix(alldatallsptrans[[sp]]$col[, ssvsenvparscol], ncol = length(ssvsenvparscol)), 2, var)
  envvarianceext = apply(as.matrix(alldatallsptrans[[sp]]$ext[, ssvsenvparsext], ncol = length(ssvsenvparsext)), 2, var)
  spvariancecol = apply(as.matrix(alldatallsptrans[[sp]]$col[, ssvsspparscolshort], ncol = length(ssvsspparscol)), 2, var)
  spvarianceext = apply(as.matrix(alldatallsptrans[[sp]]$ext[, ssvsspparsextshort], ncol = length(ssvsspparsext)), 2, var)
  
  # get parameter values
  envparacol = (outt[ssvsenvparscolsymb, "mean"])^2
  envparaext = (outt[ssvsenvparsextsymb, "mean"])^2
  spparacol = (outt[ssvsspparscolsymb, "mean"])^2
  spparaext = (outt[ssvsspparsextsymb, "mean"])^2
  
  TotVarEnv = sum(envparacol*envvariancecol) + sum(envparaext*envvarianceext)
  TotVarSp = sum(spparacol*spvariancecol) + sum(spparaext*spvarianceext)
  TotVarEnv
  TotVarSp
  TotVar = TotVarEnv + TotVarSp 
  VarianceExplainedSp = TotVarSp
  VarianceExplainedEnv = TotVarEnv
  
  VariancePartition[sp,] = c(VarianceExplainedSp, VarianceExplainedEnv, 
                             TotVar, length(c(ssvsspparscol, ssvsspparsext)), 
                             length(c(ssvsenvparscolsymb, ssvsenvparsextsymb)))
  
  
  # Partial regression plots #######################
  whichvarsModelsSSVS[mx,]
  ssvsenvparscolsymb
  ssvsenvparsextsymb
  processx = "col"
  processx = "ext"
  for(processx in processes){ ###################
    if(processx == "col"){
      envparmssymb = ssvsenvparscolsymb
      envparms = ssvsenvparscol
    } else {
      envparmssymb = ssvsenvparsextsymb
      envparms = ssvsenvparsext
    }
    envparmssymb
    envx = 1
    
    if(length(envparmssymb) == 0) next
    for(envx in 1:length(envparmssymb)){
      
      inttcept =  outt[paste0("alpha", processx), "mean"]
      param = outt[envparmssymb[envx], "mean"]
      param
     if(!grepl(envparmssymb[envx], pattern = "2")){ # if linear parameter
        
        valstoplott =  plogis(inttcept + (param * scaledEnvdata[,envparms[envx]]))
        
        PPRallspecies[[sp]][[processx]][[envparmssymb[envx]]] = valstoplott
     } else { # if quadratic, add linear response too and overwrite linear response
        
        paramlinear = outt[envparmssymb[envx-1], "mean"]
        paramlinear
        valstoplott = plogis(inttcept + 
                               (paramlinear * scaledEnvdata[,envparms[envx-1]]) +
                               (param * scaledEnvdata[,envparms[envx]])
        )
       # save under the linear name, replace linear response
        PPRallspecies[[sp]][[processx]][[envparmssymb[envx-1]]]  = valstoplott
        
      } # end if linear or quadratic
      
      # plot #################
      envparms[envx]
      envparms[envx]
      if(grepl(pattern = "2", envparms[envx])) next
       if(envparms[envx] == "decay") envnameorigplot = "Log decay (% soft wood)"
      if(envparms[envx] == "DIAM") envnameorigplot = "Log diameter (cm)"
      
      xlims = range(rawEnvdata[, envparms[envx]])
      xlims
      if(envparms[envx] == "DIAM"){
        xlims = c(min(rawEnvdata[, envparms[envx]]), 60)
      }
      if(processx == "col"){
        ylims = c(0,.4)
      } else {
        ylims =  c(0,1)
      }
      
      # plot one species #########
      plot(rawEnvdata[, envparms[envx]], rep(0,nrow(rawEnvdata)),
           type = "l", xlab= paste0(envnameorigplot), ylab = paste0("Probability of ", processeslong[processx]), #  " (",envnameorigsymb,")"
           xlim = xlims, ylim=ylims, col="white", main = splong) #,  " vs. ","P",processx))
      valstoplot = PPRallspecies[[sp]][[processx]][[envparmssymb[envx]]]
      if(processx == "ext"){
        # if extinction, invert so pro of ext, not prob of persistence
        valstoplot[] = 1 - valstoplot
      }
      lines(rawEnvdata[,envparms[envx]], valstoplot,
            col=cols[sp], lty = ltys[mx], lwd = 1.5) 
      # no confidence intervals
      
    } # end loop env variables
    
  } # end processes partial regression plots
  
  # species interactions #############################
  processx = "col"
  processx = "ext"
  for(processx in processes){ 
    if(processx == "col"){
      ssvsspparssymb = ssvsspparscolsymb
      ssvsspparsshort = ssvsspparscolshort 
    } else {
      ssvsspparssymb = ssvsspparsextsymb
      ssvsspparsshort = ssvsspparsextshort 
    }
    inttcept = outt[paste0("alpha", processx), "mean"]
    inttcept
    predAlone = plogis(inttcept)
    # loop over if each other sp present or absent
    spxx = 1
    spxx
    sp
    for(spxx in 1:length(ssvsspparsshort)){ # long loop
      print(processx)
      print(spxx)
      if(length(ssvsspparsshort) == 0) next
      spxxsymb = ssvsspparssymb[spxx]
      spxxsymb
      spxxname = ssvsspparsshort[spxx] #pars[[processx]]$otherspnamesonly[ which(spxx == pars[[processx]]$othersp)]
      spxxname
      param = outt[spxxsymb, "mean"]
      param
      predTogether = plogis(inttcept + param) # assume other covars at zero and no other sp
      
      # subtraction method
      estInteractionSub = predTogether - predAlone
      # percentage method
      estInteractionPerc = ((predTogether-predAlone) / predAlone) *100
      
      estInteractionSub = round(estInteractionSub, 3)
      estInteractionPerc = round(estInteractionPerc, 3)
      
      if(processx == "col"){
        MxCol_Sub[sp,spxxname]
        MxCol_Sub[sp,spxxname] = estInteractionSub  
        MxCol_Perc[sp,spxxname] = estInteractionPerc
      } else { # extinction
        MxExt_Sub[sp,spxxname] = -estInteractionSub
        MxExt_Perc[sp,spxxname] = -estInteractionPerc
      } 
    }
  }
  
} # end long sp loop ##########################

# plot interaction strengths #############################
dimnames(MxCol_Sub) = dimnames(MxExt_Sub) = dimnames(MxCol_Perc) = dimnames(MxExt_Perc) =
  list(paste0(":italic(", spnamesmediumlong, ")"), paste0(":italic(", spnamesmediumlong, ")"))
diag(MxCol_Sub) = diag(MxExt_Sub) = diag(MxCol_Perc) = diag(MxExt_Perc) = NA

colrampp = colorRampPalette(c("darkred","white","darkgreen"))(200)


jpeg(paste0(OUTDIR, "/_InteractionMxCol_presentation.jpg"), width = 17, height = 16, units = "cm", res = 400)
par(mfrow = c(1,1), mar = c(5,7,6,4), cex.lab = 1) #  c(bottom, left, top, right)
range(MxCol_Sub, na.rm=T)
corrplot(MxCol_Sub, main = "Colonization interactions", mar=c(0,0,1,0), method = "color", is.corr = F, col.lim = c(-.1,.25),
         col=colrampp, na.label = "square", na.label.col = "grey", tl.col = 1, tl.cex = .8, cl.cex = 1., 
         cl.align.text = "l", cl.ratio = .25)
mtext("Predecessor (species variable)", side=3, line=3, cex = 1.2) #Primary colonizers"
mtext("Succcessor (focal species)       ", side=2, line=6, cex = 1.2) # Later-arriving species
mtext("Interaction strength                      ", side=4, line=3, cex = 1.2) # Later-arriving species
dev.off()

jpeg(paste0(OUTDIR, "/_InteractionMxExt_presentation.jpg"), width = 17, height = 16, units = "cm", res = 400)
par(mfrow = c(1,1), mar = c(5,7,6,4), cex.lab = 1) #  c(bottom, left, top, right)
range(MxExt_Sub, na.rm=T)
corrplot(MxExt_Sub, main = "Colonization interactions", mar=c(0,0,1,0), method = "color", is.corr = F, col.lim = c(-.25,.25),
         col=colrampp, na.label = "square", na.label.col = "grey", tl.col = 1, tl.cex = .8, cl.cex = 1., 
         cl.align.text = "l", cl.ratio = .25)
mtext("Predecessor (species variable)", side=3, line=3, cex = 1.2) #Primary colonizers"
mtext("Succcessor (focal species)       ", side=2, line=6, cex = 1.2) # Later-arriving species
mtext("Interaction strength                      ", side=4, line=3, cex = 1.2) # Later-arriving species
dev.off()

# save for comparing results with informative and non-informative priors
# in "plots_compare_priors.R"
save(MxCol_Sub, MxExt_Sub, MxCol_Perc, MxExt_Perc, 
     file = paste0(OUTDIR, "Matrices_results.Rdata"))


# plot partial regression plots all species ############################
envparms = c("DIAM", "decay")
for(processx in processes){ 
  if(processx == "col"){
    envparmssymb = c("beta", "betadecay")
  } else {
    envparmssymb = c("gamma", "gammadecay")
  }
  
  envparmssymb
  envx = 1
  jpeg(filename = paste(OUTDIR,"_PPR_allsp_",processx,"_confint.jpeg", sep=""), 
       width = 1500, height = 600, quality = 100, res = 250)
  par(mfrow=c(1,3), mar = c(5,4,2,1)) # c(bottom, left, top, right)
  for(envx in 1:length(envparmssymb)){
    
    xlims = range(rawEnvdata[, envparms[envx]]) 
    xlims
    if(envparms[envx] == "DIAM"){
      xlims = c(min(rawEnvdata[, envparms[envx]]), 60)
    }
    if(processx == "col"){
      ylims = c(0,.4)
    } else {
      ylims =  c(0,1)
    }
    if(envparms[envx] == "decay") envnameorigplot = "Log decay (% soft wood)"
    if(envparms[envx] == "DIAM") envnameorigplot = "Log diameter (cm)"
    
    plot(rawEnvdata[, envparms[envx]], rep(0,nrow(rawEnvdata)),
         type = "l", xlab= paste0(envnameorigplot), ylab = paste0("Probability of ", processeslong[processx]), #  " (",envnameorigsymb,")"
         xlim = xlims, ylim=ylims, col="white", main = "") 
    
    for(spx in 1:length(spnames)){ #######################################################
      sp = spnames[spx]
      valstoplot = PPRallspecies[[sp]][[processx]][[envparmssymb[envx]]]
      print(valstoplot)
      if(length(valstoplot) == 0) next
      if(processx == "ext"){
        # if extinction, invert so pro of ext, not prob of persistence
        valstoplot[] = 1 - valstoplot
      }
      lines(rawEnvdata[,envparms[envx]], valstoplot,
            col=cols[sp], lty = ltys[mx], lwd = 1.5) 
    }
  }
  par(mar = c(4,0,0,0))
  plot(1, col = "white", ylab = "", xlab = "", axes = F)
  legend("bottomleft", col = cols, text.font = 3, legend = spnamesmediumlongspace, lty = ltys, title = "Species", bty = "n", lwd = 1.5) #rev(1+(1:spx/5) )
  dev.off()
  
  
}

whichvarsModelsSSVS


png(paste0(OUTDIR, "_VariancePartition.png"), width = 800, height = 400, units = "px")
par(mfrow = c(1,2), mar = c(10,5,5,1))
barplot(t(VariancePartition[,1:2]/VariancePartition[,"TotVar"]), col = c("blue", "sandybrown"), ylim = c(0,1.3), xlim = c(-1,length(spnames)),
        axisnames = F, axes = F, las = 2, space = 0, main = "", ylab = "Variance explained \n(proportion of total variance)")
#text(x = (1:12)-.5, y = 1.15, labels = paste("Sp:\n", VariancePartition[,"NvarsSp"], "\n", "Env:\n", VariancePartition[,"NvarsEnv"]), col = 1, cex = 1.2)
text(x = (1:12)-.5, y = 1.13, labels = paste0(VariancePartition[,"NvarsEnv"],"\n", VariancePartition[,"NvarsSp"]), col = 1, cex = 1.)
text(x = -.5, y = 1.13, labels = paste("N env:\nN sp:"), col = 1, cex = 1.)
#text (x = 9.5, y = 1.25,labels = "later \nsuccessional species")
text (x = 2.5, y = 1.26,labels = "Successional stage", cex = .95)
arrows(x0 = 5, y0 = 1.26, x1 = 10, y1 = 1.26, length = 0.1, angle = 30, code = 2) 
axis(1, labels = spnamesmediumlongspace, font = 3, at = 0.5:(length(spnames)-0.5), las = 2)
#axis(1, labels = paste0(":italic(", spnamesmediumlong, ")"), at = 0.5:(length(spnames)-0.5), las = 2)
axis(2, at = seq(0,1, 0.25))
par(mar = c(1,1,1,1))
plot(1, main = "", axes = F, col = "white", ylab = "", xlab = "")
legend("left", legend = c("Environmental variables", "Species variables"), fill = c("sandybrown", "blue"))
dev.off()


