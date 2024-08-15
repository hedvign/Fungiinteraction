
# compares results with and without informative priors
rm(list = ls())
setwd("~/Documents/Polypores/MODELfungiinteractions")

library(corrplot)


OUTDIR = "livefruitbodies_InfPriors/" # Informative priors, OUTDIR
OUTDIR2 = "livefruitbodies_noInfPriors/" # Uninformative priors, OUTDIR2


colrampp = colorRampPalette(c("darkred","white","darkgreen"))(200)

# Plot interaction webs #####################################################
# load ############################################
load(paste0(OUTDIR, "Matrices_results.Rdata"))
ls()
MxCol_Sub  # subtraction method
MxCol_Perc # percentage method
# interaction strengths
MxCol_SubInfpriors = MxCol_Sub
MxExt_SubInfpriors = MxExt_Sub
MxCol_PercInfpriors = MxCol_Perc
MxExt_PercInfpriors = MxExt_Perc

load(paste0(OUTDIR2, "Matrices_results.Rdata"))
MxCol_SubnoInfpriors = MxCol_Sub
MxExt_SubnoInfpriors = MxExt_Sub
MxCol_PercnoInfpriors = MxCol_Perc
MxExt_PercnoInfpriors = MxExt_Perc


fig_label <- function(text, region="figure", pos="topleft", cex=NULL, ...) {
  
  region <- match.arg(region, c("figure", "plot", "device"))
  pos <- match.arg(pos, c("topleft", "top", "topright", 
                          "left", "center", "right", 
                          "bottomleft", "bottom", "bottomright"))
  
  if(region %in% c("figure", "device")) {
    ds <- dev.size("in")
    # xy coordinates of device corners in user coordinates
    x <- grconvertX(c(0, ds[1]), from="in", to="user")
    y <- grconvertY(c(0, ds[2]), from="in", to="user")
    
    # fragment of the device we use to plot
    if(region == "figure") {
      # account for the fragment of the device that 
      # the figure is using
      fig <- par("fig")
      dx <- (x[2] - x[1])
      dy <- (y[2] - y[1])
      x <- x[1] + dx * fig[1:2]
      y <- y[1] + dy * fig[3:4]
    } 
  }
  
  # much simpler if in plotting region
  if(region == "plot") {
    u <- par("usr")
    x <- u[1:2]
    y <- u[3:4]
  }
  
  sw <- strwidth(text, cex=cex) * 180/100
  sh <- strheight(text, cex=cex) * 60/100
  
  x1 <- switch(pos,
               topleft     =x[1] + sw , 
               left        =x[1] + sw,
               bottomleft  =x[1] + sw,
               top         =(x[1] + x[2])/2,
               center      =(x[1] + x[2])/2,
               bottom      =(x[1] + x[2])/2,
               topright    =x[2] - sw,
               right       =x[2] - sw,
               bottomright =x[2] - sw)
  
  y1 <- switch(pos,
               topleft     =y[2] - sh,
               top         =y[2] - sh,
               topright    =y[2] - sh,
               left        =(y[1] + y[2])/2,
               center      =(y[1] + y[2])/2,
               right       =(y[1] + y[2])/2,
               bottomleft  =y[1] + sh,
               bottom      =y[1] + sh,
               bottomright =y[1] + sh)
  
  old.par <- par(xpd=NA)
  on.exit(par(old.par))
  
  text(x1, y1, text, cex=cex, ...)
  return(invisible(c(x,y)))
} 


# interaction strenghts #####################
cexxmain = 1.2
cexspeciestypes = .95
cexxlabel = 1.5
range(MxCol_SubInfpriors, na.rm = T)

# percentage and subtraction method
jpeg(filename = paste(OUTDIR, "_allcompareinteractions_Col_percSub.jpeg", sep=""), 
     width = 25, height = 20, units = "cm", res = 400)
par( mfrow = c(2,2), mar = c(4,4,4,4), cex.lab = 1)   # c(bottom, left, top, right)
corrplot(MxCol_SubInfpriors, main = paste0("Colonization \nInformative priors (n = ",sum(MxCol_SubInfpriors != 0, na.rm=T), ")"),
         mar=c(0,0,2,0), method = "color",is.corr = F, col.lim = c(-.05,.3), cex.main = cexxmain,
         col=colrampp, na.label = "square", na.label.col = "grey89", tl.col = 1, tl.cex =  .8, cl.cex = 1.,
         cl.align.text = "l", cl.ratio = .25)
fig_label("(a)", cex=cexxlabel, pos="topleft", region="plot")
mtext("Succcessor (focal species)        ", side=2, line=2, cex = cexspeciestypes) # Later-arriving species
mtext("Predecessor (species variable)", side=3, line=1, cex = cexspeciestypes) 
range(MxCol_SubnoInfpriors, na.rm=T)
corrplot(MxCol_SubnoInfpriors, main = paste0(" \n Non-informative priors (n = ",sum(MxCol_SubnoInfpriors != 0, na.rm=T), ")"),
         mar=c(0,0,2,0), method = "color",is.corr = F, col.lim = c(-.05,.3), cex.main = cexxmain,
         col=colrampp, na.label = "square", na.label.col = "grey89", tl.col = 1, tl.cex =  .8, cl.cex = 1.,
         cl.align.text = "l", cl.ratio = .25)
fig_label("(b)", cex=cexxlabel, pos="topleft", region="plot")
mtext("Interaction strength              ", side=4, line=2, cex = 1.2) # Later-arriving species
# percentage method
range(MxCol_PercInfpriors, na.rm = T)
range(MxCol_PercnoInfpriors, na.rm=T)
corrplot(MxCol_PercInfpriors, main = "", #paste0("Colonization \nInformative priors"),
         mar=c(0,0,2,0), method = "color",is.corr = F, col.lim = c(-100,600),
         col=colrampp, na.label = "square", na.label.col = "grey89", tl.col = 1, tl.cex =  .8, cl.cex = 1.,
         cl.align.text = "l", cl.ratio = .25)
fig_label("(c)", cex=cexxlabel, pos="topleft", region="plot")
mtext("Succcessor (focal species)        ", side=2, line=2, cex = cexspeciestypes) # Later-arriving species
mtext("Predecessor (species variable)", side=3, line=1, cex = cexspeciestypes) 
corrplot(MxCol_PercnoInfpriors, main = "",# paste0(" \n Non-informative priors"),
         mar=c(0,0,2,0), method = "color",is.corr = F, col.lim = c(-100,600),
         col=colrampp, na.label = "square", na.label.col = "grey89", tl.col = 1, tl.cex =  .8, cl.cex = 1.,
         cl.align.text = "l", cl.ratio = .25)
mtext("% colonization change               ", side=4, line=2, cex = 1.2) # Later-arriving species
fig_label("(d)", cex=cexxlabel, pos="topleft", region="plot")
dev.off()


# extinction
jpeg(filename = paste(OUTDIR, "_allcompareinteractions_Ext_percSub.jpeg", sep=""),
     width = 25, height = 20, units = "cm", res = 400)
par( mfrow = c(2,2), mar = c(4,4,4,4), cex.lab = 1)   # c(bottom, left, top, right)
range(MxExt_SubInfpriors, na.rm=T)
range(MxExt_SubnoInfpriors, na.rm=T)
# subtraction method
corrplot(MxExt_SubInfpriors, main = paste0("Extinction \nInformative priors (n = ",sum(MxExt_SubInfpriors != 0, na.rm=T), ")"),
         mar=c(0,0,2,0), method = "color",is.corr = F, col.lim = c(-0.3,.5), cex.main = cexxmain,
         col=colrampp, na.label = "square", na.label.col = "grey89", tl.col = 1, tl.cex =  .8, cl.cex = 1.,
         cl.align.text = "l", cl.ratio = .25)
mtext("Succcessor (focal species)        ", side=2, line=2, cex = cexspeciestypes) # Later-arriving species
mtext("Predecessor (species variable)", side=3, line=1, cex = cexspeciestypes) 
fig_label("(a)", cex=cexxlabel, pos="topleft", region="plot")
corrplot(MxExt_SubnoInfpriors, main = paste0(" \n Non-informative priors (n = ",sum(MxExt_SubnoInfpriors != 0, na.rm=T), ")"),
         mar=c(0,0,2,0), method = "color",is.corr = F, col.lim = c(-0.3,.5), cex.main = cexxmain,
         col=colrampp, na.label = "square", na.label.col = "grey89", tl.col = 1, tl.cex =  .8, cl.cex = 1.,
         cl.align.text = "l", cl.ratio = .25)
mtext("Interaction strength              ", side=4, line=2, cex = 1.2) # Later-arriving species

fig_label("(b)", cex=cexxlabel, pos="topleft", region="plot")
# percentage method
range(MxExt_PercInfpriors, na.rm=T)
range(MxExt_PercnoInfpriors, na.rm=T)
corrplot(MxExt_PercInfpriors, main = "",
         mar=c(0,0,2,0), method = "color",is.corr = F, col.lim = c(-100,100),
         col=colrampp, na.label = "square", na.label.col = "grey89", tl.col = 1, tl.cex =  .8, cl.cex = 1.,
         cl.align.text = "l", cl.ratio = .25)
mtext("Succcessor (focal species)        ", side=2, line=2, cex = cexspeciestypes) # Later-arriving species
mtext("Predecessor (species variable)", side=3, line=1, cex = cexspeciestypes) 
fig_label("(c)", cex=cexxlabel, pos="topleft", region="plot")
corrplot(MxExt_PercnoInfpriors, main = "",#paste0(" \n Non-informative priors"),
         mar=c(0,0,2,0), method = "color",is.corr = F,col.lim = c(-100,100),
         col=colrampp, na.label = "square", na.label.col = "grey89", tl.col = 1, tl.cex =  .8, cl.cex = 1.,
         cl.align.text = "l", cl.ratio = .25)
mtext("% extinction change              ", side=4, line=2, cex = 1.2) # Later-arriving species
fig_label("(d)", cex=cexxlabel, pos="topleft", region="plot")
dev.off()

