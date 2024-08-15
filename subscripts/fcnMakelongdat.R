
# get unscaled data for all sp
makelongdat <- function(sparray, spnames, othersp, colparms, extparms, 
                        envdata, decaydat, decaydat2){ # for one sp data...
  
  processes = processes = c("col", "ext")
  alldatallsp = list()
  spx=1
  for(spx in 1:length(spnames)){
    sp = spnames[spx]
    sp
    
    alldat = list()
    
    # make species-specific data
    spdat = sparray[, sp, , drop = T]
    spdat = matrix(spdat, ncol = dim(sparray)[[3]], dimnames = list(NULL, dimnames(sparray)[[3]]))
    
    otherspdat = sparray[, othersp, , drop = T]# array 
    otherspdat
    # decaydat not needed to subset for sp, same for allsp
    processx = "ext"
    for(processx in processes){
      alldatone = NA
      yx=1
      for(yx in 1:(ncol(spdat)-1)){
        oneyr = spdat[,yx]
        #print(yx)
        #print(sum(oneyr))
        #if(is.na(sum(oneyr))) {next}
        if(processx == "col"){
          whichlines = oneyr == 0 # when colonization
          envnamessp = colparms
          envv = envdata[,envnamessp]
          
        } else {
          whichlines = oneyr == 1 # when extinction    
          envnamessp = extparms
          envv = envdata[,envnamessp]
          
        }
        whichlines[is.na(whichlines)] = F
        table(names(whichlines))
        table(whichlines)
        dim(envv)
        if(length(envnamessp) == 1){
          envv
          #envvsp = envv[whichlines,]
          envvsp = as.matrix(envv, ncol = length(envnamessp))
          colnames(envvsp) = envnamessp
        }
        if(length(envnamessp) > 1){
          envv
          #envvsp = envv[whichlines,]
          envvsp = as.matrix(envv, ncol = length(envnamessp))
          colnames(envvsp) = envnamessp
        }
        envvsp # should always be matrix 
        envvsp = data.frame(envvsp)
        # subset
        envvsp2 = envvsp[whichlines,, drop = F]
        if(length(envvsp)>0) { 
          colnames(envvsp) = envnamessp 
        }
        otherspdatt = matrix(otherspdat[whichlines,,yx], ncol = length(othersp))
        colnames(otherspdatt) = othersp
        decaydat[whichlines,]
        dim(decaydat)
        onedat = cbind(t1 = spdat[whichlines,yx], t2 = spdat[whichlines,yx+1],
                       envvsp2,
                       decay = decaydat[whichlines,yx],
                       decay2 = decaydat2[whichlines,yx],
          #             disp = dispdat[whichlines,yx],
                       otherspdatt
        )
        head(onedat)
        if(nrow(onedat)>0){
          alldatone = rbind(alldatone, onedat)
        }
      }
      if(length(dim(alldatone[1])) == 2){
        if(is.na(alldatone[1,1])){
          alldatone = alldatone[-1,] # remove first row
        }
      }
      head(alldatone)
      alldat[[processx]] = alldatone
    }
    alldatallsp[[sp]] = alldat
  } #end sp
  #str(alldat)
  return(alldatallsp)
}
