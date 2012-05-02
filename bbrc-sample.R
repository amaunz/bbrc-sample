# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# Author: Andreas Maunz
# Year: 2012


# BBRC bootstrapping
# Sampling BBRC descriptors with non-parametric, stratified bootstrapping
# Each bootstrap sample is processed via a webservice (see SERVER) and parallel processing
# NOTE: the server is loaded with as many jobs in parallel as there are CPU cores on the client
# Probabilistic chi-square test is performed on each sampled pattern using a poisson MLE estimate

# packages installed
if (sum(installed.packages()[,1]=="RCurl")==0) install.packages('RCurl')
if (sum(installed.packages()[,1]=="doMC")==0) install.packages('doMC')

# load packages and libs
source("bbrc-sample-psqerr.R")
suppressPackageStartupMessages(library('RCurl'))
suppressPackageStartupMessages(library('doMC'))
registerDoMC()

# set global variables here
SERVER <- "http://toxcreate3.in-silico.ch:8082"
HAMSTER <- "http://toxcreate3.in-silico.ch:8082/dataset/5092"
BBRC <- paste(SERVER,"algorithm","fminer","bbrc",sep='/')
BOOTS <- 100
MIN_FREQUENCY_PER_SAMPLE <- 8
MIN_SAMPLING_SUPPORT <- 10


# # # Main function
bootBbrc = function(datasetUri, # dataset to process (URI)
                    numboots=BOOTS, # nr of BS samples
                    min.frequency.per.sample=MIN_FREQUENCY_PER_SAMPLE, # min freq inside each sample
                    min.sampling.support=MIN_SAMPLING_SUPPORT, # min nr of BS samples that have this pattern
                    del=NULL) {

  set.seed(1)

  # load dataset
  ds <- getDataset(datasetUri)
  ds.endpoint.type = class(ds[,2])
  (ds.endpoint.type == "numeric" || ds.endpoint.type == "character") || return("Wrong endpoint type")
  ds.factors <- factor(ds[,2])
  ds.levels <<- levels(ds.factors)
  ds.table <- table(ds.factors)

  n <- dim(ds)[1]
  cat(paste("Full set size:", n, "\n"))

  # main loop
  if(is.null(del)) {
    bb <- foreach(j=1:numboots, .combine=mergeLists) %dopar% {
    #for (j in 1:numboots) {
      
      idx <- c()
      for (fac in ds.levels) {
        if (ds.endpoint.type == "numeric") fac = as.numeric(fac)
        draw <- drawSample(which(ds[,2] == fac))
        idx <- c(idx, draw)
      }
      sample <- ds[idx,]
  
      task <- postDataset(sample, tempFilePrefix=paste("boot_bbrc_sample_",j,"_", sep=""))
      sampleUri <- getResult(task)
  
      task <- postRequest(BBRC, list( dataset_uri=sampleUri, 
                                      min_frequency=as.character(min.frequency.per.sample)))
      sampleFeaturesUri <- getResult(task)
  
      sampleFeatures <- getDataset(sampleFeaturesUri)
      sampleFeatures["SMILES"] <- NULL
      nr.features <- dim(sampleFeatures)[2]
      class.support <- apply(sampleFeatures,2,function(x) supportPerFactor(x,sample[,2],ds.levels))
  
      deleteRequest(sampleFeaturesUri)
      deleteRequest(sampleUri)
     
      as.list(as.data.frame(class.support))
    }
    bb <- data.frame(bb,check.names=F)
    bb$levels <- rep(ds.levels, dim(bb)[1]/length(ds.levels))
  }
  else
    bb <- del
  
  # Filter patterns with enough sampling support
  cat("\nFiltering\n")
  enough.data=rep(F,length(names(bb)))
  for (l in ds.levels) { 
    mask <- apply(bb[bb$levels==l,,drop=F], 2, function(x) { sum(complete.cases(x)) > min.sampling.support } )
    names(mask) <- NULL
    enough.data <- enough.data | mask
  }
  cat(paste("Stripped",dim(bb)[2]-sum(enough.data),"patterns, kept",sum(enough.data),"\n"))
  bb <- bb[,enough.data]
  

  # Get chi-square
  cat("\nChisq\n")
  for (p in names(bb)[names(bb) != "levels"]) {
    cat("\n")
    
    # sp
    sp <- rep(0,length(bb$levels) / length(ds.levels))
    for (l in ds.levels) sp <- sp + bb[bb$levels==l,p]
    sp <- sp[complete.cases(sp)] 
    
    # dsp
    dsp <- prop.table(table(sp)) # Categorical distribution for p
    
    chisq <- 0.0
    for (l in ds.levels) {
      
      spl <- bb[bb$levels==l,p]
      spl <- spl[complete.cases(spl)]
      
      # lambda.spl: Distribution parameters for p on level l
      lambda.spl <- list()
      for (idx in names(dsp)) {
        lambda.spl[[idx]] <- mean(spl[sp==idx])
        p.weight <- as.numeric(idx)/n # weight for p
        spx <- ds.table[l] * p.weight # eXpected support for p
        chisq <- chisq + dsp[idx] * rectintegrate(f=wSquaredErr,from=0,to=10000,lambda=lambda.spl[[idx]],spx=spx)
      }
      
    }
    cat(paste(p,":",chisq,"\n"))
  }
  cat("\nDone\n")
  bb

}

# merges 2nd list to 1st and returns result
# suitable for .combine
mergeLists = function(x, xn, levels=length(ds.levels)) {

  x<-data.frame(x,check.names=F)
  xn<-data.frame(xn,check.names=F)

  xus <- names(x)[names(x) %in% names(xn)]
  xnn <- names(xn)[!names(xn) %in% names(x)]

  if (length(xus)>0) {
    x.bel <- data.frame(matrix(NA,levels,dim(x)[2]))
    names(x.bel) <- names(x)
    x.bel[1:levels,xus] <- xn[1:levels,xus]
    x <- rbind(x,x.bel)
  }
  
  if (length(xnn)>0) {
    x.rig <- data.frame(matrix(NA,dim(x)[1],length(xnn)))
    names(x.rig) <- xnn
    x.rig[c(dim(x)[1]-levels+1,dim(x)[1]), xnn] <- xn[1:levels, xnn]
    x <- cbind(x,x.rig)
  }

  as.list(x)

}

# draw a bootstrap sample
drawSample = function(x, size=length(x)) {
  sample(x, size=size, replace=T)
}

# get dataset from task
getResult = function(uri) {
  while(1) { 
    Sys.sleep(0.5)
    uri <- getRequest(uri,"text/uri-list")
    if (length(grep("dataset",uri))>0) break
  }
  uri
}

# get per-factor support in x from y, with levels from z
supportPerFactor = function(x, y, levels) {
  pattern <- factor(y, levels=levels)
  table(pattern[x>0]) # table sums across the factors
}


# # # REST library

# GET
getRequest = function(uri, accept) {
  getURL(uri, httpheader=c(Accept=accept))
}

# POST (assumes a task is returned as uri-list)
postRequest = function(uri, params) {
  postForm(uri, .params=params)
}

# DELETE a dataset by URI
deleteRequest = function(uri) {
  httpDELETE(uri)
}

# GET a dataset by URI
getDataset = function(uri) {
  read.csv(paste(uri, ".csv", sep=""), header=T, check.names=F, stringsAsFactors=F)
}

# POST a dataset (CSV)
postDataset = function(x,tempFilePrefix="R_") {
  tf <- tempfile(pattern=paste(tempFilePrefix, Sys.getpid(), sep=""))
  cat(paste(gsub(".*boot", "boot", gsub(gsub(".*_",  "_", tf),"", tf)),"\n"))
  tryCatch({
    write.csv(x=x, file=tf, row.names=F, quote=F, na='')
    task <- postForm(paste(SERVER,"dataset",sep='/'), file=fileUpload(filename=tf, contentType="text/csv"), .checkParams=F )    
  }, finally = {
    unlink(tf)
  })
  task
}
