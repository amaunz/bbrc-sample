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

# Class-correlated subgraph descriptors are calculated repeatedly on bootstrap samples of a database of molecular graphs, where each graph is associated with one from a finite set of classes.
# Each sample is processed by the graph mining webservice architecture. In addition to using the efficient BBRC algorithm for mining the graphs, the server is loaded with several jobs in parallel, which further increases performance.
# A probabilistic chi-square test is performed on the most frequently sampled descriptors (patterns) using a poisson MLE estimate.
# Patterns surviving the test are mapped back on the database molecules, using a parallelized molecular fragment matching method.


# packages installed
if (sum(installed.packages()[,1]=="RCurl")==0) install.packages('RCurl')
if (sum(installed.packages()[,1]=="doMC")==0) install.packages('doMC')

# load packages and libs
if (file.exists("bbrc-sample-psqerr.R")) { 
  source("bbrc-sample-psqerr.R") 
} else { 
  source("bbrc-sample/bbrc-sample-psqerr.R") 
}
suppressPackageStartupMessages(library('RCurl'))
suppressPackageStartupMessages(library('doMC'))
registerDoMC(4)

# set global variables here
SERVER <- "http://toxcreate3.in-silico.ch:8082"
HAMSTER <- "http://toxcreate3.in-silico.ch:8082/dataset/5092"
BBRC <- paste(SERVER,"algorithm","fminer","bbrc",sep='/')
BOOTS <- 100
MIN_FREQUENCY_PER_SAMPLE <- 8
MIN_SAMPLING_SUPPORT <- 10


# # # Main function
bootBbrc = function(dataset.uri, # dataset to process (URI)
                    prediction.feature.uri=NULL,
                    num.boots=BOOTS, # nr of BS samples
                    min.frequency.per.sample=MIN_FREQUENCY_PER_SAMPLE, # min freq inside each sample
                    min.sampling.support=MIN_SAMPLING_SUPPORT, # min nr of BS samples that have this pattern
                    del=NULL,
                    bbrc.service=BBRC,
                    dataset.service=paste(SERVER,"dataset",sep='/'),
                    do.ot.log=F,
                    random.seed=1) {

  do.ot.log <<- do.ot.log # should be global
  merge.time.start <<- NULL
  merge.time.end <<- NULL

  # load dataset
  ds <- getDataset(dataset.uri)
  names(ds) = curlUnescape(names(ds)) # names unencoded in R

  if (is.null(prediction.feature.uri)) {
    if (length(names(ds))>2) { 
      return("Too many columns")
    }
    else {
      ds.endpoint <- names(ds)[2]
      if (do.ot.log) otLog(paste("Endpoint (guessed from dataset):",ds.endpoint))
    }  
  } 
  else {
    ds.endpoint <- curlUnescape(gsub(".*/","",prediction.feature.uri)) # select endpoint
    if (do.ot.log) otLog(paste("Endpoint (obtained as parameter):",ds.endpoint))
  }

  ds.endpoint.type <- class(ds[,ds.endpoint])
  if (ds.endpoint.type != "numeric" && ds.endpoint.type != "character") return("Wrong endpoint type")
  ds.factors <- factor(ds[,ds.endpoint])
  ds.levels <<- levels(ds.factors)
  ds.table <- table(ds.factors)

  ds.n <- dim(ds)[1]

  # main loop
  if(is.null(del)) {
    bb <<- foreach(j=1:num.boots, .combine=mergeLists) %dopar% {
    #for (j in 1:num.boots) {
      set.seed(j+random.seed) 
      idx <- c()
      for (fac in ds.levels) {
        if (ds.endpoint.type == "numeric") fac = as.numeric(fac)
        draw <- drawSample(which(ds[,ds.endpoint] == fac))
        idx <- c(idx, draw)
      }
      ds.sample <- ds[idx,,drop=F]
      ds.oob <- ds[! c(1:ds.n) %in% unique(idx),,drop=F]
      task <- postDataset(ds.sample, dataset.service, tempFilePrefix=paste("boot_bbrc_sample_",j,"_", sep=""))
      sampleUri <- getResult(task)
      task <- postDataset(ds.oob, dataset.service, tempFilePrefix=paste("boot_bbrc_oob",j,"_", sep=""))
      oobUri <- getResult(task)

      prediction.feature.uri <- paste(sampleUri, "feature", curlEscape(ds.endpoint), sep="/")
      bbrc.params <- list( dataset_uri=sampleUri, prediction_feature=prediction.feature.uri, min_frequency=as.character(min.frequency.per.sample))
      task <- postRequest(bbrc.service, bbrc.params)
      sampleFeaturesUri <- getResult(task)
      bbrc.params <- list( dataset_uri=oobUri, feature_dataset_uri=sampleFeaturesUri)
      task <- postRequest(paste(bbrc.service,"match",sep="/"),bbrc.params)
      oobFeaturesUri <- getResult(task)
  
      sampleFeatures <- getDataset(oobFeaturesUri)
      sampleFeatures <- as.data.frame(t(apply(ds.oob, 1, sampleFeatureRow, fds=sampleFeatures)), stringsAsFactors=F)
      sampleFeatures["SMILES"] <- NULL
      apply( sampleFeatures[,names(sampleFeatures) != "SMILES"], 2, function(x) as.numeric(x) )
      class.support <- apply(sampleFeatures,2,function(x) supportPerFactor(x,ds.oob[,ds.endpoint],ds.levels))
  
      deleteRequest(sampleFeaturesUri)
      deleteRequest(sampleUri)
      deleteRequest(oobUri)
     
      as.list(as.data.frame(class.support))
    }

    merge.time.end <<- Sys.time()
    merge.time <<- as.numeric(merge.time.end - merge.time.start, units="secs")
    data.frame(bb,check.names=F)
  }
  else
    bb <<- del

  levels <- rep(ds.levels, dim(bb)[1]/length(ds.levels))
  # Filter patterns with enough sampling support
  if (do.ot.log) otLog("Filtering")
  enough.data <- rep(F,length(names(bb)))
  for (l in ds.levels) { 
    mask <- apply(bb[levels==l,,drop=F], 2, function(x) { sum(complete.cases(x)) > min.sampling.support } )
    names(mask) <- NULL
    if (length(enough.data) != length(mask)) otLog("\n\nERROR1! Send email to andreas@maunz.de that you have seen this error!\n\n")
    enough.data <- enough.data | mask
  }
  n.stripped.mss <<- dim(bb)[2]-sum(enough.data)
  n.kept <- sum(enough.data)
  if (do.ot.log) otLog(paste("Stripped",n.stripped.mss,"patterns, kept",n.kept))
  bb <- bb[,enough.data]
  bb$levels <- levels
  

  # Get chi-square
  if (do.ot.log) otLog("Chisq")
  ans.patterns <<- c()
  ans.p.values <<- c()
  for (p in names(bb)[names(bb) != "levels"]) {
    sp <- rep(0,length(bb$levels) / length(ds.levels))
    for (l in ds.levels) sp <- sp + bb[bb$levels==l,p]
    sp <- sp[complete.cases(sp)] 
    sp.zero <- sp == 0 # AM: remove zero entries from failed matches
    if (do.ot.log) if (sum(sp.zero>0)) otLog(paste("Pattern",p,": removed",sum(sp.zero),"zero matches"))
    sp <- sp[!sp.zero]
    dsp <- prop.table(table(sp)) # Categorical distribution for p
    sp.values <- as.numeric(names(dsp))
    chisqv <- 0.0
    for (l in ds.levels) {
      spl <- bb[bb$levels==l,p]
      spl <- spl[complete.cases(spl)]
      spl <- spl[!sp.zero]
      for (spv in sp.values) {
        chisqv <- chisqv + dsp[as.character(spv)] * rectIntegrate(f=wSquaredErr,from=0,to=10000,lambda=mean(spl[sp==spv]),spx=(ds.table[l]*spv/ds.n))
      }
    }
    ans.patterns <<- c(ans.patterns, p)
    ans.p.values <<- c(ans.p.values, pchisq(chisqv,length(ds.levels)-1))

  }

  # if (do.ot.log) otLog(ans.patterns)
  # if (do.ot.log) otLog(ans.p.values)
  n.stripped.cst <<- sum(is.nan(ans.p.values) | ans.p.values<=0.95)
  n.kept <- sum(!is.nan(ans.p.values) & ans.p.values>0.95)
  if (do.ot.log) otLog(paste("Stripped",n.stripped.cst,"patterns, kept",n.kept))
  ans.patterns <<- ans.patterns[!is.nan(ans.p.values) & ans.p.values>0.95]
  ans.p.values <<- ans.p.values[!is.nan(ans.p.values) & ans.p.values>0.95]

  if (do.ot.log) otLog("Done")
}

# merges 2nd list to 1st and returns result
# suitable for .combine
mergeLists = function(x, xn, levels=length(ds.levels)) {
  
  if (is.null(merge.time.start)) merge.time.start <<- Sys.time()
  
  if (do.ot.log) otLog("Merging")

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


# get a single sample feature row
sampleFeatureRow = function(x,fds) {
  ans = as.matrix( fds[ fds["SMILES"]==x["SMILES"],,drop=F ][1,] )[1,] # as.matrix and [1,] drops 'ans' to named vector
  ans[is.na(ans)] = 0
  ans["SMILES"]=x["SMILES"]
  ans
}

# get per-factor support in x from y, with levels from 'levels'
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
postDataset = function(x,server,tempFilePrefix="R_") {
  tf <- tempfile(pattern=paste(tempFilePrefix, Sys.getpid(), sep=""))
  if (do.ot.log) otLog(paste(gsub(".*boot", "boot", gsub(gsub(".*_",  "_", tf),"", tf))))
  tryCatch({
    write.csv(x=x, file=tf, row.names=F, quote=F, na='')
    task <- postForm(server, file=fileUpload(filename=tf, contentType="text/csv"), .checkParams=F )    
  }, finally = {
    unlink(tf)
  })
  task
}

# Emit a log message
otLog = function(text) {
  cat(paste("D,",format(Sys.time(), "[%Y-%m-%dT%H:%M:%S]"),"R DEBUG -- : bbrc-sample        ::", text, "\n"))
}
