sum(installed.packages()[,1]=="RCurl")==0 && (install.packages('RCurl')||1)
sum(installed.packages()[,1]=="doMC")==0 && (install.packages('doMC')||1)

suppressPackageStartupMessages(library('RCurl'))
suppressPackageStartupMessages(library('doMC'))
registerDoMC()

source("bbrc-sample-psqerr.R")

SERVER="http://toxcreate3.in-silico.ch:8082"
HAMSTER="http://toxcreate3.in-silico.ch:8082/dataset/5092"
BBRC=paste(SERVER,"algorithm","fminer","bbrc",sep='/')
BOOTS=100
MIN_FREQUENCY_PER_SAMPLE=8
MIN_SAMPLING_SUPPORT=10


# # # Main function

bootBbrc = function(datasetUri, numboots=BOOTS, min_frequency_per_sample=MIN_FREQUENCY_PER_SAMPLE, min_sampling_support=MIN_SAMPLING_SUPPORT,del=NULL) {

  set.seed(1)

  # load dataset
  ds = getdataset(datasetUri)
  ds_factors = factor(ds[,2])
  ds_levels <<- levels(ds_factors)
  ds_table = table(ds_factors)

  n = dim(ds)[1]
  cat(paste("Full set size:",n,"\n"))

  # main loop
  if(is.null(del)) {
    bb = foreach(j=1:numboots, .combine=mergelists) %dopar% {
    #for (j in 1:numboots) {
      idx = drawsample(1:n)
      sample = ds[idx,]
  
      task = postdataset(sample,tempfileprefix=paste("boot_bbrc_sample_",j,"_",sep=""))
      sampleUri = getresult(task)
  
      task = postrequest(BBRC,list(dataset_uri=sampleUri,min_frequency_per_sample=as.character(min_frequency_per_sample)))
      sampleFeaturesUri = getresult(task)
  
      sampleFeatures = getdataset(sampleFeaturesUri)
      sampleFeatures["SMILES"]=NULL
      nr_features=dim(sampleFeatures)[2]
      class_support = apply(sampleFeatures,2,function(x) supportperfactor(x,sample[,2],ds_levels))
  
      deleterequest(sampleFeaturesUri)
      deleterequest(sampleUri)
     
      as.list(as.data.frame(class_support))
    }
    bb=data.frame(bb,check.names=F)
    bb$levels=rep(ds_levels, dim(bb)[1]/length(ds_levels))
  }
  else
    bb=del
  
  # Filter patterns with enough sampling support
  cat("\nFiltering\n")
  enough_data=rep(F,length(names(bb)))
  for (l in ds_levels) { 
    mask = apply(bb[bb$levels==l,,drop=F], 2, function(x) { sum(complete.cases(x))>min_sampling_support } )
    names(mask) = NULL
    enough_data = enough_data | mask
  }
  cat(paste("Stripped",dim(bb)[2]-sum(enough_data),"patterns, kept",sum(enough_data),"\n"))
  bb=bb[,enough_data]
  

  # Get chi-square
  cat("\nChisq\n")
  for (p in names(bb)[names(bb)!="levels"]) {
    cat("\n")
    sp = rep(0,length(bb$levels)/length(ds_levels))
    for (l in ds_levels) {
      sp = sp + bb[bb$levels==l,p]
    }
    sp = sp[complete.cases(sp)] 
    chisq=0.0
    for (l in ds_levels) {
      spl = bb[bb$levels==l,p]
      spl = spl[complete.cases(spl)]
      p_weight = mean(sp)/n # weight for p
      spx = ds_table[l] * p_weight # eXpected support for p
      chisq = chisq + rectintegrate(f=wsquarederr,from=0,to=10000,lambda=mean(spl),spx=spx) / spx
    }
    cat(paste(p,":",chisq,"\n"))
  }
  cat("\nDone\n")
  bb

}

# merges 2nd list to 1st and returns result
# suitable for .combine
mergelists = function(x,xn,levels=length(ds_levels)) {

  x=data.frame(x,check.names=F)
  xn=data.frame(xn,check.names=F)

  xus = names(x)[names(x) %in% names(xn)]
  xnn = names(xn)[!names(xn) %in% names(x)]

  if (length(xus)>0) {
    x_bel = data.frame(matrix(NA,levels,dim(x)[2]))
    names(x_bel) = names(x)
    x_bel[1:levels,xus] = xn[1:levels,xus]
    x = rbind(x,x_bel)
  }
  
  if (length(xnn)>0) {
    x_rig = data.frame(matrix(NA,dim(x)[1],length(xnn)))
    names(x_rig) = xnn
    x_rig[c(dim(x)[1]-levels+1,dim(x)[1]), xnn] = xn[1:levels, xnn]
    x = cbind(x,x_rig)
  }

  as.list(x)

}

# draw a bootstrap sample
drawsample = function(x,size=length(x)) {
  sample(x,size=size,replace=T)
}

# get dataset from task
getresult= function(uri) {
  while(1) { 
    Sys.sleep(0.5)
    uri = getrequest(uri,"text/uri-list")
    if (length(grep("dataset",uri))>0) break
  }
  uri
}

# get per-factor support in x from y, with levels from z
supportperfactor= function(x,y,levels) {
  pattern = factor(y,levels=levels)
  table(pattern[x>0]) # table sums across the factors
}


# # # REST library

# GET
getrequest = function(uri,accept) {
  getURL(uri,httpheader=c(Accept=accept))
}

# POST (assumes a task is returned as uri-list)
postrequest = function(uri,params) {
  postForm(uri,.params=params)
}

# DELETE a dataset by URI
deleterequest = function(uri) {
  httpDELETE(uri)
}

# GET a dataset by URI
getdataset = function(uri) {
  read.csv(paste(uri,".csv",sep=""),header=T,check.names=F,stringsAsFactors=F)
}

# POST a dataset (CSV)
postdataset = function(x,tempfileprefix="R_") {
  tf=tempfile(pattern=paste(tempfileprefix, Sys.getpid(), sep=""))
  cat(paste(gsub(".*boot","boot",gsub(gsub(".*_", "_", tf),"",tf)),"\n"))
  tryCatch({
    write.csv(x=x,file=tf,row.names=F,quote=F,na='')
    task=postForm(paste(SERVER,"dataset",sep='/'), file=fileUpload(filename=tf,contentType="text/csv"), .checkParams=F )
  }, finally = {
    unlink(tf)
  })
  task
}
