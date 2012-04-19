suppressPackageStartupMessages(library('RCurl'))
suppressPackageStartupMessages(library('doMC'))
registerDoMC()

SERVER="http://toxcreate3.in-silico.ch:8082"
BBRC=paste(SERVER,"algorithm","fminer","bbrc",sep='/')
BOOTS=100


# # # Main function

bootBbrc = function(datasetUri, numboots=BOOTS) {

  # load dataset
  ds = getdataset(datasetUri)
  n = dim(ds)[1]
  cat(paste("Full set size:",n,"\n"))

  # main loop
  bb = foreach(j=1:numboots, .combine=mergelists) %dopar% {
    idx=drawsample(1:n)
    sample=ds[idx,]

    task=postdataset(sample,tempfileprefix=paste("boot_bbrc_sample_",j,"_",sep=""))
    sampleUri = getresult(task)

    task=postrequest(BBRC,list(dataset_uri=sampleUri))
    sampleFeaturesUri = getresult(task)
    print(sampleFeaturesUri)

    sampleFeatures = getdataset(sampleFeaturesUri)
    print(dim(sampleFeatures))
    sampleFeatures["SMILES"]=NULL
    sampleFeatures=apply(sampleFeatures,2,function(x) sum(x))
    as.list(sampleFeatures)
  }
  bb

}

# merges 2nd list to 1st and returns result
# suitable for .combine
mergelists = function(x,xn) {
  padlen = length(x[[1]])
  print("PAD UNK")
  for (n in names(x)[!names(x) %in% names(xn)])  xn[[n]] = 0 
  print("PAD NEW")
  for (n in names(xn)[!names(xn) %in% names(x)]) xn[[n]] = c(rep(0,padlen), xn[[n]])
  print("MERGE")
  for (idx in names(xn)) { x[[idx]] = c( x[[idx]], xn[[idx]] ) } 
  print("DONE")
  x
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




# # # REST library

# GET
getrequest = function(uri,accept) {
  getURL(uri,httpheader=c(Accept=accept))
}

# POST (assumes a task is returned as uri-list)
postrequest = function(uri,params) {
  postForm(uri,.params=params)
}

# GET a dataset by URI
getdataset = function(uri) {
  read.csv(paste(uri,".csv",sep=""),header=T,check.names=F,stringsAsFactors=F)
}

# POST a dataset (CSV)
postdataset = function(x,tempfileprefix="R_") {
  tf=tempfile(pattern=paste(tempfileprefix, Sys.getpid(), sep=""))
  print(tf)
  tryCatch({
    write.csv(x=x,file=tf,row.names=F,quote=F,na='')
    task=postForm(paste(SERVER,"dataset",sep='/'), file=fileUpload(filename=tf,contentType="text/csv"), .checkParams=F )
  }, finally = {
    unlink(tf)
  })
  task
}


# # # start
bootBbrc("http://toxcreate3.in-silico.ch:8082/dataset/923",numboots=100)
