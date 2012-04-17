library('RCurl')
library('doMC')
registerDoMC()

SERVER="http://toxcreate3.in-silico.ch:8082"
BBRC=paste(SERVER,"algorithm","fminer","bbrc",sep='/')
num.boots=100


# # # Main function

bootBbrc = function(datasetUri, numboots=num.boots) {

  # load dataset
  ds = getdataset(datasetUri)
  n = dim(ds)[1]

  # main loop
  foreach(j=1:numboots, .combine=mergelists) %dopar% {
    idx=drawsample(1:n)
    sample=ds[idx,]
    task=postdataset(sample,tfprefix=paste("boot_bbrc_sample",j,sep=""))
    while(1) { sampleUri = gettask(task); if (sampleUri != task) break; }
    task=postrequest(BBRC,list(dataset_uri=sampleUri))
    while(1) { sampleFeaturesUri = gettask(task); if (sampleFeaturesUri != task) break; }
    sampleFeatures = getdataset(sampleFeaturesUri)
    sampleFeatures["SMILES"]=NULL
    apply(sampleFeatures,2,function(x) sum(x))
    as.list(sampleFeatures)
  }

}

# merges 2nd list to 1st and returns result
# suitable for .combine
mergelists = function(x,xn) {
  for(idx in names(xn)) { x[[idx]] = c( x[[idx]], xn[[idx]] ) }
  x
}

# draw a bootstrap sample
drawsample = function(x,size=length(x)) {
  sample(x,size=size,replace=T)
}

# get task
gettask = function(uri) {
  Sys.sleep(0.75)
  uri = getrequest(uri,"text/uri-list")
  #substr(uri, 1, nchar(uri)-1) # remove trailing newline
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
postdataset = function(x,tfprefix="R_") {
  tf=tempfile(pattern=paste(tfprefix, Sys.getpid(), sep=""))
  print(tf)
  tryCatch({
    write.csv(x=x,file=tf,row.names=F,quote=F,na='')
    task=postForm(paste(SERVER,"dataset",sep='/'), file=fileUpload(filename=tf,contentType="text/csv"), .checkParams=F )
  }, finally = {
    #unlink(tf)
  })
  task
}


# # # start
bootBbrc("http://toxcreate3.in-silico.ch:8082/dataset/923",1)
