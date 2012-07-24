# require package xtable to create latex table
packs=c("xtable","lattice")
for (pack in packs) {
  if (!suppressPackageStartupMessages(suppressWarnings(require(pack,character.only=T)))) {
      install.packages(pack)
    suppressPackageStartupMessages(require(pack,character.only=T))
  }
}

# library
source("anal-lib.R")


#' Analyse the results of bbrc-sample runs
#' The following measures are assessed:
#' Error measures E1 and E2
#' The number of subgraphs returned
#' The number of subgraphs stripped by either minimum sample size
#' or chi-squared test
#' @param assays Specify an assay vector
#' @param outputFile Specify Latex output file
#' @return Produces a latex table
#' @example {
#'   anal(assays=c("SAL", "RAT", "MCC", "KAZ"), outputFile="anal.tex")
#' }
anal <- function (assays=c(), outputFile="anal.tex") {
  res = NULL
  for (i in 1:length(assays)) {
     add = NULL
     
     # Errors
     for (j in 1:2) { 
       file <- paste(assays[i],"_E",j,".csv",sep="")
       data <- na.omit( read.csv( file ) )
       meanData <- data.frame(apply(data,2,function(x) paste(sprintf("%.3f",mean(x)),"(",sprintf("%.4f",sd(x)),")")))
       add <- if (is.null(add)) meanData else cbind(add,meanData)
     }
  
     # Nr features
     file <- paste(assays[i],"_nr_features",".csv",sep="")
     data <- read.csv( file )
     meanData <- data.frame(apply(data,2,mean))
     add <- if (is.null(add)) meanData else cbind(add,meanData)
  
     # Nr stripped mss
     file <- paste(assays[i],"_n_stripped_mss",".csv",sep="")
     data <- read.csv( file )
     meanData <- data.frame(apply(data,2,mean))
     add <- if (is.null(add)) meanData else cbind(add,meanData)
  
     # Nr stripped cst
     file <- paste(assays[i],"_n_stripped_cst",".csv",sep="")
     data <- read.csv( file )
     meanData <- data.frame(apply(data,2,mean))
     add <- if (is.null(add)) meanData else cbind(add,meanData)
  
     assayNames=matrix(rep(assays[i],3),3,1)
     methodNames=rownames(add)
     add <- cbind(assayNames, methodNames, add)
     rownames(add)=NULL
     res = if (is.null(res)) add else rbind(res,add)
  }
  
  if (!is.null(res)) {
    colnames(res) = c("Assay", "Method", "E1", "E2", "Subgraphs", "Stripped MSS", "Stripped CST")
    print(
          xtable(
                 res,
                 digits=c(3,3,3,3,3,1,1,1),
                 label="t:anal",
                 caption="Bias and Accuracy",
                 ), 
          file=outputFile,
          table.placement="t"
         )
  }
  else {
    stop("result not set")
  }

}


#' Generate significance test table
#' @param assays data to consider
#' @param pairsList list of pairs of methods to compare
#' @param alpha significance level
#' @param outputFile Latex code
#' @example {
#'  tests (assays=c("SAL", "RAT", "MCC", "KAZ"))
#' }
tests <- function(assays, pairsList=list(c("MLE","BBRC"),c("MEAN","BBRC"),c("MLE","MEAN")), alpha=0.0001, outputFile="sign.tex") {
  res = NULL
  if (length(assays)>0) {
    results=list()
    for (i in 1:length(assays)) {
      assayRes=list()
      for (error in c("E1", "E2")) {
        errorRes=c()
        file <- paste(assays[i],"_",error,".csv",sep="")
        data <- na.omit ( read.csv( file ) )
        for (pairs in pairsList) {
          add <- runPairedTests(data[,pairs[1]], data[,pairs[2]], alpha)
          errorRes<-c(errorRes, add)
        }
        assayRes[[error]]=errorRes
      }
      results[[assays[i]]] = assayRes
    }
    res=data.frame(results)
    rowNames=unlist(lapply(pairsList,function(x) paste(x,collapse=" vs. ")))
    row.names(res) = rowNames
    names(res) = sapply(names(res), function(x) gsub("\\."," ", x))

  }
  if (!is.null(res)) {
    print(
          xtable(
                 res,
                 label="t:sign",
                 caption=paste("Significance ($p$=", alpha, ")" , sep="")
                 ), 
          file=outputFile,
          table.placement="t"
         )
  }
}


#' Generate boxplots for assays
#' Plots are stacked by methods
#' Plotgroups are stacked by assays
#' @param assays data to consider (methods are fixed)
#' @param error Measure to plot
#' @return bwplot from lattice package
#' @example {
#'   plots (assays=c("SAL", "RAT", "MCC", "KAZ"))
#' }
plots <- function(assays, error="E1", layout=c(1,length(assays))) {
  res = NULL
  methods=c("MLE", "MEAN", "BBRC")
  if (length(assays)>0) {
    results=NULL
    for (i in 1:length(assays)) {
      assayRes=list()
      file <- paste(assays[i],"_",error,".csv",sep="")
      data <- na.omit( read.csv( file ) )
      errorRes=NULL
      mLabls=NULL
      mCLabls=NULL
      for (method in methods) {
        if (is.null(errorRes)) errorRes=data[,method]
        else errorRes=c(errorRes,data[,method])
        if (is.null(mLabls)) mLabls=rep(method,NROW(data))
        else mLabls=c(mLabls,rep(method,NROW(data)))
        if (is.null(mCLabls)) mCLabls=rep(assays[i],NROW(data))
        else mCLabls=c(mCLabls,rep(assays[i],NROW(data)))
      }
      errorRes=data.frame(plotCollectionLabel=mCLabls,plotLabel=mLabls,values=errorRes)
      
      if (is.null(results)) results=errorRes
      else results=rbind(results,errorRes)
    }
    res=results
  }
  if (!is.null(res)) {
    bpCollection( data=res, layout=layout, xlab=error )
  }
}


#' Main
tests (assays=c("SAL", "RAT", "MCC", "KAZ"))
anal  (assays=c("SAL", "RAT", "MCC", "KAZ"))
postscript(file="bp.eps",horizontal=F,paper="special",width=12, height=5)
plots (assays=c("SAL", "RAT", "MCC", "KAZ"), layout=c(2,2))
dev.off()
