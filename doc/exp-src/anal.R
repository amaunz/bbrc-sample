# require package xtable to create latex table
if (!suppressPackageStartupMessages(suppressWarnings(require("xtable",character.only=T)))) {
    install.packages("xtable")
  suppressPackageStartupMessages(require("xtable",character.only=T))
}


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
       data <- read.csv( file )
       meanData <- data.frame(apply(data,2,mean))
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


runPairedTests <- function(x,y,alpha=.05) {
  add <- ""
  tP <- t.test(x,y,paired=T)$p.value
  wP <- wilcox.test(x,y,paired=T)$p.value
  if (tP < alpha) add <- paste(add,"t",sep="")
  if (wP < alpha) add <- paste(add,"w",sep="")
  add
}


tests <- function(assays, pairsList=list(c("MLE","BBRC"),c("MEAN","BBRC"),c("MLE","MEAN")), alpha=0.001, outputFile="sign.tex") {
  res = NULL
  if (length(assays)>0) {
    results=list()
    for (i in 1:length(assays)) {
      assayRes=list()
      for (error in c("E1", "E2")) {
        errorRes=c()
        file <- paste(assays[i],"_",error,".csv",sep="")
        data <- read.csv( file )
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

tests(assays=c("SAL", "RAT", "MCC", "KAZ"))
anal(assays=c("SAL", "RAT", "MCC", "KAZ"))
