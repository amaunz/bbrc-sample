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
#' @param dir Directory to process
#' @return Produces a latex table
#' @example {
#'   anal(assays=c("SAL", "RAT", "MCC", "KAZ"), outputFile="anal.tex")
#' }
anal <- function (assays=c(), outputFile="anal.tex", dir) {
  res = NULL
  for (i in 1:length(assays)) {
     add = NULL
     
     # Errors
     for (j in 1:5) { 
       file <- paste(dir,"/",assays[i],"/",assays[i],"_E",j,".csv",sep="")
       data <- na.omit( read.csv( file ) )
       meanData <- data.frame(apply(data,2,function(x) paste(sprintf("%.4f",mean(x))," (",sprintf("%.4f",sd(x)),")",sep="")))
       add <- if (is.null(add)) meanData else cbind(add,meanData)
     }
     add$AssaySize=assaySizes()[[assays[i]]]
  
     # Nr features
     # SAL_bbrc_ds_nr_com.csv  SAL_bbrc_ds_nr_f.csv
     file <- paste(dir,"/",assays[i],"/",assays[i],"_","bbrc_ds_nr_f.csv",sep="")
     data <- read.csv( file )
     meanData <- data.frame(apply(data,2,function(x) mean(x,na.rm=T)))
     add <- if (is.null(add)) meanData else cbind(add,meanData)
  
     ## Nr stripped mss
     #file <- paste(dir,"/",assays[i],"/",assays[i],"_n_stripped_mss",".csv",sep="")
     #data <- read.csv( file )
     #meanData <- data.frame(apply(data,2,function(x) mean(x,na.rm=T)))
     #add <- if (is.null(add)) meanData else cbind(add,meanData)
  
     ## Nr stripped cst
     #file <- paste(dir,"/",assays[i],"/",assays[i],"_n_stripped_cst",".csv",sep="")
     #data <- read.csv( file )
     #meanData <- data.frame(apply(data,2,function(x) mean(x,na.rm=T)))
     #add <- if (is.null(add)) meanData else cbind(add,meanData)

     assayNames=matrix(rep(assays[i],3),3,1)
     methodNames=rownames(add)
     add <- cbind(assayNames, methodNames, add)
     rownames(add)=NULL
     res = if (is.null(res)) add else rbind(res,add)
  }


  
  if (!is.null(res)) {
    
    colnames(res) = c("Assay", "Method", "E1", "E2", "E3", "E4", "E5", "AssaySize", "Subgraphs")
    res=res[with(res, order(AssaySize)), ]
    res$AssaySize=NULL
    rownames(res)=NULL

    print(
          xtable(
                 res,
                 label="t:anal",
                 caption="Bias and accuracy",
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
#' @param dir Directory to process
#' @example {
#'  tests (assays=c("SAL", "RAT", "MCC", "KAZ"))
#' }
tests <- function(assays, errors=c("E1", "E2", "E3", "E4", "E5"), pairsList=list(c("MLE","BBRC"),c("MEAN","BBRC"),c("MLE","MEAN")), alpha=0.0001, outputFile="sign.tex", dir) {
  res = NULL
  if (length(assays)>0) {
    results=list()
    for (i in 1:length(assays)) {
      assayRes=list()
      for (error in errors) {
        errorRes=c()
        file <- paste(dir,"/",assays[i],"/",assays[i],"_",error,".csv",sep="")
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

    resNew<-NULL
    for (error in errors) {
      resBroken<-NULL
      for (assay in assays) {
        colName<-paste(assay,error,sep=".")
        if (is.null(resBroken))  { 
          resBroken <- res[,colName,drop=F]
          names(resBroken) = assay  
        } 
        else { 
          rBNames = c(names(resBroken), assay)
          resBroken <- cbind(resBroken, res[,colName,drop=F])
          names(resBroken) = rBNames
        } 
      }
      row.names(resBroken) <- paste(row.names(res), error)
      if (is.null(resNew)) resNew<-resBroken else resNew<-rbind(resNew,resBroken)
    }
    res<-resNew
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
#' @param dir Directory to process
#' @return bwplot from lattice package
#' @example {
#'   plots (assays=c("SAL", "RAT", "MCC", "KAZ"))
#' }
boxplots <- function(assays, error="E1", layout=c(1,length(assays)), dir) {
  res = NULL
  methods=c("MLE", "MEAN", "BBRC")
  if (length(assays)>0) {
    results=NULL
    for (i in 1:length(assays)) {
      assayRes=list()
      file <- paste(dir,"/",assays[i],"/",assays[i],"_",error,".csv",sep="")
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
    bpCollection( bpdata=res, layout=layout, xlab=error )
  }
}

#' Generate boxplots for assays
#' Plots are stacked by methods
#' Plotgroups are stacked by assays
#' @param assays data to consider (methods are fixed)
#' @param error Measure to plot
#' @param dir Directory to process
#' @return bwplot from lattice package
#' @example {
#'   plots (assays=c("SAL", "RAT", "MCC", "KAZ"))
#' }
lineplots <- function(assays, error="E1", dir, yOffset) {

  if (length(assays)>0) {
    results=NULL
    for (i in 1:length(assays)) {
      file <- paste(dir,"/",assays[i],"/",assays[i],"_",error,".csv",sep="")
      data <- na.omit( read.csv( file ) )
      meanData <- data.frame(apply(data,2,function(x) mean(x)))
      #sdData <- data.frame(apply(data,2,function(x) sd(x)))
      #myData <- cbind(meanData,sdData)

      myData <- meanData
      #names(myData) <- c(paste("Mean ",error,sep=""), 
      #                   paste("SD ",error,sep=""))
      names(myData) <- paste("Mean ",error,sep="")
      myData$Assay = rep(assays[i],dim(myData)[1])
      myData$AssaySize = rep(assaySizes()[[assays[i]]], dim(myData)[1])
      myData$Method = row.names(myData)
      if (is.null(results)) results=myData else results=rbind(results,myData)
    }

    row.names(results)=NULL

    #max_y=max(
    #          c(results[[paste("Mean ",error,sep="")]], 
    #            results[[paste("SD ",error,sep="")]]
    #           ))
    max_y=max(results[[paste("Mean ",error,sep="")]])

    #min_y=min(
    #          c(results[[paste("Mean ",error,sep="")]], 
    #            results[[paste("SD ",error,sep="")]]
    #           ))
    min_y=min(results[[paste("Mean ",error,sep="")]])


    resultsMLE <- results[results$Method=="MLE",]
    resultsMLE <- resultsMLE[with(resultsMLE, order(AssaySize)), ]
    resultsMLE$AssaySize <- log(resultsMLE$AssaySize)

    resultsMEAN <- results[results$Method=="MEAN",]
    resultsMEAN <- resultsMEAN[with(resultsMEAN, order(AssaySize)), ]
 
    resultsBBRC <- results[results$Method=="BBRC",]
    resultsBBRC <- resultsBBRC[with(resultsBBRC, order(AssaySize)), ]

    max_x=max(resultsMLE$AssaySize)
    plot_colors=c('blue','red','forestgreen')
    plot_points=c(1,0,2)
    par(mar=c(4.2, 3.8, 0.2, 0.2))

    plot(resultsMLE$AssaySize,resultsMLE[[paste("Mean ",error,sep="")]],type='o',col=plot_colors[1],axes=F,ann=F,ylim=c(min_y,max_y),lty=2,pch=plot_points[1])
    axis(2, las=1, cex.axis=1.0)
    axis(1, lab=F, at=resultsMLE$AssaySize)
    text(x=resultsMLE$AssaySize, y=min_y-yOffset, srt=90, adj=1, labels=paste(resultsMLE$Assay, sep=""),xpd=T, cex=1.0)
    box()
    #lines(resultsMLE$AssaySize,resultsMLE[[paste("SD ",error,sep="")]],type='o',col=plot_colors[2],pch=22,lty=2)

    lines(resultsMLE$AssaySize,resultsMEAN[[paste("Mean ",error,sep="")]],type='o',col=plot_colors[2],lty=2,pch=plot_points[2])
    #lines(resultsMLE$AssaySize,resultsMEAN[[paste("SD ",error,sep="")]],type='o',col=plot_colors[4],pch=22,lty=2)

    lines(resultsMLE$AssaySize,resultsBBRC[[paste("Mean ",error,sep="")]],type='o',col=plot_colors[3],lty=2,pch=plot_points[3])
    #lines(resultsMLE$AssaySize,resultsBBRC[[paste("SD ",error,sep="")]],type='o',col=plot_colors[6],pch=22,lty=2)

    title(xlab=error,ylab="Error")
    #legend(max_x+log(0.5), max_y, c(
    #                                  paste("Mean MLE",sep=""),
    #                                  paste("SD MLE",sep=""),
    #                                  paste("Mean MEAN",sep=""),
    #                                  paste("SD MEAN",sep=""),
    #                                  paste("Mean BBRC",sep=""),
    #                                  paste("SD BBRC",sep="")
    #                                ), cex=1.0, col=plot_colors, pch=21:22, lty=1:2, bty="n")

    legend(max_x+log(0.5), max_y, c(
                                      paste("MLE",sep=""),
                                      paste("MEAN",sep=""),
                                      paste("BBRC",sep="")
                                    ), cex=1.0, col=plot_colors, pch=plot_points, lty=2, bty="n")



  }

}



#' Main

dir="exp12"
assays=c("INT", "MCC", "RAT", "MUL", "KAZ", "MOU")
alpha=0.01

# statistical tests and comparison table
tests (assays=assays, errors=c("E1","E2","E3","E4","E5"), dir=dir, alpha=alpha)
#anal  (assays=assays, dir=dir)


# boxplots (for-loop produces empty plots for reasons unknown)
postscript(file=paste("bp1.eps",sep=""),horizontal=F,paper="special",width=8, height=5)
boxplots (assays=assays, error="E1",layout=c(3,2), dir=dir)
dev.off()

postscript(file=paste("bp2.eps",sep=""),horizontal=F,paper="special",width=8, height=5)
boxplots (assays=assays, error="E2",layout=c(3,2), dir=dir)
dev.off()

postscript(file=paste("bp3.eps",sep=""),horizontal=F,paper="special",width=8, height=5)
boxplots (assays=assays, error="E3",layout=c(3,2), dir=dir)
dev.off()

postscript(file=paste("bp4.eps",sep=""),horizontal=F,paper="special",width=8, height=5)
boxplots (assays=assays, error="E4",layout=c(3,2), dir=dir)
dev.off()

postscript(file=paste("bp5.eps",sep=""),horizontal=F,paper="special",width=8, height=5)
boxplots (assays=assays, error="E5",layout=c(3,2), dir=dir)
dev.off()


# lineplots, ordered by dataset size (for-loop works)
for (e in seq(1,5)) {
  yOffset=0.012
  if (e==4) yOffset = 0.056
  if (e==5) yOffset = 0.07
  #if (e==1) yOffset = 0.016
  #if (e==2 || e==3) yOffset = 0.014
  postscript(file=paste("lp",e,".eps",sep=""),horizontal=F,paper="special",width=4, height=3)
  lineplots (assays=assays, error=paste("E",e,sep=""), dir=dir, yOffset=yOffset)
  dev.off()
}
