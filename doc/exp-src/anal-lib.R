
runPairedTests <- function(x,y,alpha=.05) {
  add <- ""
  tP <- t.test(x,y,paired=T)$p.value
  wP <- wilcox.test(x,y,paired=T)$p.value
  if (!is.nan(tP) && tP < alpha) add <- paste(add,"t",sep="")
  if (!is.nan(wP) && wP < alpha) add <- paste(add,"w",sep="")
  add
}


# Produces a layout[0] x layout[1] boxplot-collection
# Each boxplot-collection consists of levels(data$plotLabel) many boxplots
# @param data data frame with columns plotLabel (factor), plotCollectionLabel (factor), values (numeric)
bpCollection <- function ( bpdata=NULL, layout=c(0,0), xlab="xlab") {
  if (!is.null(bpdata)) {
    cols = list(col="black",cex=0.75,fill="gray")
    bwplot( plotLabel ~ values | plotCollectionLabel, 
            data = bpdata, 
            xlab=xlab,
            layout=layout, 
            pch="|", 
            scales=list(alternating=3),
            col=c("black"),
            par.settings = list(
              plot.symbol=cols,
              box.rectangle =  cols,
              box.umbrella = cols,
              box.symbol = cols
            )
          )
  }
}
 
assaySizes <- function () {
  list(KAZ=4069, RAT=1128, MCC=1051, MOU=914, SAL=800, MUL=678, INT=458)
}

