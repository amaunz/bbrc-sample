
runPairedTests <- function(x,y,alpha=.05) {
  add <- ""
  tP <- t.test(x,y,paired=T)$p.value
  wP <- wilcox.test(x,y,paired=T)$p.value
  if (tP < alpha) add <- paste(add,"t",sep="")
  if (wP < alpha) add <- paste(add,"w",sep="")
  add
}


# Produces a layout[0] x layout[1] boxplot-collection
# Each boxplot-collection consists of levels(data$plotLabel) many boxplots
# @param data data frame with columns plotLabel (factor), plotCollectionLabel (factor), values (numeric)
bpCollection <- function ( data=NULL, layout=c(0,0), xlab="xlab") {
  if (!is.null(data)) {
    cols = list(col="black",cex=0.75,fill="gray")
    bwplot( plotLabel ~ values | plotCollectionLabel, 
            data = data, 
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
 
