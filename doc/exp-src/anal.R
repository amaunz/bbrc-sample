library(xtable)
assays <- c("SAL", "MOU", "RAT", "MCC")

res = NULL
for (i in 1:length(assays)) {
   add = NULL
   for (j in 1:2) { 
     file <- paste(assays[i],j,".csv",sep="")
     data <- read.csv( file )
     meanData <- data.frame(apply(data,2,mean))
     add <- if (is.null(add)) meanData else cbind(add,meanData)
   }
   assayNames=matrix(rep(assays[i],3),3,1)
   methodNames=rownames(add)
   add <- cbind(assayNames, methodNames, add)
   rownames(add)=NULL
   res = if (is.null(res)) add else rbind(res,add)
}
colnames(res) = c("Assay", "Method", "E1", "E2")
print(
      xtable(
             res,
             digits=3,
             label="t:anal",
             caption="Analysis of BBRC sampling",
             ), 
      file="anal.tex",
      table.placement="t"
     )

