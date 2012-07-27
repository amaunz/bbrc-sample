cpdb=read.csv("CPDB.csv", na.strings="")
mycpdb=data.frame(cpdb$STRUCTURE_SMILES, cpdb$ActivityOutcome_CPDBAS_MultiCellCall, cpdb$ActivityOutcome_CPDBAS_Mutagenicity)
mycpdb=mycpdb[complete.cases(mycpdb),]
newcol=character(length=dim(mycpdb)[1])
newcol[mycpdb[,2]=="active"   & mycpdb[,3]=="active"] = "cymy"
newcol[mycpdb[,2]=="inactive" & mycpdb[,3]=="inactive"] = "cnmn"
newcol[mycpdb[,2]=="active"   & mycpdb[,3]=="inactive"] = "cymn"
newcol[mycpdb[,2]=="inactive" & mycpdb[,3]=="active"] = "cnmy"
mycpdb=data.frame(mycpdb[,1], newcol)
names(mycpdb)=c("SMILES", "CM")
write.csv(mycpdb, file="cpdb_carc_mutag.csv",row.names=F)
