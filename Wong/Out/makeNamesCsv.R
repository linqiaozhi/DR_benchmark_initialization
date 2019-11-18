
if (!require(flowCore)) { 
  source("http://bioconductor.org/biocLite.R")
  biocLite("flowCore")
} 


inFileN <- "test.fcs"        #FCS file to make names file for 
outNamesCsv <- "names.csv"   #Output csv file name 





FF <- read.FCS(inFileN)

#Fixup column names
colNames <- FF@parameters$desc
empties <- which(is.na(colNames) | colNames== " ")
colNames[empties] <- FF@parameters$name[empties]
write.csv(colNames, outNamesCsv, row.names = F)