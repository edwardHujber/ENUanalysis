

missense <- read.delim("~/Lab/R/ENU analysis/guntherAAchanges.txt", header=FALSE)

missense <- data.frame(as.matrix(sapply(missense,toupper),ncol=2))


colnames(missense) <- c("wt","mut")


conserved <- missense[which(allAA$class[match(missense$wt,allAA$AA)] == allAA$class[match(missense$mut,allAA$AA)]),]
unconserved <- missense[which(allAA$class[match(missense$wt,allAA$AA)] != allAA$class[match(missense$mut,allAA$AA)]),]

unconChangeClasses <- data.frame(allAA$class[match(missense$wt,allAA$AA)],allAA$class[match(missense$mut,allAA$AA)])
colnames(unconChangeClasses)[1:2] <- c("wtClass","mutClass") 
unconChangeClasses <- as.data.frame(table(unconChangeClasses))
# unconChangeClasses <- unconChangeClasses[which(unconChangeClasses$wtClass != unconChangeClasses$mutClass),]
unconChangeClasses <- unconChangeClasses[order(unconChangeClasses$wtClass, decreasing=FALSE),]

