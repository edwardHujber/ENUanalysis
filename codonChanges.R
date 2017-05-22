
## AA classes
hydrophobic <- c("A","I","L","V")
aromatic <- c("F","W","Y")
polarNeutral <- c("N","C","Q","M","S","T")
acidic <- c("D","E")
basic <- c("R","H","K")
unique <- c("G","P")

allAA <- data.frame("AA" = c(hydrophobic,aromatic,polarNeutral,acidic,basic,unique), 
                    "class"= c(rep("hydrophobic",length(hydrophobic)),
                               rep("aromatic",length(aromatic)),
                               rep("polarNeutral",length(polarNeutral)),
                               rep("acidic",length(acidic)),
                               rep("basic",length(basic)),
                               rep("unique",length(unique))))
charged <- c(polarNeutral,acidic,basic)
uncharged <- c(hydrophobic,aromatic)

## this is the important processing stuff

effectClasses <- data.frame(table(full$V16))
effectClasses <- effectClasses[which(effectClasses$Freq >0 ),]


codonChanges <- data.frame(table(full$V18))
codonChanges <- codonChanges[which(codonChanges$Freq >0 ),]
codonChanges <- cbind(codonChanges,full$V17[match(codonChanges$Var1,full$V18)])
codonChanges <- codonChanges[,c(1,3,2)]
colnames(codonChanges)[1:2] <- c("CodonChange","AAchange")
AAchanges <- t(data.frame(strsplit((as.character(codonChanges$AA)),"/")))
codonChanges <- cbind(codonChanges,AAchanges)
colnames(codonChanges)[4:5] <- c("wt","mut")
codonChanges <- codonChanges[,c(1,2,4,5,3)]

missensecodon <- codonChanges[which(codonChanges$wt != codonChanges$mut & 
                                 codonChanges$wt != "*" &
                                 codonChanges$mut != "*"),]

missense <- full[which(full$V16 == "NON_SYNONYMOUS_CODING"),]

AAchanges <- t(data.frame(strsplit((as.character(missense$V17)),"/")))
colnames(AAchanges) <- c("wt","mut")
missense <- cbind(missense,AAchanges)


conserved <- missense[which(allAA$class[match(missense$wt,allAA$AA)] == allAA$class[match(missense$mut,allAA$AA)]),]
unconserved <- missense[which(allAA$class[match(missense$wt,allAA$AA)] != allAA$class[match(missense$mut,allAA$AA)]),]

unconChangeClasses <- data.frame(allAA$class[match(missense$wt,allAA$AA)],allAA$class[match(missense$mut,allAA$AA)])
colnames(unconChangeClasses)[1:2] <- c("wtClass","mutClass") 
unconChangeClasses <- as.data.frame(table(unconChangeClasses))
# unconChangeClasses <- unconChangeClasses[which(unconChangeClasses$wtClass != unconChangeClasses$mutClass),]
unconChangeClasses <- unconChangeClasses[order(unconChangeClasses$wtClass, decreasing=FALSE),]

