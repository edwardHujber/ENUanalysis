library("seqinr")
library("ggplot2")


upsidedown2 <- function(bp){
   # "G>A" becomes "C>T"
   paste0(  
      basePair(substr(bp,1,1)),">",
      basePair(substr(bp,3,3))
   )
   
}
humanGC <- 0.41

denovo.db.variants.v.1.2 <- read.delim("Z:/programming/R/ENU analysis/denovo-db/denovo-db.variants.v.1.2.tsv", skip=1,header=TRUE)
denovoSNPs <-head(sort(table(denovo.db.variants.v.1.2$Variant),decreasing=TRUE),n=12)

i<-1
while(i <length(denovoSNPs)){
   samePair <- match(upsidedown2(names(denovoSNPs)[i]),names(denovoSNPs))
   if(!is.na(samePair)){
      denovoSNPs[i] <- denovoSNPs[i]+denovoSNPs[samePair]
      denovoSNPs <- denovoSNPs[-samePair]
      i<-0
   }
   i<-i+1
}

gcnormal <- rep(NA,6)
names(gcnormal)<-names(denovoSNPs)
for( i in 1:length(denovoSNPs)){
   if(grepl("A>|T>",names(denovoSNPs[i]))){gcnormal[i] <- (denovoSNPs[i]/(1-humanGC))}
   if(grepl("C>|G>",names(denovoSNPs[i]))){gcnormal[i] <- (denovoSNPs[i]/(humanGC))}
}

  
   
denovoSNPpercents <- gcnormal/sum(gcnormal)

sortByParamter <- AA_ASA$Pi
AAandParameter <- data.frame("AA"=AA_ASA$X[order(sortByParamter)],"param"=sortByParamter[order(sortByParamter)])
AAbyparam<-as.character(AAandParameter$AA)



AAchanges <-data.frame(matrix(NA,length(allAAwStop),length(allAAwStop)),row.names = allAAwStop,stringsAsFactors = FALSE)
colnames(AAchanges) <- allAAwStop
translatedCodons<- apply(as.array(codons),1,function(x){translate(s2c(x))} )
## put Zeros into the cells that are possible with a single nt change
for(i in 1:length(allAAwStop)){
   codonSet <- codons[which(translatedCodons==allAAwStop[i])]
   for(j in 1:length(codonSet)){
      for(k in 1:3){
         thisNT <- substr(codonSet[j],k,k)
         otherNTs <- nt[-which(nt==thisNT)]
         for(l in 1:length(otherNTs)){
            testCodon<-codonSet[j]
            substr(testCodon,k,k)<-otherNTs[l]
            AAchanges[allAAwStop[i],translate(s2c(testCodon))]<-0
         }
      }
   }
}

for(i in 1:length(allAAwStop)){
   codonSet <- codons[which(translatedCodons==allAAwStop[i])]
   for(j in 1:length(codonSet)){
      for(k in 1:3){
         thisNT <- substr(codonSet[j],k,k)
         otherNTs <- nt[-which(nt==thisNT)]
         for(l in 1:length(otherNTs)){
            testCodon<-codonSet[j]
            substr(testCodon,k,k)<-otherNTs[l]
            AAchanges[allAAwStop[i],translate(s2c(testCodon))]<-AAchanges[allAAwStop[i],translate(s2c(testCodon))]+1
         }
      }
   }
}
AAchanges["W","W"]<-0
AAchanges["M","M"]<-0
AAchangePaths<-as.matrix(AAchanges[1:20,1:20])
class(AAchangePaths)<-"numeric"
AAchanges[!is.na(AAchanges)]<-0



# successCount<-0
# while(successCount<100000){
#    start<-sample(codons,size=1)
#    startAA <- translate(s2c(start))
#    position<-sample(1:3,size=1)
#    base2change <- substr(start,position,position)
#    change<-sample(names(denovoSNPpercents),size=1,prob=denovoSNPpercents, replace=TRUE)
#    changeFrom <- substring(change,1,1)
#    changeTo<-substring(change,3,3)
#    if(base2change == changeFrom ){
#       substr(start,position,position)<-changeTo
#       endAA <- translate(s2c(start))
#       AAchanges[startAA,endAA]<-1+as.numeric(AAchanges[startAA,endAA])
#       successCount<-successCount+1
#    }
#    else{
#       if(base2change==basePair(changeFrom)){
#       substr(start,position,position)<-basePair(changeTo)
#       endAA <- translate(s2c(start))
#       AAchanges[startAA,endAA]<-1+as.numeric(AAchanges[startAA,endAA])
#       successCount<-successCount+1
#       }
#       else{
#          AAchanges[startAA,startAA]<-1+as.numeric(AAchanges[startAA,startAA])
#          successCount<-successCount+1
#       }
#    }
# }


CCnums<-matrix(0,dim(CodonChanges)[1],dim(CodonChanges)[2]+1)

for(i in 0:((64*15000)-1)){
   start<-end<-codons[1+i-floor(i/64)*64]
   
   position<-sample(1:3,size=1)
   base2change <- substr(start,position,position)
   change<-sample(names(denovoSNPpercents),size=1,prob=denovoSNPpercents, replace=TRUE)
   changeFrom <- substring(change,1,1)
   changeTo<-substring(change,3,3)
   if(base2change == changeFrom ){
      substr(end,position,position)<-changeTo
      r<-which(rownames(CodonChanges)==start)
      c<-which(CodonChanges[r,]==end)+1
      CCnums[r,c]<-1+CCnums[r,c]

   }
   else{
      if(base2change==basePair(changeFrom)){
         substr(end,position,position)<-basePair(changeTo)
         r<-which(rownames(CodonChanges)==start)
         c<-which(CodonChanges[r,]==end)+1
         CCnums[r,c]<-1+CCnums[r,c]

      }
      else{
         r<-which(rownames(CodonChanges)==start)
         c<-1
         CCnums[r,c]<-1+CCnums[r,c]

      }
   }
}
rownames(CCnums) <-codons
heatmap.2( CCnums[,2:10],
           col = colorRampPalette( c("blue","red","yellow","white"), space="rgb")(100),
           cellnote = CodonChanges,
           notecol="black",
           trace = "none", 
           na.color="gray20",
           Rowv = FALSE,
           Colv = FALSE,
           dendrogram = "none",
           offsetRow = -46,
           offsetCol = -45,
           srtCol = 0,
           main="denovo codon"
)


## translate codons to AA change table
for(i in 1:nrow(CodonChanges)){
   for(j in 1:ncol(CodonChanges)){
      from <- translate(s2c(rownames(CodonChanges[i,])))
      to <- translate(s2c(as.character(CodonChanges[i,j])))
      AAchanges[from,to]<-AAchanges[from,to]+CCnums[i,j+1]
   }
}



AAchanges<-AAchanges[AAbyparam,AAbyparam]
AAchanges_numeric<-as.matrix(AAchanges)
class(AAchanges_numeric)<-"numeric"
print(sum(AAchanges_numeric,na.rm=TRUE))

AAchanges_numeric


## All changes
plot(x = rep(AAandParameter$param,each=20),
     y = rep(AAandParameter$param,20),
     cex = as.numeric(AAchanges_numeric[1:20,1:20])/max(as.numeric(AAchanges_numeric[1:20,1:20]),na.rm=TRUE)*4.9+0.5, 
     pch = 20,
     ylim = rev(range(AAandParameter$param)),
     axes=FALSE)

axis(2, at=AAandParameter$param,labels=AAandParameter$AA,
      las=1, tck=-.01)

axis(3, at=AAandParameter$param,labels=AAandParameter$AA,
     las=1, tck=-.01)

## Missense
AAmissense_numeric<-AAchanges_numeric[1:20,1:20]
diag(AAmissense_numeric)<-NA
plot(x = rep(AAandParameter$param,each=20),
     y = rep(AAandParameter$param,20),
     cex = as.numeric(AAmissense_numeric[1:20,1:20])/max(as.numeric(AAmissense_numeric[1:20,1:20]),na.rm=TRUE)*4.9+0.5, 
     pch = 20,
     ylim = rev(range(AAandParameter$param)),
     axes=FALSE,
     main="Missense by random number",
     xlab="",
     ylab=""
)

axis(2, at=AAandParameter$param,labels=AAandParameter$AA,
     las=1, tck=-.01)

axis(3, at=AAandParameter$param,labels=AAandParameter$AA,
     las=1, tck=-.01)





deltaASA<-matrix(nrow=20,ncol=20,rep(AAandParameter$param,each=20))-matrix(nrow=20,ncol=20,rep(AAandParameter$param,20))

diag(AAchangePaths)<-NA
obsCountTab <- table(rep(as.numeric(deltaASA)[!is.na(AAmissense_numeric)],as.numeric(AAmissense_numeric)[!is.na(AAmissense_numeric)]))

codonCountTab <- table(rep(as.numeric(deltaASA)[!is.na(AAmissense_numeric)],as.numeric(AAchangePaths)[!is.na(AAchangePaths)]))

plot(obsCountTab)
plot(codonCountTab)

log2diff <- log2((obsCountTab/sum(obsCountTab))/(codonCountTab/sum(codonCountTab)))
plot(log2diff,ylim=1.1*range(log2diff))
log2diffDF<-data.frame("ASA"=as.numeric(names(log2diff)),"ratio"=as.numeric(log2diff))
lw1<-loess(ratio~ASA,log2diffDF,span=0.75)

j <- order(log2diffDF$ASA)
lines(log2diffDF$ASA[j],lw1$fitted[j],col="red",lwd=3)


rep(as.numeric(deltaASA)[!is.na(AAchanges_numeric)],as.numeric(AAchanges_numeric)[!is.na(AAchanges_numeric)])


heatmap.2( log(AAchanges_numeric),
           col = colorpanel(100,"red","yellow","green"),
           #cellnote = log(AAchanges_numeric/sum(AAchanges_numeric,na.rm=TRUE)),
           notecol="black",
           trace = "none", 
           na.color="gray20",
           Rowv = FALSE,
           Colv = FALSE,
           dendrogram = "none",
           offsetRow = -46,
           offsetCol = -45,
           srtCol = 0,
           main="log thing"
)

heatmap.2( (AAchanges_numeric),
           col = colorpanel(100,"red","yellow","green"),
           cellnote = round(100*AAchanges_numeric/sum(AAchanges_numeric,na.rm=TRUE),2),
           notecol="black",
           trace = "none", 
           na.color="gray20",
           Rowv = FALSE,
           Colv = FALSE,
           dendrogram = "none",
           offsetRow = -46,
           offsetCol = -45,
           srtCol = 0,
           main="thing"
)



