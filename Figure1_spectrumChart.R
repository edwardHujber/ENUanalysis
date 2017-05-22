library("abind")
library("ggplot2")
library("RColorBrewer")

#### 'dat' should already have been defined by one of the import scripts. 
SEM<-function(x,na.rm=FALSE){
   if(na.rm){x<-x[!is.na(x)]}
   sd(x)/sqrt(length(x))}

barPlotwSEM<-function(table, ylab=NULL, xlab=NULL){
   if(length(dim(table))==0){meanCount<-table;SEMCount<-0}
   else{meanCount<- rowMeans(table)
   SEMCount<-apply(table,1,SEM)}
   print(c("means:",meanCount))
   print(c("SEMs:",SEMCount))
   bars<-barplot(meanCount, ylim=c(0,max(meanCount,na.rm = TRUE)*1.2),xlab=xlab, ylab=ylab)
   arrows(x0=bars,y0=meanCount-SEMCount,y1=meanCount+SEMCount, angle=90, code=3, length=0)
}

repFunctionOverStrains<-function(DAT,FUN){
   strains <- unique(DAT$ID)
   print(strains)
   nStrains <-length(strains)
   outMatrix<-NULL
   for(i in 1:nStrains){
      temp<-DAT[which(DAT$ID==strains[i]),]
      print(nrow(temp))
      funOut<-FUN(temp)
      print(funOut)
      # print(funOut)
      outMatrix<-abind(outMatrix,funOut,along=(1+length(dim(funOut))))
   }
   return(outMatrix)
}


## count variants, transitions, transversion
codingOnly <- FALSE
oneLinePerVar <- TRUE

datList <- c("f1EMS","f1ENU","denovo")
genomeList <- c("w","w","h")

for(k in 1:length(datList)){
   


useThisDat<-get(datList[k])
GCbias<-switch(genomeList[k],
   w = wormGCbias,
   h = humanGC)
      


if(codingOnly){useThisDat <- useThisDat[which(nchar(as.character(useThisDat$Old_codon.New_codon))!=0),]}
if(oneLinePerVar){
   CPRCs <- levels(factor(useThisDat$CPRC))
   for(i in 1:length(CPRCs)){
      if(length(which(useThisDat$CPRC==CPRCs[i]))>1) { useThisDat<- useThisDat[-which(useThisDat$CPRC==CPRCs[i])[-1],] }
   }
   
}


GAper<-repFunctionOverStrains(useThisDat,function(x){

useThisDat <- x

varTypes <- data.frame("Ref"=NA,"Var"=NA,"Count"=NA)

for(i in ntnum){
  for(j in ntnum[-i]){    
    
    cnt<-length(which(useThisDat$wt_dna==nt[i] & useThisDat$mut_dna==nt[j]))
    varTypes <- rbind(varTypes,c(nt[i],nt[j],cnt))
  }
}

varTypes <- varTypes[-1,]
varTypes$Percent <- apply(varTypes,1,function(x) as.numeric(x[3])/sum(as.numeric(varTypes$Count))   )


transitions<- varTypes[which((is.na(match(varTypes$Ref,purines))==FALSE & 
                                is.na(match(varTypes$Var,purines))==FALSE) | 
                               (is.na(match(varTypes$Ref,pyrimidines))==FALSE & 
                                  is.na(match(varTypes$Var,pyrimidines))==FALSE)),]

transversions <- varTypes[which((is.na(match(varTypes$Ref,purines))==FALSE & 
                                   is.na(match(varTypes$Var,pyrimidines))==FALSE) | 
                                  (is.na(match(varTypes$Ref,pyrimidines))==FALSE & 
                                     is.na(match(varTypes$Var,purines))==FALSE)),]




## Move varTypes to the cleaner mutTypes
mutType <- data.frame("mutType"=mutRows[c(6,4,5,1,2,3)])
mutType$Count <- apply(mutType,1,function(x)  
  sum(as.numeric(varTypes$Count[which(  (varTypes$Ref==substring(x,1,1)&varTypes$Var==substring(x,7,7)) | (varTypes$Ref==substring(x,3,3)&varTypes$Var==substring(x,9,9))   )]))
)
mutType$countGCnormal <- apply(mutType,1,function(x){
   if(grepl("A/T >",x['mutType'])){return(as.numeric(x['Count'])/(1-GCbias))}
   if(grepl("C/G >",x['mutType'])){return(as.numeric(x['Count'])/(GCbias))}
})
mutType$Percent <- apply(mutType,1,function(x) as.numeric(x['Count'])/sum(as.numeric(mutType$Count))   )
mutType$PercentGCnormal <- apply(mutType,1,function(x) as.numeric(x['countGCnormal'])/sum(as.numeric(mutType$countGCnormal))   )
 return(mutType)

})



mutSprectrum <- rowMeans(matrix(nrow=dim(GAper[,5,])[1],ncol=dim(GAper[,5,])[2],as.numeric(GAper[,5,])))


if (datList[k] == "denovo"){
   denovoSNPpercents
   RMs <- rbind(RMs,data.frame("Mutagen"=datList[k],"Mutation"=names(denovoSNPpercents), "Percent"=denovoSNPpercents))
   
}



if(k == 1){
   RMs <- data.frame("Mutagen"=datList[k],"Mutation"=GAper[,1,1], "Percent"=mutSprectrum)
}else{
   RMs <- rbind(RMs,data.frame("Mutagen"=datList[k],"Mutation"=GAper[,1,1], "Percent"=mutSprectrum))
}
}

ggplot(RMs, aes(Mutagen, Percent, fill=Mutation)) + 
   geom_bar(position="stack",stat="identity", width=0.5) +
   scale_fill_brewer(palette="Set3")


rm(list=setdiff(ls(), c("mutType",VARS2KEEP))) 
gc()

## Figure 1 prototype ##################
# barplot(as.matrix(mutType$PercentGCnormal))#####
########################################