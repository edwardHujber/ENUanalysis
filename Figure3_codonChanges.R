library("plyr") ## required for ddply()
library("seqinr")
library("abind")
library("gplots")
####### functions ########

SEM<-function(x,na.rm=FALSE){
   if(na.rm){x<-x[!is.na(x)]}
   sd(x)/sqrt(length(x))}

barPlotwSEM<-function(table, ylab=NULL){
   if(length(dim(table))==0){meanCount<-table;SEMCount<-0}
   else{meanCount<- rowMeans(table)
   SEMCount<-apply(table,1,SEM)}
   print(c("means:",meanCount))
   print(c("SEMs:",SEMCount))
   bars<-barplot(meanCount, ylim=c(0,max(meanCount)*1.2),xlab="Position in codon", ylab=ylab)
   arrows(x0=bars,y0=meanCount-SEMCount,y1=meanCount+SEMCount, angle=90, code=3, length=0)
}

repFunctionOverStrains<-function(DAT,FUN){
   strains <- unique(DAT$ID)
   nStrains <-length(strains)
   outMatrix<-NULL
   for(i in 1:nStrains){
      temp<-DAT[which(DAT$ID==strains[i]),]
      funOut<-FUN(temp)
      # print(funOut)
      outMatrix<-abind(outMatrix,funOut,along=(1+length(dim(funOut))))
   }
   return(outMatrix)
}

####### codon position #########
##### NOT CURRENTLY NORMALIZED FOR CODON FREQUENCY ########
synonymousOnly=FALSE
codonPosTable<-repFunctionOverStrains(dat,function(x){
   if(synonymousOnly){x<-x[which(substr(x$old_AA.new_AA,1,1)==substr(x$old_AA.new_AA,3,3) & c(x$old_AA.new_AA)!=which(nchar(levels(x$old_AA.new_AA))==0)),]}
   if(nrow(x)>0){
      positions <-sapply(as.character(x$Old_codon.New_codon),function(aa){
         if(!grepl("/",aa) | (nchar(aa)!=7)){0} 
         else {
            which(sapply(1:3,function(y){
               substr(aa,y,y)==toupper(substr(aa,y,y))  ## which position in codon "aa" is capital (ie the mut)
            }))
         }
      })
      tab<-table(factor(positions, levels=0:3)) ## do it this way to include 0s
   }
})


barPlotwSEM(codonPosTable[-1,],ylab="Counts")

codonPosTable<-cbind(0,codonPosTable) ## This juggle necessary for if there's only 1 line, like in CB4856. Removed 2 lines below
codonPosTablePct<-100*codonPosTable[-1,]/rbind(apply(codonPosTable[-1,],2,sum),apply(codonPosTable[-1,],2,sum),apply(codonPosTable[-1,],2,sum))
codonPosTablePct<-codonPosTablePct[,-1]
barPlotwSEM(codonPosTablePct,ylab="Percent")


######## codon changes #########


OldNewCodon <-as.character(dat$Old_codon.New_codon[which(grepl("/",dat$Old_codon.New_codon))])
if(length(which(nchar(OldNewCodon)!=7)!=0)){OldNewCodon<-OldNewCodon[-which(nchar(OldNewCodon)!=7)]}

CCnums<-matrix(0,dim(CodonChanges)[1],dim(CodonChanges)[2])
a<-sapply(strsplit(as.character(OldNewCodon),"/"),function(x){
   CCrow<-which(rownames(CodonChanges)==toupper(x[1]))
   CCcol<-which(CodonChanges[CCrow,]==toupper(x[2]))
   CCnums[CCrow,CCcol]<<-CCnums[CCrow,CCcol]+1
})
CCnums<-data.frame(CCnums)
rownames(CCnums)<-codons
colnames(CCnums)<-posChanges

CCnumsNormalized <-CCnums/codonUsage

Per<-1000
CCPer <- CCnums/(sum(CCnums,na.rm=TRUE)/Per)
heatmap.2( as.matrix(CCPer),
           col = colorRampPalette( c("blue","red","yellow","white"), space="rgb")(100),
           cellnote = CodonChanges,
           notecol="black",
           trace = "none", 
           na.color="gray",
           Rowv = FALSE,
           Colv = FALSE,
           dendrogram = "none",
           offsetRow = -46,
           offsetCol = -55,
           srtCol = 0,
           main=paste0(datName,": # mutations per ",Per)
       
)



CCnumsNormalizedPer <- CCnumsNormalized/(sum(CCnumsNormalized,na.rm=TRUE)/Per)
# F1ENUWRONG <- CCnumsNormalizedPer
# F1ENUWRONG["TAG","P2_I"]<-0
# CCnumsNormalizedPer<-F1ENUWRONG
assign(paste0("CCnumsNormalizedPer_",datName),CCnumsNormalizedPer)


#  CCnumsNormalizedPer<-CCnumsNormalizedPer_f1ENU
#  CCnumsNormalizedPer<-CCnumsNormalizedPer_f1EMS

heatmap.2( as.matrix(CCnumsNormalizedPer),
           col = colorRampPalette( c("blue","red","yellow","white"), space="rgb")(100),
           cellnote = CodonChanges,
           notecol="black",
           trace = "none", 
           na.color="gray",
           Rowv = FALSE,
           Colv = FALSE,
           dendrogram = "none",
           offsetRow = -45,
           offsetCol = -55,
           srtCol = 0,
           main=paste0(datName,": # mutations per ",Per, ". Normalized by codon freq."),
           breaks=seq(0,25,length=101)
)

####### AA changes #####
AAchanges <-data.frame(matrix("",length(allAAwStop),length(allAAwStop)),row.names = allAAwStop,stringsAsFactors = FALSE)

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
missenseChanges<-do.call(rbind,AAchanges[1:20,1:20])
colnames(missenseChanges)<-rownames(missenseChanges)
class(missenseChanges)<-"numeric"
diag(missenseChanges)<-NA
## counts each change
AAchangeTable<-repFunctionOverStrains(dat,function(x){
   changesObserved <- x$old_AA.new_AA[which(nchar(as.character(x$old_AA.new_AA))>0)]
   thisAAchange <- AAchanges
   for(i in 1:length(changesObserved)){
      thisAAchange[substr(changesObserved[i],1,1), substr(changesObserved[i],3,3)] <- 
         as.numeric(thisAAchange[substr(changesObserved[i],1,1), substr(changesObserved[i],3,3)])+1
   }
   return(thisAAchange)
})
class(AAchangeTable)<-"numeric"

## summary sheets
AAchangeSum <- apply(AAchangeTable,1:2,sum)

assign(paste0("AAchangePerc_",datName),(AAchangeSum/sum(AAchangeSum,na.rm = TRUE)))

AAchangeMean <- apply(AAchangeTable,1:2,mean)
AAchangeSEM <- apply(AAchangeTable,1:2,SEM)

synPerAll<-apply(AAchangeTable,3,function(x){sum(diag(x),na.rm = TRUE)/sum(x,na.rm = TRUE)})
print(paste0(mean(synPerAll)*100," ± ",SEM(synPerAll)*100," % synonymous"))

## heatmap
Per<-1000
orderByBLO <- TRUE
if(orderByBLO){AAchangeSum<-AAchangeSum[AAorderedbyBLO,AAorderedbyBLO]}
changePer <- AAchangeSum/(sum(AAchangeSum,na.rm=TRUE)/Per)

#  changePer<-AAchangePer_mmpEMS
#  changePer<-AAchangePer_mmpENU



heatmap.2( changePer,
           col = colorpanel(100,"#241E1E","#FA758A","#FFFB5D"),
           cellnote = round(changePer, digits=2),
           notecol="black",
           trace = "none", 
           na.color="gray20",
           Rowv = FALSE,
           Colv = FALSE,
           dendrogram = "none",
           offsetRow = -46,
           offsetCol = -45,
           srtCol = 0,
           main=paste0(datName,": # mutations per ",Per),
           breaks=seq(0,40,length=101)
)

assign(paste0("AAchangePer_",datName),changePer)

###### comparison heatmap ####

## regorup into classes?
## Reorder by BLOSUM?


groupByClass <- FALSE


numerator   <-   data.matrix(CCnumsNormalizedPer_f1EMS)
demoninator <-   data.matrix(CCnumsNormalizedPer_f1ENU)

# numerator   <-AAchangePer_mmpEMS
# demoninator <-AAchangePer_mmpENU

   
diff<-diff1<-log(numerator/demoninator,2)
if(groupByClass){diff<-diff1<-diff[allAA$AA[order(allAA$class)],allAA$AA[order(allAA$class)]]}
if(orderByBLO){diff<-diff1<-diff[AAorderedbyBLO,AAorderedbyBLO]}

CPrange<-100
CP<-colorpanel(CPrange,"#c51b7d","gray","#4d9221")
rangeMax <- max(diff[which(is.finite(diff))],na.rm=TRUE)
rangeMin <- min(diff[which(is.finite(diff))],na.rm=TRUE)
range <-rangeMax-rangeMin
rangeABS <- max(abs(c(rangeMax,rangeMin)))

if(sum((diff<0)&(is.infinite(diff)))!=0){
   CP<-c("#c51b7d",CP)
   diff[(diff<0)&(is.infinite(diff))]<- -rangeABS-(range/CPrange)
      }
if(sum((diff>0)&(is.infinite(diff)))!=0){
   CP<-c(CP,"#4d9221")
   diff[(diff>0)&(is.infinite(diff))]<- rangeABS+(range/CPrange)
   }

heatmap.2( diff,
           col = CP,
           cellnote = round(diff1, digits=2),
           notecol="black",
           trace = "none", 
           symbreaks = min(diff, na.rm=TRUE),
           na.color="gray55",
           Rowv = FALSE,
           Colv = FALSE,
           dendrogram = "none",
           offsetRow = -45,
           offsetCol = -53,
           srtCol = 0,
           main="MMP: log2 ratio of EMS to ENU mutational frequency"
)
barplot(
   table(
      factor(
         round(diff1),
         levels=
            c(-Inf,
              min(round(diff),na.rm = T) : max(round(diff),na.rm = T),
              Inf)
         )),
   col=
      colorpanel(
         3+max(round(diff),na.rm = T)-min(round(diff),na.rm = T),
         "#c51b7d","gray","#4d9221")
)
        
####### class changes (missense only) #####
classChanges <- data.frame(matrix("",length(classes),length(classes)),row.names = classes,stringsAsFactors = FALSE)
colnames(classChanges) <- classes

AAchangeTableNonsyn<-array(apply(AAchangeTable,3,function(x){  ## Put zeros along the synonymous diagonals
   diag(x)<-NA
   return((x))
}),dim=dim(AAchangeTable))

classChangeList<-apply(AAchangeTableNonsyn,3,function(x){
   CCTemp<-classChanges
   for(i in 1:length(classes)){
      for(j in 1:length(classes)){
         CCTemp[i,j] <- sum(x[which(allAA$class==classes[i]),which(allAA$class==classes[j])],na.rm=TRUE)
      }
   }
   return(CCTemp)
})
classChangeArray<-do.call(abind,c(classChangeList,along=3))
class(classChangeArray)<-"numeric"

## summary sheets
CCSum <- apply(classChangeArray,1:2,sum)
CCMean <- apply(classChangeArray,1:2,mean)
CCSEM <- apply(classChangeArray,1:2,SEM)

inClassPerAll<-apply(classChangeArray,3,function(x){sum(diag(x),na.rm = TRUE)/sum(x,na.rm = TRUE)})
print(paste0(mean(inClassPerAll,na.rm=TRUE)*100," ± ",SEM(inClassPerAll,na.rm=TRUE)*100," % stayed in class"))

####### PAM score #########




### Including synonymous
# from the change counts, shave off the column,row corrosponding to Stop
AAchangesAAonly <- AAchangeTable[,-which(is.na(match(colnames(AAchanges),allAA$AA))),]
AAchangesAAonly <- AAchangesAAonly[-which(is.na(match(colnames(AAchanges),allAA$AA))),,]

# multiply the AA change counts by the PAM, then get the mean
PAMscoresWsyn <-apply(AAchangesAAonly,3,function(x){
   sum((x*PAM250),na.rm = TRUE)/sum(x,na.rm = TRUE)
})
print(paste0("PAM score with syn: ",mean(PAMscoresWsyn)," ± ",SEM(PAMscoresWsyn)))

### Missense only
# from the change counts, shave off the column,row corrosponding to Stop
AAchangesAAonlyNonsyn <- AAchangeTableNonsyn[,-which(is.na(match(colnames(AAchanges),allAA$AA))),]
AAchangesAAonlyNonsyn <- AAchangesAAonlyNonsyn[-which(is.na(match(colnames(AAchanges),allAA$AA))),,]

# multiply the AA change counts by the PAM, then get the mean
PAMscoresWoutsyn <-apply(AAchangesAAonlyNonsyn,3,function(x){
   sum((x*PAM250),na.rm = TRUE)/sum(x,na.rm = TRUE)
})
print(paste0("PAM score without syn: ",mean(PAMscoresWoutsyn,na.rm=TRUE)," ± ",SEM(PAMscoresWoutsyn,na.rm=TRUE)))

## Make a histogram of scores.
AAdat <- sample(dat$old_AA.new_AA[which(  ## select only AA changes. sample will randomize the order, making the exploration curve (next section) smoother
   nchar(as.character(dat$old_AA.new_AA))>0 &
      substr(dat$old_AA.new_AA,1,1)!="*" &
      substr(dat$old_AA.new_AA,3,3)!="*" )]) 
   
PAMhist<-table(do.call(c,lapply(AAdat,function(x){
   PAM250[substr(x,1,1),substr(x,3,3)]
})))


PAMhist_missense<-table(do.call(c,lapply(AAdat[which(substr(AAdat,1,1)!=substr(AAdat,3,3))],function(x){
   PAM250[substr(x,1,1),substr(x,3,3)]
})))

assign(paste0("PAMhist_",datName),PAMhist)
assign(paste0("PAMhistMis_",datName),PAMhist_missense)


histTable_mis <- rbind((PAMhistMis_mmpEMS/max(PAMhistMis_mmpEMS)),(PAMhistMis_mmpENU/max(PAMhistMis_mmpENU)))
barplot(histTable_mis, main="title",
        xlab="PAM value", col=c("darkblue","red"),legend = c("mmpEMS","mmpENU"),
      beside=TRUE)

# 

####### Codon Table exploration ######
iter <- 10
AAdat<-AAdat[sample(1:length(AAdat))] ## randomizes the array again
AAdat<-AAdat[1:1000]

## each var sits in its own change table. synonymous vars are untouched tables
indyAAchange<-lapply(AAdat,function(x){
   thisAAchange <- missenseChanges
   if(substr(x,1,1)!=substr(x,3,3)){thisAAchange[substr(x,1,1), substr(x,3,3)] <- 1}
   return(thisAAchange)
})
indyAAchangeArray<-do.call(abind,c(indyAAchange,along=3))
rm(indyAAchange)
gc()
class(indyAAchangeArray)<-"numeric"

pb<-txtProgressBar(min=0,max=(dim(indyAAchangeArray)[3]),style=3)
itered<-NULL
for(j in 1:iter){
   indyAAchangeArray<-indyAAchangeArray[,,sample(1:dim(indyAAchangeArray)[3])] ## randomizes the array again
   print(paste0("Iter ",j))
   
   
   tableExploration<-NULL
   composite <- as.matrix(missenseChanges)
   class(composite)<-"numeric"
   for(i in 1:(dim(indyAAchangeArray)[3])){
      setTxtProgressBar(pb, i)
      composite <- composite+indyAAchangeArray[,,i]
      explored<-sum(apply(composite,1:2,sum)!=0,na.rm=TRUE)  /   sum(!is.na(apply(composite,1:2,sum)),na.rm=TRUE)
      tableExploration<-c(tableExploration,explored)
   }
   itered <- cbind(itered,tableExploration)
}
close(pb)

iterMean<-apply(itered,1,mean)

assign(paste0("tableExplore_",datName),itered)
assign(paste0("tableExploreMean_",datName),iterMean)

plot(1, type="n", ylim=c(0,1), xlim=c(0,20000), xlab="", ylab="")

# apply(tableExplore_f1EMS,2,lines,col="gray80")
# points(tableExploreMean_f1EMS,pch=20,cex=0.75,col="black")


apply(tableExplore_mmpEMS,2,lines,col="lightgreen")
#apply(tableExplore_mmpENU,2,lines,col="lightpink")
 apply(tableExplore_f1EMS,2,lines,col="cadetblue1")
 apply(tableExplore_f1ENU,2,lines,col="gray80")


points(tableExploreMean_mmpEMS,pch=20,cex=0.75,col="green")
#points(tableExploreMean_mmpENU,pch=20,cex=0.75,col="red")
 points(tableExploreMean_f1EMS,pch=20,cex=0.75,col="blue")
 points(tableExploreMean_f1ENU,pch=20,cex=0.75,col="black")




plot(tableExplore_mmpEMS,
    pch=20,
    cex=0.75,
    ylab="cumulative % \n possible AA changes observed",
    xlab="# Variants",
    mgp=c(2.2,.8,0),
    las=1,
    xlim=c(0,20000)
    )
points(tableExplore_mmpENU,col="red", pch=20, cex=0.75)
points(tableExplore_labellaENU,col="blue", pch=20, cex=0.75)

