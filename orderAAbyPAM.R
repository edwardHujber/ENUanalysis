calcScore<-function(order,scoreMatirx,choseNumber=1){
   orderPairs<-cbind(order[-length(order)],order[-1])
   scores<-apply(orderPairs,1,function(x){
      return(scoreMatirx[x[1],x[2]])
   })
   return(c(sum(scores),which(scores==scores[order(scores,decreasing =LowScore)][choseNumber])))
}
cleanMatix<-function(x){ # clean up: remove ambiguous and stop, then sort
   x<-x[,-which(is.na(match(colnames(x),allAA$AA)))]
   x<-x[-which(is.na(match(rownames(x),allAA$AA))),]
   x<-x[,order(colnames(x))]
   x<-x[order(rownames(x)),]
   return(x)
}

# http://www.ncbi.nlm.nih.gov/Class/FieldGuide/BLOSUM62.txt
BLOSUM62 <- read.table(paste0(workPath,"BLOSUM62.txt"), header=TRUE, quote="\"", stringsAsFactors=FALSE)
BLOSUM62_cln <- cleanMatix(BLOSUM62)

# order by PAM?
# https://raw.githubusercontent.com/noporpoise/seq-align/master/scoring/PAM250.txt
PAM250 <- read.table(paste0(workPath,"PAM250.txt"), header=TRUE, quote="\"", stringsAsFactors=FALSE)
PAM250_cln <- cleanMatix(PAM250)


# PAM250nosyn<-PAM250
# diag(PAM250nosyn)<-NA
# invOdds<-1/(2^PAM250nosyn)


scoringMatrix<-BLOSUM62_cln
LowScore <- FALSE

order<-as.character(allAA$AA)
Best<-100*LowScore
k<- -1

while(TRUE){
   k<-k+1
   cat('\r',format(Sys.time(),'%H:%M:%S'), k*100, " tries" )
   gc()
   for(j in 1:100){
      order<- sample(order)
      checkEdgeNum<-1
      for(i in 1:10000){
         Scores <- calcScore(order,scoringMatrix,checkEdgeNum)
          # print(Scores[1])
         scorePos<-Scores[sample(length(Scores)-1,1)+1]
         swap<-scorePos+sample(0:1,1)
         newOrders<-sapply(1:20,function(x){
            replace(order, c(x, swap), order[c(swap, x)])
         })
         newScores<-lapply(split(t(newOrders),seq(NROW(t(newOrders)))),calcScore,scoringMatrix)
         newLowIdx<-which(sapply(newScores,function(x){x[1]})[order(sapply(newScores,function(x){x[1]}),decreasing =!LowScore)[1]]==sapply(newScores,function(x){x[1]}))
         choseNewLow <- newLowIdx[sample(length(newLowIdx),1)]
         # print(newScores[[choseNewLow]][1])
         if(choseNewLow==swap){checkEdgeNum<-checkEdgeNum+1}else{checkEdgeNum <- 1}
         if(checkEdgeNum>19){break}
         order<-newOrders[,choseNewLow]
      }
      if(calcScore(order,scoringMatrix)[1]>=Best){
         Best<-calcScore(order,scoringMatrix)[1]
         cat("\n",(order),calcScore(order,scoringMatrix)[1], "\n")
      }
   }
}