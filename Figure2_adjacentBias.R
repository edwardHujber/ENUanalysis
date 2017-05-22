# library("Hmisc") ## required for inc() 
library("seqinr") ## required for read.fasta()
library("R.utils") ## required for extract()
library("abind") ## required for abind()

contextFolder <- "contexts/"

#### 'dat' should already have been defined by one of the import scripts. 


#######
getLocalNT <- function(rome, pos, back, forward){
  pos <- as.numeric(pos)
  if(rome=="MtDNA"){rome<-"M"}
  latin <- as.numeric(as.roman(rome))
  if(latin==10){latin<-6}
  if(latin==1000){latin<-7}
  local <- WS220[[latin]][(pos-back):(pos+forward)]
  local[back+1] <- "*"
  local <- paste0(local, collapse="")
  local <- toupper(local)
  return(local)
}

adjacent <- function(datRow, packed){
  chr<-datRow["chr"]
  pos<-as.numeric(datRow["pos"])
  ref<-datRow["wt_dna"]
  var<-datRow["mut_dna"]
  three <- packed[["three"]]
  five <- packed[["five"]]
  plus <- length(dim(five))-1   ## rederive "adj"
  localBP <- getLocalNT(chr,pos, plus,plus) ## get the local region. ref +/- plus bps
  if(ref == "C" | ref == "T"){localBP <- revcomp(localBP)} ## revcomp if we arent looking at G or A
  mut <- paste0(ref,"/",basePair(ref)," > ",var,"/",basePair(var), collapse="") ## build mut to directly compare with mutrows
  coord5 <-coord3<- which(mutRows == mut | mutRows == upsidedown(mut)) ## begin building coordinates to increment
  refpos <- nchar(localBP)-plus
  for(i in 1:plus){
    coord5 <- paste0(coord5,",",BP2num(substr(localBP,refpos-i,refpos-i)))
    coord3 <- paste0(coord3,",",BP2num(substr(localBP,refpos+i,refpos+i)))
  }
  eval(parse(text=paste0("five[",coord5,"]<-",eval(parse(text=paste0("five[",coord5,"]")))+1)))    ## this is gross. But increments at the defined coords
  eval(parse(text=paste0("three[",coord3,"]<-",eval(parse(text=paste0("three[",coord3,"]")))+1)))
  ret <- list("five"=five, "three"=three)
  return(ret)
}

revcomp <- function(oligo){
  nbase <- nchar(oligo)
  revcmp <- NULL
  while(nbase>0){
    revcmp <- paste0(revcmp,basePair(substring(oligo,nbase,nbase)) )
    nbase <- nbase-1
  }
  return(revcmp)
}

turn <- function(bp){
  paste0( substr(bp,7,7),"/",
          substr(bp,5,5),"|",
          substr(bp,3,3),"/",
          substr(bp,1,1)
  )
}


BP2num <- function(p){
  if(p == "A"){return(1)}
  if(p == "T"){return(2)}
  if(p == "C"){return(3)}
  if(p == "G"){return(4)}
}

local2BP <- function(loc){
  split <- unlist(strsplit(loc,""))
  pair<-unlist(lapply(split,basePair))
  mat<-cbind(split,pair)
  out <- paste0(mat[1,],collapse="/")
  for(i in 2:nrow(mat)){
    out<-paste0(c(out,paste0(mat[i,],collapse="/")),collapse="|")
  }
  return(out)
}

lookAtAdj <- function(mat, dims){
  apply(mat, c(1,(dims+1)), sum)
}

discardDups <- function(str){
  q <- 1
  while (q <= length(str)){
    str <- str[which(upsidedown(str) != str[q])]
    q <- q+1
  }
  return(str)
}

upsidedown <- function(bp){
  # "a/b > c/d" becomes "b/a > d/c"
  paste0(  substr(bp,3,3),"/",
           substr(bp,1,1)," > ",
           substr(bp,9,9),"/",
           substr(bp,7,7)
  )
}

countContextWithRC <- function(ctVec, gen){
  ctVec_RC <- sapply(ctVec,revcomp)
  print("Counting forward")
  countsF <- countContext(ctVec,gen)
  print("Counting reverse")
  countsR <- countContext(ctVec_RC,gen)
  totalCounts<-countsF+countsR
  return(totalCounts)
}

countContext <- function(ctVec, gen){   
    sapply(ctVec,function(SEARCHSEQ){
      print(paste0("Context ",which(ctVec==SEARCHSEQ), " of ",length(ctVec)))
      Reduce("+",lapply(gen, function(genome){
        print(getName(genome))
        length(gregexpr(SEARCHSEQ,getSequence(genome,as.string=TRUE)[[1]],ignore.case = TRUE)[[1]])
      }))
    })
}

adjacentBias <- function(varDat,WS220,pos5,pos3) {
   pos5<-pos5[order(pos5)]
   pos3<-pos3[order(pos3)]
   max3Distance <- max(pos3)
   max5Distance <- max(pos5)
   maxDistance <- max(pos3,pos5)
   
   
   contextName <- paste0("ctx_",paste(rev(pos5),collapse="n"),"x",paste(pos3,collapse="n"),".txt")
   if(file.exists(paste0(workPath,contextFolder,contextName))){
      print(paste0("Found a saved ",contextName," file, will use that one..."))
      searchGrid<-read.table(paste0(workPath,contextFolder,contextName))
   }else{
      print(paste0("No saved context file found, must generate ",contextName," now..."))
      
      contextSize <- 1+max3Distance+max5Distance
      contextList <- as.list(rep(paste0("[",paste(nt,collapse=","),"]"),contextSize)) ## rep [A,G,T,C] (ie, N)
      contextList[max5Distance+1+pos3] <- list(nt) ## nt list to 3' positions
      contextList[max5Distance+1-pos5] <- list(nt) ## nt list to 5' positions
      contextList[max5Distance+1] <- "[*]"         ## * to 0 postion
      contextGrid <-expand.grid(contextList)
      contextVector <- apply(contextGrid,1,function(x) paste0(x, collapse=""))
      
      searchList <- list("wt_dna"=unique(substr(mutRows,1,1)), "mut_dna"=nt, "context"=contextVector)
      searchGrid <- expand.grid(searchList)
      searchGrid <- cbind(searchGrid,"count"=NA, "contextCount"=NA, "Pmut"=NA ,"Pctx"=NA,"P"=NA, "expect"=NA, "pctDiff"=NA)
      searchGrid <- searchGrid[order(searchGrid[,1],searchGrid[,2],searchGrid[,3]),]
      searchGrid$wt_dna <- factor(searchGrid$wt_dna, levels=levels(searchGrid$mut_dna))
      dropRows <- which(searchGrid$wt_dna == searchGrid$mut_dna)
      searchGrid<-searchGrid[-dropRows,]
      
      ## expected will need to recieve adjustment for biases. Easiset way will probably be to build a context vector straight from WS220
      
      contextVectorWith0 <- c(apply(  as.array( unique(substr(mutRows,1,1))) ,1,function(x){ sub("[*]",x, contextVector)}))
      contextCounts <- countContextWithRC(contextVectorWith0,WS220) ## This is the big call to countContext
      
      ## puts context counts into searchGrid
      for (i in 1:length(contextCounts)){
         focusRow<-which(apply(searchGrid,1,function(x){sub("[*]",x['wt_dna'],x['context'] )}) == names(contextCounts[i]))
         searchGrid$contextCount[focusRow] <- contextCounts[i]
      }
      print(paste0("Saving to ",paste0(workPath,contextFolder,contextName)))
      write.table(searchGrid,paste0(workPath,contextFolder,contextName))
   }
   
   print("Retrieving variant contexts...")
   DNAWithContext <-  data.frame("wt_dna"=varDat$wt_dna, "mut_dna"=varDat$mut_dna,"context"=apply(varDat, 1, function(x) getLocalNT(x['chr'],x['pos'],back=maxDistance,forward=maxDistance)))
   
   ## counts the variants that fit each context
   print("Counting variant contexts...")
   searchGrid$count<- apply(searchGrid,1,function(x){
      sum(
         (DNAWithContext$wt_dna == x['wt_dna'] &  DNAWithContext$mut_dna == x['mut_dna'] & grepl(x['context'], DNAWithContext$context)) |
            ( DNAWithContext$wt_dna == basePair(x['wt_dna']) &  DNAWithContext$mut_dna == basePair(x['mut_dna']) & grepl(revcomp(x['context']), DNAWithContext$context) )
      ) 
   })
   
   ## how to group?? (or do I need to?)
   groups <- split(searchGrid,list(searchGrid$wt_dna,searchGrid$mut_dna),drop=TRUE) ## Mutation prefers a context?
   groups <- split(searchGrid,list(searchGrid$wt_dna,searchGrid$context),drop=TRUE) ## Context influences mutation?n
   
   ## calc Pcontext,expect, other previously used numbers
   for(j in 1:nrow(searchGrid)){
      Pmut <- ( sum(searchGrid$count[which(searchGrid$wt_dna == searchGrid$wt_dna[j] & searchGrid$mut_dna == searchGrid$mut_dna[j] )])) / sum(searchGrid$count[which(searchGrid$wt_dna== searchGrid$wt_dna[j])])
      
      Pctx <- ( sum(searchGrid$contextCount[which(searchGrid$wt_dna == searchGrid$wt_dna[j] & searchGrid$context == searchGrid$context[j] )]))/ sum(searchGrid$contextCount[ which(searchGrid$wt_dna == searchGrid$wt_dna[j])            ])
      
      searchGrid$Pmut[j]<-Pmut
      searchGrid$Pctx[j]<-Pctx
      searchGrid$P[j]<-Pmut*Pctx
      
      ## expect = total count of *this*Wt->mut * proportion found in a context.
      searchGrid$expect[j] <- Pctx * ( sum(searchGrid$count[         which(searchGrid$wt_dna == searchGrid$wt_dna[j] & searchGrid$mut_dna == searchGrid$mut_dna[j] )         ]))  
   }
   
   ## dividing P by 2 right now becasue its summing to 1 from Pcontext
 #  searchGrid$expect <-(searchGrid$P/2)*sum(searchGrid$count)  ## old expect formula. I think this is wrong. newer version in j loop immediatly above

   searchGrid$pctDiff <- 100*(searchGrid$count-searchGrid$expect)/searchGrid$expect
   
   return(searchGrid)
   
}

pos3 <- c(0)           ## define which nucelotide positions to keep. <= adj
pos5 <- c(1)

biasGrid <- adjacentBias(dat,WS220,pos5,pos3)

changes <- unique(paste0(biasGrid$wt_dna,">",biasGrid$mut_dna))
changeContext <- unique(biasGrid$context)

biasGridNicer <- data.frame(matrix(ncol=length(changeContext),nrow=length(changes)))
rownames(biasGridNicer)<-changes
colnames(biasGridNicer)<-changeContext

for( i in 1:nrow(biasGrid)){
   pct <-biasGrid$pctDiff[i]
   r <- paste0(biasGrid$wt_dna[i],">",biasGrid$mut_dna[i])
   c <- as.character(biasGrid$context[i])
   biasGridNicer[r,c]<-pct
}

#### this for for if you want to test a whole series of contexts. set it up and go to lunch
# meanPctDiff <- NA
# for(i in 1:10){
#   pos3 <- as.numeric(i)
#   pos5 <- as.numeric(0)
#   biasGrid <- adjacentBias(WS220,pos3,pos5)
#   meanPctDiff[i]<-mean(abs(biasGrid$pctDiff))
# }

a<-list(

   list(c(1,2),c(0)),
   list(c(0),c(1,2)),
   list(c(1,2,3),c(0)),
   list(c(0),c(1,2,3)),
   list(c(1,2,3),c(1,2,3))
)
b<-lapply(a,function(x){
   adjacentBias(WS220,x[[1]],x[[2]])
   }
   )


RMSfive5 <- mean(searchGrid$pctDiff)

# system.time(contextMatrix3<-buildContextMatrix(WS220,pos3,FALSE))
# system.time(contextMatrix5<-buildContextMatrix(WS220,pos5,TRUE))


contextFraction3 <- contextMatrix3/sum(contextMatrix3)


expected3 <- sum(threePrimeFocused)*contextFraction3

pctdiff5 <- (expected5-fivePrimeFocused)/expected5
pctdiff3 <- 100*(expected3-threePrimeFocused)/expected3


## maybe this shuould be a chi-squared metric? I would like a more intuitive "% bias" metic though.
RMSpctdiff5 <- sqrt(sum(pctdiff5^2)/prod(dim(threePrimeFocused)))
RMSpctdiff3 <- sqrt(sum(pctdiff3^2)/prod(dim(fivePrimeFocused)))




