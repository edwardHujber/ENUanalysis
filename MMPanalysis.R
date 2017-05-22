library("seqinr", lib.loc="~/R/win-library/3.0") ## required for read.fasta()
library("Hmisc", lib.loc="~/R/win-library/3.0") ## required for inc()
library("plyr", lib.loc="~/R/win-library/3.0") ## required for ddply()

cutoff <- 2

getLocalNT <- function(chr, pos, back, forward){
  if(chr=="MtDNA"){chr<-"M"}
  latin <- as.numeric(as.roman(chr))
  if(latin==10){latin<-6}
  if(latin==1000){latin<-7}
  local <- WS220[[latin]][(pos-back):(pos+forward)]
  local[back+1] <- "*"
  local <- paste0(local, collapse="")
  local <- toupper(local)
  return(local)
}

searchDinucs <- function(chr, pos, ref, var, mat){
  priorNT <- getLocalNT(chr,pos, 1,0)
  col2incr <- which(dinucCols == priorNT)
  row2incr <- which(changeRows == paste0(ref,">",var, collapse=""))
  inc(mat[row2incr,col2incr]) <- 1
  
  priorNT <- getLocalNT(chr,pos, 0,1)
  col2incr <- which(dinucCols == priorNT)
  row2incr <- which(changeRows == paste0(ref,">",var, collapse=""))
  inc(mat[row2incr,col2incr]) <- 1
  
  return(mat)
}

searchTrinucs <- function(chr, pos, ref, var, mat){
  priorNT <- getLocalNT(chr,pos, 2,0)
  col2incr <- which(trinucCols == priorNT)
  row2incr <- which(changeRows == paste0(ref,">",var, collapse=""))
  inc(mat[row2incr,col2incr]) <- 1
  
  priorNT <- getLocalNT(chr,pos, 1,1)
  col2incr <- which(trinucCols == priorNT)
  row2incr <- which(changeRows == paste0(ref,">",var, collapse=""))
  inc(mat[row2incr,col2incr]) <- 1
  
  priorNT <- getLocalNT(chr,pos, 0,2)
  col2incr <- which(trinucCols == priorNT)
  row2incr <- which(changeRows == paste0(ref,">",var, collapse=""))
  inc(mat[row2incr,col2incr]) <- 1
  
  return(mat)
}

## read in data
WS220 <- read.fasta("~/Lab/R/ENU analysis/WS220.fasta")

mmp <- read.delim("~/Lab/R/ENU analysis/mmp_mut_strains_data_Mar14.txt",  stringsAsFactors=FALSE)

mmp <- mmp[which((nchar(mmp$wt_dna)==1) & nchar(mmp$mut_dna)==1),]  ## drop indels

## discard nonunique variants
mmp$CPRC <- apply(mmp,1,function(x) paste(x[4],x[5],x[6],x[7]))
CPRC <- data.frame("CPRC" = mmp$CPRC )
CPRCcount <- data.frame(table(CPRC))
parental <- CPRCcount[which(CPRCcount$Freq>=cutoff),]      ## "parental" variants must show up >= cutoff number of times
uniqueVars <- mmp[which(is.na(match(mmp$CPRC,parental$CPRC))),]   ## remove variants that match a parental



## split by VC number
UV_TMP <- uniqueVars[grep("VC10", uniqueVars$strain),]
EMS <- uniqueVars[grep("VC20", uniqueVars$strain),]
ENU <- uniqueVars[grep("VC30", uniqueVars$strain),]
EMS_ENU <- uniqueVars[grep("VC40", uniqueVars$strain),]

mutagens <- c("ENU")

purines <- c("A","G")
pyrimidines <- c("T","C")
nt <- c(purines,pyrimidines)
ntnum <- 1:length(nt)

for (l in 1:length(mutagens)){
  vcf <- get(mutagens[l])

  
  ## count variants, transitions, transversion
  
  varTypes <- data.frame("Ref"=NA,"Var"=NA,"Count"=NA)
  
  for(i in ntnum){
    for(j in ntnum[-i]){    
      
      cnt<-length(which(vcf$wt_dna==nt[i] & vcf$mut_dna==nt[j]))
      varTypes <- rbind(varTypes,c(nt[i],nt[j],cnt))
    }
  }
  varTypes <- varTypes[-1,]
  
  transitions<- varTypes[which((is.na(match(varTypes$Ref,purines))==FALSE & 
                                  is.na(match(varTypes$Var,purines))==FALSE) | 
                                 (is.na(match(varTypes$Ref,pyrimidines))==FALSE & 
                                    is.na(match(varTypes$Var,pyrimidines))==FALSE)),]
  
  transversions <- varTypes[which((is.na(match(varTypes$Ref,purines))==FALSE & 
                                     is.na(match(varTypes$Var,pyrimidines))==FALSE) | 
                                    (is.na(match(varTypes$Ref,pyrimidines))==FALSE & 
                                       is.na(match(varTypes$Var,purines))==FALSE)),]
  
  
  
  ## define rows for dinuc and trinuc
  itions <- mapply(paste0, transitions$Ref, ">", transitions$Var)
  versions <- mapply(paste0, transversions$Ref, ">", transversions$Var)
  changeRows <- c(itions,versions)
  
  
  
  ## define columns, matrix for dinuc
  dinucCols <- c(mapply(paste0, nt, "*", collapse=""),mapply(paste0,"*",nt, collapse=""))
  dinucMat <- matrix(0,length(changeRows), length(dinucCols))
  
  ## define columns, matrix for trinuc
  c <- list(f=nt,f2=nt, m="*",e=nt,e2=nt)
  pentanuc <- expand.grid(c)
  pentanuc <- paste0(pentanuc$f,pentanuc$f2,pentanuc$m,pentanuc$e,pentanuc$e2)
  f3 <- unique(substr(pentanuc,1,3))
  m3 <- unique(substr(pentanuc,2,4))
  l3 <- unique(substr(pentanuc,3,5))
  trinucCols <- c(f3,m3,l3)
  trinucMat <- matrix(0,length(changeRows), length(trinucCols))
  
  ## search local
  for ( k in 1:nrow(vcf)){
    print(k)
    dinucMat <- searchDinucs(vcf$chr[k],vcf$pos[k],vcf$wt_dna[k],vcf$mut_dna[k],dinucMat)
    trinucMat<- searchTrinucs(vcf$chr[k],vcf$pos[k],vcf$wt_dna[k],vcf$mut_dna[k], trinucMat)
  }
  
  diNucs <- data.frame(dinucMat, row.names=changeRows)
  names(diNucs) <- dinucCols
  diNucs <- rbind(diNucs, total = lapply(diNucs, sum))
  diNucs <- cbind(diNucs, total = rowSums(diNucs))
  
  triNucs <- data.frame(trinucMat, row.names=changeRows)
  names(triNucs) <- trinucCols
  triNucs <- rbind(triNucs, total = lapply(triNucs, sum))
  triNucs <- cbind(triNucs, total = rowSums(triNucs))
  
  ## assign everything
  assign(paste0(mutagens[l],"transitions"),transitions)
  assign(paste0(mutagens[l],"transversions"),transversions)
  assign(paste0(mutagens[l],"diNucs"),diNucs)
  assign(paste0(mutagens[l],"triNucs"),triNucs)
  
}