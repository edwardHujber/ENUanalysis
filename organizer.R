library("seqinr")
###################  define basic stuff we'll need  #######################
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
basePair <- function(p){
   if(p == "A"){return("T")}
   if(p == "T"){return("A")}
   if(p == "C"){return("G")}
   if(p == "G"){return("C")}
   if(p == "*"){return("*")}
   if(p == "["){return("]")}
   if(p == "]"){return("[")}
   if(p == ","){return(",")}
}



workPath <- "~/Lab/R/ENU analysis/"
# workPath <- "Z:/programming/R/ENU analysis/"


CDS<-read.fasta(paste0(workPath,"Caenorhabditis_elegans.WBcel235.cds.all.fa"))
codonCounts<-colSums(do.call(rbind,lapply(CDS,uco)))
codonUsage<-codonCounts/sum(codonCounts)
codonUsage[order(codonUsage,decreasing = T)]


cutoff <- 2 ## varient is called parental is it shows up in 'cutoff' or more strains 
WS220 <- read.fasta(paste0(workPath,"WS220.fasta"))
wormGCbias<-0.354396697367839796921629158532596193253993988037109375 ## full precision of expression below.
# GCbias<- 
#   Reduce("+",lapply(WS220, function(genome){
#     print(getName(genome))
#     length(gregexpr("G|C",getSequence(genome,as.string=TRUE)[[1]],ignore.case = TRUE)[[1]])
#   })) / sum(sapply(WS220,length))  ## number of occurences of G or C divided by total nts

##### various AA stuff #####
AA_ASA <- read.csv(paste0(workPath,"/AA_ASA.csv"), comment.char="#")
AA_ASA$random <- rnorm(n = nrow(AA_ASA))

hydrophobic <- c("A","I","L","V")
aromatic <- c("F","W","Y")
polarNeutral <- c("N","C","Q","M","S","T")
acidic <- c("D","E")
basic <- c("R","H","K")
unique <- c("G","P")

allAAwClass <- data.frame("AA" = c(hydrophobic,aromatic,polarNeutral,acidic,basic,unique), 
                    "class"= c(rep("hydrophobic",length(hydrophobic)),
                               rep("aromatic",length(aromatic)),
                               rep("polarNeutral",length(polarNeutral)),
                               rep("acidic",length(acidic)),
                               rep("basic",length(basic)),
                               rep("unique",length(unique))))
allAAwClass<-allAAwClass[order(allAAwClass$AA),]
allAA <- allAAwClass
allAAwStop <- c(as.character(allAA$AA),"*")
AAorderedbyBLO <- c(unlist(strsplit("P K R Q E D N H Y F W M L I V T S G A C"," ")),"*")
# charged <- c(polarNeutral,acidic,basic)
# uncharged <- c(hydrophobic,aromatic)
classes <- factor(levels(allAAwClass$class))

# order by PAM?
# https://raw.githubusercontent.com/noporpoise/seq-align/master/scoring/PAM250.txt
PAM250 <- read.table(paste0(workPath,"PAM250.txt"), header=TRUE, quote="\"", stringsAsFactors=FALSE)
# clean up PAM250: remove ambiguous and stop, then sort
PAM250<-PAM250[,-which(is.na(match(colnames(PAM250),allAA$AA)))]
PAM250<-PAM250[-which(is.na(match(rownames(PAM250),allAA$AA))),]
PAM250<-PAM250[,order(colnames(PAM250))]
PAM250<-PAM250[order(rownames(PAM250)),]

####### varios NT stuff ######
purines <- c("A","G")
pyrimidines <- c("T","C")
nt <- sort(c(purines,pyrimidines))
ntnum <- 1:length(nt)
AT <- c("A/T","T/A")
GC <- c("G/C","C/G")
allPairs <- sort(c(AT,GC))
changeRows <- expand.grid(list(allPairs," > ", allPairs))
changeRows <- changeRows[ order(changeRows[,1]),]
dropRows <- which(changeRows[,1]==changeRows[,3])
mutRows <- discardDups(apply(changeRows[-dropRows,],1,function(x) paste0(x,collapse="")))



VARS2KEEP <- ls()


########################## IMPORT ############################
## must generate dataframe 'dat' with columns "ID" ,"chr","pos","wt_dna","mut_dna","Old_codon.New_codon","old_AA.new_AA"

## import MMP EMS. 
getDat<- "EMS"
VARS2KEEP <- ls()
source(paste0(workPath,"importData_MMP.R")) ## returns 'dat'
mmpEMS <- dat ## store it if you want to come back fast
dat <- mmpEMS; datName<-"mmpEMS" ## recall it to use it for figures

## import MMP ENU. 
getDat<- "ENU"
VARS2KEEP <- ls()
source(paste0(workPath,"importData_MMP.R")) ## returns 'dat'
mmpENU <- dat ## store it if you want to come back fast
dat <- mmpENU; datName<-"mmpENU" ## recall it to use it for figures

## import F1 EMS. 
getDat<- "EMS"
VARS2KEEP <- ls()
source(paste0(workPath,"importData_F1_mut.R")) ## returns 'dat'
f1EMS <- dat ## store it if you want to come back fast
dat <- f1EMS; datName<-"f1EMS" ## recall it to use it for figures

## import F1 ENU. 
getDat<- "ENU"
VARS2KEEP <- ls()
source(paste0(workPath,"importData_F1_mut.R")) ## returns 'dat'
f1ENU <- dat ## store it if you want to come back fast
dat <- f1ENU; datName<-"f1ENU" ## recall it to use it for figures

# ## import Hawaiian CB4856. path to data is set inside the import script.
# VARS2KEEP <- ls()
# source(paste0(workPath,"importData_CB4856.R")) ## returns 'dat'
# CB4856 <- dat ## store it if you want to come back fast
# dat <- CB4856; datName<-"CB4856"## recall it to use it for figures

# ## import Labella ENU. path to data is set inside the import script.
# VARS2KEEP <- ls()
# source(paste0(workPath,"importData_matt.R")) ## returns 'dat'
# Labella <- dat ## store it if you want to come back fast
# dat <- Labella; datName<-"labellaENU"## recall it to use it for figures

########################## FIGURES ############################

## Figure 1
## Figure 1 outputs GC normalized base-level mutational spectrum.
VARS2KEEP <- ls()
source(paste0(workPath,"Figure1_spectrumChart.R"))
mutType  ## numbers
barplot(as.matrix(mutType$PercentGCnormal)) ## as a chart


## Figure 2
VARS2KEEP <- ls()
source(paste0(workPath,"Figure2_adjacentBias.R"))


## Figure 3
VARS2KEEP <- ls()
source(paste0(workPath,"Figure3_codonChanges.R"))















