
doc = scan("~/Lab/R/ENU analysis/wormbase/ENUsearch.html",what="list")

WBvars <- doc[grep("variation/WBVar",doc)]
allelespan <- doc[grep("variation/WBVar",doc)+3] ## contains allele (or +2 if )
alls <- sapply(strsplit(allelespan,"[><]"), "[", 2)
toFix <- grep("/tr",alls)
fixalls <- sapply(strsplit(doc[grep("variation/WBVar",doc)+2],"[><]"), "[", 2)
alls[toFix] <- fixalls[toFix]


toSkip <- c("gk","ok","tm")


allstoSkip <- unique (grep(paste(toSkip,collapse="|"), alls, value=TRUE))

WBvars <- sapply(strsplit(WBvars,"[/\"]"), "[", 8)

WBvars <- WBvars[-match(allstoSkip,alls)]
# alls[-match(allstoSkip,alls)]


head <- "http://www.wormbase.org/rest/widget/variation/"

overview <- "/overview"
details <- "/molecular_details"
isolation <- "/isolation"

mutations <- data.frame(allele = 0, species = 0, type = 0, mutagen = 0, effect=0)

numvars <- length(WBvars)

for(i in 1:numvars){
    
  url <-paste0(head,WBvars[i],overview)
  e406 = scan(url,what="list")
  alleleblock <- e406[grep("locus",e406)[1]]
  allele <- strsplit(alleleblock, "[><]")[[1]][2]
  species <- sapply(strsplit(e406[grep("Species:",e406)+7],"<"), "[", 1)
  
  url <-paste0(head,WBvars[i],isolation)
  e406 <- scan(url,what="list")
  if(length(grep("evidence",e406))==0){ mutagen <- e406[grep("Mutagen:",e406)+5]}
  if(length(grep("evidence",e406))>0){  mutagen <- e406[grep("Mutagen:",e406)+8]}
  
  url <-paste0(head,WBvars[i],details)
  e406 <- scan(url,what="list")
  type <- e406[grep("Allele:",e406)+1]
  if (length(type)==0){type <- "-"}
  effect <- e406[grep("protein:",e406)+1]
  if (length(effect)>0) {
    effect <- sapply(strsplit(effect,"<"), "[", 1)
    effect <- paste0(effect, collapse=",")
  }
  if (length(effect)==0){effect <- "-"}
  
  scrape <- data.frame(allele = allele, species = species, type = type, mutagen = mutagen, effect=effect)
  print(c(i,allele,species,type,mutagen,effect))
  mutations<-rbind(mutations,scrape)
}


substitutions <- mutations[which(mutations$species=="elegans" & mutations$type == "substitution"),]
subswEffect <- substitutions[which(substitutions$effect != "-"),]

N_missense <- length(grep("Missense",subswEffect$effect))
N_nonsense <- length(grep("Nonsense",subswEffect$effect))
N_splice <- length(grep("Splice",subswEffect$effect))






