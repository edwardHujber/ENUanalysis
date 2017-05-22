library("readr")
library("seqinr")


## read in data

mmp <- read_delim(paste0(workPath,"mmpfixedsnpEff.snpeff"), delim="\t")


# ## discard nonunique variants. I dont think I want to do this since there are dup rows from snpeff. 
# mmp$CPRC <- apply(mmp,1,function(x) paste(x['CHROM'],x['POS'],x['REF'],x['ALT']))
# CPRC <- data.frame("CPRC" = mmp$CPRC )
# CPRCcount <- data.frame(table(CPRC))
# parental <- CPRCcount[which(CPRCcount$Freq>=cutoff),]      ## "parental" variants must show up >= cutoff number of times
# uniqueVars <- mmp[which(is.na(match(mmp$CPRC,parental$CPRC))),]   ## remove variants that match a parental


## split by VC number
# UV_TMP <- uniqueVars[grep("VC10", uniqueVars$strain),]
EMS <- mmp[grep("VC20", mmp$ID),]
ENU <- mmp[grep("VC30", mmp$ID),]

# EMS_ENU <- uniqueVars[grep("VC40", uniqueVars$strain),]

if(getDat=="ENU") {alldat <- ENU}
if(getDat=="EMS") {alldat <- EMS}

## fix the old_AA.new_AA column to jive better with Fig3
print("Fixing columns...")
AAchange <-apply(alldat,1,function(s){
   if(!grepl("/",s['EFF_Codon_Change'])| (nchar(s['EFF_Codon_Change'])!=7)){NA}  ## if EFF_Codon_Change isnt of the format '***/***', put NA in the AA change col
   else{
      codonChange <-strsplit(s['EFF_Codon_Change'],"/")[[1]]
      paste(sapply(codonChange,function(x){translate(s2c(x))}),collapse="/")     ## Otherwise, translate the codons to AA, return the */* format.
      }
})

## get things into uniform structure, with only the necessary columns
dat<-data.frame("ID"=alldat$ID,"chr"=alldat$CHROM,"pos"=alldat$POS,"wt_dna"=alldat$REF,"mut_dna"=alldat$ALT, "Old_codon.New_codon"=alldat$EFF_Codon_Change, "old_AA.new_AA"=AAchange)


rm(list=setdiff(ls(), c("dat",VARS2KEEP))) 
gc()


# #### THis is all done in Fig3 now
# effectClasses <- data.frame(table(dat$effect))
# effectClasses <- effectClasses[which(effectClasses$Freq >0 ),]
# 
# missense <- dat[which(dat$effect=="missense"),]
# 
# 
# ## which missense mutations fall into the same AA class?
# conserved <- missense[which(allAA$class[match(missense$wt_prot,allAA$AA)] == allAA$class[match(missense$mut_prot,allAA$AA)]),]
# unconserved <- missense[which(allAA$class[match(missense$wt_prot,allAA$AA)] != allAA$class[match(missense$mut_prot,allAA$AA)]),]
# 
# unconChangeClasses <- data.frame(allAA$class[match(missense$wt_prot,allAA$AA)],allAA$class[match(missense$mut_prot,allAA$AA)])
# colnames(unconChangeClasses)[1:2] <- c("wtClass","mutClass") 
# unconChangeClasses <- as.data.frame(table(unconChangeClasses))
# # unconChangeClasses <- unconChangeClasses[which(unconChangeClasses$wtClass != unconChangeClasses$mutClass),]
# unconChangeClasses <- unconChangeClasses[order(unconChangeClasses$wtClass, decreasing=FALSE),]


