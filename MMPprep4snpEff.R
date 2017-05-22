library("readr")
library("data.table")

mmp <- read.delim(paste0(workPath,"mmp_mut_strains_data_Mar14.txt"),  stringsAsFactors=FALSE)

mmp <- mmp[which((nchar(mmp$wt_dna)==1) & nchar(mmp$mut_dna)==1),]  ## drop indels

vcf <- data.frame("CHROM"=mmp$chr,"POS"=mmp$pos,"ID"=mmp$strain,"REF"=mmp$wt_dna,"ALT"=mmp$mut_dna,"QUAL"="30","FILTER"="pass","INFO"=".")
names(vcf)[1]<-paste0("#",names(vcf)[1])

write.table(vcf,paste0(workPath,"mmpVCF.vcf"),sep="\t",quote=FALSE,row.names=FALSE)


# javaCall <- paste0("java -Xms1024M -jar C:\\snpEff\\snpEff.jar eff -v -formatEff -no-downstream -no-intergenic -no-intron -no-upstream  WBcel235.82 \"Z:\\R\\ENU analysis\\mmpVCF.vcf\" > \"Z:\\R\\ENU analysis\\mmp.snpeff\"")
# system(javaCall)

Ncomments <- sum(substr(readLines(paste0(workPath,"mmp.snpeff"),n=100),1,1)=="#")
comment <- readLines(paste0(workPath,"mmp.snpeff"), n=Ncomments)


eff <- read_delim(paste0(workPath,"mmp.snpeff"), delim="\t", comment="", skip = Ncomments-1)
colnames(eff)[1]<-substr(colnames(eff)[1],2,nchar(colnames(eff)[1])) ## gets rid of the # at begin of the first col name

infoLine <-grep("##INFO",comment) ## index of lines with annotation info

## makes a list of ann names for each anngroup
annNames <- list()
j<-0
for (i in infoLine){
j<-j+1
annID <- strsplit(comment[i], "[=,]")[[1]][3] ## = 'EFF' or whatever
annFormatpre <- strsplit(comment[i], "Format: '")[[1]][2]
annFormatpro<- strsplit(annFormatpre, "[()]")[[1]][which(nchar(strsplit(annFormatpre, "[()]")[[1]])==max(nchar(strsplit(annFormatpre, "[()]")[[1]])))] ## pulls out the stuff between parenthesis (if theyre there) by assuming its thhe longest thing
annFormat<-gsub(" ","",annFormatpro)
annNames[[j]] <- c(annID,paste0(annID,"_",unlist(strsplit(annFormat,"\\|")))) ## chr array of effects with prefix
}




##################
# effsmall <- eff[1:1000,]
eff1<-eff[2,]
eff<-eff1
gc()
pb<-txtProgressBar(max=nrow(eff),style=3,width=100)
effectedList<-apply(eff,1,function(x){
   i<-which(eff$ID==x['ID'])[1]
   setTxtProgressBar(pb,i)
   splitInfo <-lapply(strsplit(x['INFO'],";")[[1]],function(y){y})
   print(splitInfo)
   
   splitEff <-strsplit(splitInfo[[1]],",")
   lineList <-lapply(X=splitEff[[1]],FUN=function(X,y,z){
      effLine <- X
      if(as.logical(length(grep("=",effLine)))){ ## true if = is there
         effLine <- substring(effLine,5)
      }
      effects <-strsplit(effLine,"[(\\|)]")
      if(class(effects)=="list"){effects<-unlist(effects)}
      while(length(effects)!=length(annNames[[1]])){effects<-c(effects,"")} 
      names(effects)<-annNames[[1]]
      y <- c(z[-which(names(z)=="INFO")],effects)
   },y=effout,z=x)
   lineDF<-data.frame(matrix(unlist(lineList), nrow=length(lineList), byrow=T),stringsAsFactors=FALSE)
})
close(pb)
## Took 2.224077 hours

effectedDF<- rbindlist(effectedList)
## Took 8.920444 secs

newNames <- c(names(eff)[-which(names(eff)=="INFO")],annNames[[1]])
names(effectedDF)<-newNames

write.table(effectedDF,paste0("C:/Users/image_analysis/Desktop/","mmpfixedsnpEff.snpeff"),sep="\t",quote=FALSE,row.names=FALSE)
## Took 18.51829 secs


