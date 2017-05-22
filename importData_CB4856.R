## Thompson OA, Snoek LB, Nijveen H, et al. Remarkably Divergent Regions Punctuate the Genome Assembly of the Caenorhabditis elegans Hawaiian Strain CB4856. Genetics. 2015;200(3):975-989. doi:10.1534/genetics.115.175950.

CB4856 <- read.table(paste0(workPath,"File_S2_Celegans_N2_CB4856_SNVs_features.txt"), quote="\"")

fixedCodons<-apply(cbind(as.character(CB4856[,11]),CB4856[,12]),1,function(x){
   if(is.na(x[1])){return(x[1])}
   else{
      varPos<-as.numeric(x[2])+1
      wtCodon <- tolower(substr(x[1],1,3))
      mutCodon<- tolower(substr(x[1],6,8))
      substr(wtCodon,varPos,varPos)<-toupper(substr(wtCodon,varPos,varPos))
      substr(mutCodon,varPos,varPos)<-toupper(substr(mutCodon,varPos,varPos))
      return(paste0(wtCodon,"/",mutCodon))
   }
})

fixedAA<-gsub("->","/",CB4856[,13])

dat <- data.frame("ID"="CB4856","chr"=CB4856[,1],"pos"=CB4856[,2],"wt_dna"=CB4856[,3],"mut_dna"=CB4856[,5],"Old_codon.New_codon"=fixedCodons,"old_AA.new_AA"=fixedAA)

rm(list=setdiff(ls(), c("dat",VARS2KEEP))) 
gc()