# combine the data raw RNAseq counts
logFile<-"/sc/orga/projects/changr04b/F980/marcRawCounts/GENESIPS_IRIS_validation_rawCounts_R.log"

cat(file=logFile,append=F,"This is GENESIPS_IRIS_RNAseq_Validation_combineData.R.\n\nCombining the feature count files below into the following raw count file: /sc/orga/projects/changr04b/F980/marcRawCounts/GENESIPS_IRIS_validation_rawCounts.tab.\n\t")

for(i in 1:24){
      cat(file=logFile,append=T,paste(sep="",paste(sep="","/sc/orga/projects/changr04b/F980/QC_B243.F980_Genesips_24_KDA_validation.SE.RNASeqPolyA.RAPiD.Human/Raw/RNA.IlluminaHiSeq2500.RNASeqPolyA/Sample_Genesips_",i,"/processed/featurecount/Genesips_",i,".primary.txt"),"\n\t"))
      featCountFile<-paste(sep="","/sc/orga/projects/changr04b/F980/QC_B243.F980_Genesips_24_KDA_validation.SE.RNASeqPolyA.RAPiD.Human/Raw/RNA.IlluminaHiSeq2500.RNASeqPolyA/Sample_Genesips_",i,"/processed/featurecount/Genesips_",i,".primary.txt")
      featCountTmp<-read.table(featCountFile,sep="\t",header=T,stringsAsFactors=F,comment.char="#")
      if(i==1){
	rawCounts<-featCountTmp
	colnames(rawCounts)[ncol(rawCounts)]<-"Genesips_1"
      }else{
	rawCounts<-cbind(rawCounts,featCountTmp[match(rawCounts[,1],featCountTmp[,1]),ncol(featCountTmp)])
	colnames(rawCounts)[ncol(rawCounts)]<-paste(sep="","Genesips_",i)
	if(sum(!is.element(el=featCountTmp[,1],set=rawCounts[,1]))>0){
		print("NEED TO WRITE EXTRA CODE - SORRY...")
	}
      }
}

write.table(rawCounts,sep="\t",file="/sc/orga/projects/changr04b/F980/marcRawCounts/GENESIPS_IRIS_validation_rawCounts.tab",row.names=F,col.names=T,quote=F)
cat(file=logFile,append=T,"\nThis is the end.\n\n")
