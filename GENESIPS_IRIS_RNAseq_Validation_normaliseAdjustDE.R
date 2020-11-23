##########
# set-up #
##########

rm(list = ls())
library(limma)
library(edgeR)
library(ggplot2)
library(variancePartition)
library(RColorBrewer)
library(doParallel)

source("/hpc/users/xhenrim01/VRC/scripts/R_scripts/usefulFunctionsForRNAseqAnalysis.R")
source("/hpc/users/xhenrim01/VRC/scripts/R_scripts/rnaSeqLimmaPlotPCAFunction_general.R")
source("/hpc/users/xhenrim01/VRC/scripts/R_scripts/msigDB_enrichment_modified_subroutineOnly.R")

registerDoParallel(10)

args<-commandArgs(TRUE)
date2use<-args[1]
doVarPartPlots<-as.logical(args[2])
skip2GSEA<-as.logical(args[3])

if(is.na(date2use)){date2use<-Sys.Date()}
if(is.na(doVarPartPlots)){doVarPartPlots<-FALSE}
if(is.na(skip2GSEA)){skip2GSEA<-FALSE}

# filtering thresholds set below in the definition of function getGeneFilteredGeneExprMatrix
minGeneCpm=1
minSamplePercentWithMinGeneCpm=0.3
# same parameters as what was used for the main / discovery analysis

workDir<-"/sc/orga/projects/changr04b/F980/"
outDir<-paste(sep="/",workDir,"marcNormalisedAdjustedDE",date2use)
outPrefix<-paste(sep="/",outDir,"genesipsIrisValidationRnaseqAtorvastatin")
#logFile<-paste(sep="",outPrefix,"_R.log")
logFile<-""

setwd(workDir)
system2(command="mkdir",args=c("-p",outDir))

cat(file=logFile,append=F,"This is GENESIPS_IRIS_RNAseq_Validation_normaliseAdjustDE.R\n\nInput parameters:\n\tdate2use = < ",date2use," >\n\tdoVarPartPlots = < ",doVarPartPlots," >\n\tminGeneCpm = < ",minGeneCpm," >\n\tminSamplePercentWithMinGeneCpm = < ",minSamplePercentWithMinGeneCpm," >\n\toutPrefix = < ",outPrefix," >\n\tskip2GSEA = < ",skip2GSEA," >.\n\n")

#hgncTable<-read.table("/sc/orga/projects/changr04a/marc_VRC/hgnc_complete_set.txt",header=T,sep="\t",stringsAsFactors=F,comment.char="",quote="")
#ensembl70<-read.table(sep="\t",header=T,stringsAsFactors=F,"/sc/orga/projects/changr04a/marc_VRC/Ensembl/Homo_sapiens.GRCh37.70.gtf.flat.withBiotype.tab")

countFile<-"marcRawCounts/GENESIPS_IRIS_validation_rawCounts.tab"
covFile<-"ivanPaigeCovariates/covariateTableIpscIrisAtorvastatinFinalWithRin.csv"
msigdb_C5BP_File<-"/sc/orga/projects/changr04b/F980/mSigDB_gmts/c5.bp.v6.1.symbols.gmt"
msigdb_C5ALL_File<-"/sc/orga/projects/changr04b/F980/mSigDB_gmts/c5.all.v6.1.symbols.gmt"
outliers<-character(0)

## NOTE: script from the discovery data analysis
## /hpc/users/xhenrim01/VRC/scripts/RuiCode/R/RNAseqDataExplore.GENESIPS.IRvsIS.v6.317.88.short.R


###############################
# START PROCESSING & ANALYSIS #
###############################

if(!skip2GSEA){
# read the count data
cat(file=logFile,append=T,"Reading the count data and removing lowly expressed genes.\n")
rawcounts<-read.table(sep="\t",stringsAsFactor=F,header=T,countFile)
genenames<-rawcounts[,1]
geneinfo<-rawcounts[,1:6]
sample<-rawcounts[,-c(1:6)]
rownames(sample)<-genenames
idx<-match(as.character(outliers),colnames(sample))
if(length(idx)>0){sample<-sample[,-idx]}
filteredSample<-getGeneFilteredGeneExprMatrix(sample)
save(filteredSample,file=paste(sep="",outPrefix,"_filteredRNAseqCounts.RData"))
cat(file=logFile,append=T,paste(sep="","Done - read ",ncol(sample)," samples and ",nrow(sample)," genes.\n\n"))

# read the covariate data
cat(file=logFile,append=T,"Reading the covariate data.\n")
info<-read.csv(covFile,header=T,stringsAsFactor=F)
newinfo=info[match(as.character(colnames(filteredSample)),info$SampleName),]
rownames(newinfo)<-newinfo$SampleName
newinfo$State<-ifelse(newinfo$State=="Insulin resistant","IR","IS")
newinfo$State<-factor(newinfo$State)
newinfo$Indiv<-paste(sep="","indiv",newinfo$Indiv)
newinfo$Indiv<-factor(newinfo$Indiv)
newinfo$Atorvastatin<-factor(gsub(pattern="[0-9]+-[0-9]-",replacement="",newinfo$NoNameCovariate))
save(newinfo,file=paste(sep="",outPrefix,"_all_covariates_newinfo.RData"))
cat(file=logFile,append=T,"Done - read ",nrow(newinfo)," samples and ",ncol(newinfo)," covariates.\n\n")
	
# normalise data - no covariates except State specified in design
cat(file=logFile,append=T,"Normalising data without taking any covariates into account for now (but with blocking on Indiv).\n")
vobj<-doNorma(expDat=filteredSample,saveNorma=TRUE,outPrefix=outPrefix,blockVar=newinfo$Indiv)
cat(file=logFile,append=T,"Done.\n\n")

# plot of mislabelled samples
cat(file=logFile,append=T,"Checking recorded and measured sex.\n")
#UTY="ENSG00000183878"
#XIST="ENSG00000229807"
UTY="UTY"
XIST="XIST"
pdf(paste(sep="",outPrefix,"_normalisedNotAdjusted_SexCheck.pdf"))
plot(vobj$E[UTY,],vobj$E[XIST,])
points(vobj$E[UTY,newinfo$SEX=="Male"],vobj$E[XIST,newinfo$SEX=="Male"],col="red")
points(vobj$E[UTY,newinfo$SEX=="Female"],vobj$E[XIST,newinfo$SEX=="Female"],col="blue") # the outlying samples at the bottom left quadrant are from clones iE231_3 and iE255_1; manually checked the samples; they are NOT failed samples but the low sex genes expression is probably due to the process of making them pluripotent; they will be retained
dev.off()
	
idx.ambiguous<-which((vobj$E[UTY,]<0 & vobj$E[XIST,]<0) | (vobj$E[UTY,]>0 & vobj$E[XIST,]>0))
print("The following samples have ambiguous sex gene expression levels:")
print(idx.ambiguous)
print(newinfo[idx.ambiguous,])
cat(file=logFile,append=T,"Done.\n\n")

# hiearchical clustering
cat(file=logFile,append=T,"Clustering the samples to check for outliers.\n")
doHierarchClust(expDat=vobj$E,colLabels=T,colInfo=newinfo$State,outFile=paste(sep="",outPrefix,"_normalisedNotAdjusted_hierarchClust.pdf"))
cat(file=logFile,append=T,"Done.\n\n")

# do PCA plots
cat(file=logFile,append=T,"PCA plots to decided what covariates to adjust for.\n")
doPcaPlots(expDat=vobj$E,covDat=newinfo,outFile=paste(sep="",outPrefix,"_normalisedNotAdjusted_PCA.pdf"),covarList=c("State","SEX","Age","Atorvastatin","BMI","Indiv","Sequencing.Batch","Reprogramming.Batch","RNA.Harvest.Date"),covarTypeList=c("factor","factor","numeric","factor","numeric","factor","factor","factor","factor"))
# Note: RACE, Total.Reads, Reprogramming.Source.Cell, Sendai.Virus.Lot, RIN, RNA.Preparation.Date either all had the same value, or so few values that these would have been captured by Indiv or were completely missing
# Note (after PCA plotting): NoNameCovariate has a different value for each sample (hence dropped), Reprogramming.Batch has 1 distinct value for every individuals (hence dropped because captured by Indiv) 
cat(file=logFile,append=T,"Done.\n\n")

# variance partition on normalised, non-adjusted data
cat(file=logFile,append=T,"Variance partition analysis to decide covariates to adjust for.\n")
form.varPart <- ~ (1|State) + (1|Indiv) + (1|Sequencing.Batch) + (1|RNA.Harvest.Date) + (1|SEX) + (1|Atorvastatin) + Age + BMI
if(doVarPartPlots){
	doVarPart(expDat=vobj,covDat=newinfo,form=form.varPart,outPrefix=paste(sep="",outPrefix,"_normalisedNotAdjusted_varPart"))
}
cat(file=logFile,append=T,"Done.\n\n")

# adjust for covariates / compute expression residuals -- variancePartition for Sequencing.Batch and RNA.Harvest.Date first, then lmFit for SEX, Age, BMI [same adjustment procedure as for the discovery data]

## A: using variancePartition for Sequencing.Batch and RNA.Harvest.Date, lmFit for SEX, Age, CMI and accounting for repeated obs. per patient as random effect for Indiv (blocking)
cat(file=logFile,append=T,"Adjusting for covariates and computing residuals -- method A.\n")
cat(file=logFile,append=T,".. step1: variancePartition for Sequencing.Batch and RNA.Harvest.Date\n")
resid.vp<-getResiduals(expDat=vobj,covDat=newinfo,form=~(1|Sequencing.Batch)+(1|RNA.Harvest.Date),outPrefix=outPrefix,adjCovs="Sequencing.Batch-RNA.Harvest.Date",voomFirst=FALSE,useVarPar=TRUE,blockVar=NULL,returnFit=TRUE,saveNorma=F,dupCorDoneAlready=TRUE)
cat(file=logFile,append=T,".. step2: lmFit for SEX, Age, BMI - blocking on Indiv\n")
resid.vplm<-getResiduals(expDat=resid.vp$resid,covDat=newinfo,form=~SEX+Age+BMI,outPrefix=paste(sep="",outPrefix,"_residuals-VP_Sequencing.Batch-RNA.Harvest.Date"),adjCovs="SEX-Age-BMI",voomFirst=FALSE,useVarPar=FALSE,blockVar=newinfo$Indiv,returnFit=TRUE,saveNorma=F,dupCorDoneAlready=TRUE)
cat(file=logFile,append=T,"Done.\n\n")

## B: using lmFit for Sequencing.Batch, RNA.Harvest.Date, SEX, Age, CMI and accounting for repeated obs. per patient as random effect for Indiv (blocking)
cat(file=logFile,append=T,"Adjusting for covariates and computing residuals -- method B.\n")
cat(file=logFile,append=T,".. lmFit for Sequencing.Batch, RNA.Harvest.Date, SEX, Age, BMI - blocking on Indiv\n")
resid.lm<-getResiduals(expDat=vobj,covDat=newinfo,form=~Sequencing.Batch+SEX+Age+BMI,outPrefix=outPrefix,adjCovs="Sequencing.Batch-RNA.Harvest.Date-SEX-Age-BMI_blockIndiv",voomFirst=FALSE,useVarPar=FALSE,blockVar=newinfo$Indiv,returnFit=TRUE,saveNorma=F,dupCorDoneAlready=TRUE)
cat(file=logFile,append=T,"Done.\n\n")

## C: using lmFit for Sequencing.Batch, RNA.Harvest.Date, SEX, Age, CMI WITHOUT accounting for repeated obs. per patient
cat(file=logFile,append=T,"Adjusting for covariates and computing residuals -- method C.\n")
cat(file=logFile,append=T,".. lmFit for Sequencing.Batch, RNA.Harvest.Date, SEX, Age, BMI - NO blocking on Indiv\n")
resid.lm.noblock<-getResiduals(expDat=vobj,covDat=newinfo,form=~Sequencing.Batch+SEX+Age+BMI,outPrefix=outPrefix,adjCovs="Sequencing.Batch-RNA.Harvest.Date-SEX-Age-BMI",voomFirst=FALSE,useVarPar=FALSE,blockVar=NULL,returnFit=TRUE,saveNorma=F,dupCorDoneAlready=TRUE)
cat(file=logFile,append=T,"Done.\n\n")

# redo PCA, variancePartition, clustering
# do PCA plots
cat(file=logFile,append=T,"PCA plots to check residuals\n")
doPcaPlots(expDat=resid.vplm$resid,covDat=newinfo,outFile=paste(sep="",outPrefix,"_residuals-VP_Sequencing.Batch-RNA.Harvest.Date_residuals-lmFit_SEX-Age-BMI_PCA.pdf"),covarList=c("State","SEX","Age","Atorvastatin","BMI","Indiv","Sequencing.Batch","Reprogramming.Batch","RNA.Harvest.Date"),covarTypeList=c("factor","factor","numeric","factor","numeric","factor","factor","factor","factor"))

doPcaPlots(expDat=resid.lm$resid,covDat=newinfo,outFile=paste(sep="",outPrefix,"_residuals-lmFit_Sequencing.Batch-RNA.Harvest.Date-SEX-Age-BMI_blockIndiv_PCA.pdf"),covarList=c("State","SEX","Age","Atorvastatin","BMI","Indiv","Sequencing.Batch","Reprogramming.Batch","RNA.Harvest.Date"),covarTypeList=c("factor","factor","numeric","factor","numeric","factor","factor","factor","factor"))

doPcaPlots(expDat=resid.lm.noblock$resid,covDat=newinfo,outFile=paste(sep="",outPrefix,"_residuals-lmFit_Sequencing.Batch-RNA.Harvest.Date-SEX-Age-BMI_PCA.pdf"),covarList=c("State","SEX","Age","Atorvastatin","BMI","Indiv","Sequencing.Batch","Reprogramming.Batch","RNA.Harvest.Date"),covarTypeList=c("factor","factor","numeric","factor","numeric","factor","factor","factor","factor"))
cat(file=logFile,append=T,"Done.\n\n")

# variance partition on normalised, non-adjusted data
cat(file=logFile,append=T,"Variance partition analysis to check residuals.\n")
form.varPart <- ~ (1|State) + (1|Indiv) + (1|Sequencing.Batch) + (1|RNA.Harvest.Date) + (1|SEX) + (1|Atorvastatin) + Age + BMI
if(doVarPartPlots){
	doVarPart(expDat=resid.vplm$resid,covDat=newinfo,form=form.varPart,outPrefix=paste(sep="",outPrefix,"_residuals-VP_Sequencing.Batch-RNA.Harvest.Date_residuals-lmFit_SEX-Age-BMI_varPart"))
        doVarPart(expDat=resid.lm$resid,covDat=newinfo,form=form.varPart,outPrefix=paste(sep="",outPrefix,"_residuals-lmFit_Sequencing.Batch-RNA.Harvest.Date-SEX-Age-BMI_blockIndiv_varPart"))
        doVarPart(expDat=resid.lm.noblock$resid,covDat=newinfo,form=form.varPart,outPrefix=paste(sep="",outPrefix,"_residuals-lmFit_Sequencing.Batch-RNA.Harvest.Date-SEX-Age-BMI_varPart"))
}
cat(file=logFile,append=T,"Done.\n\n")

# hiearchical clustering
cat(file=logFile,append=T,"Clustering the samples to check for outliers.\n")
doHierarchClust(expDat=resid.vplm$resid,colLabels=T,colInfo=newinfo$State,outFile=paste(sep="",outPrefix,"_residuals-VP_Sequencing.Batch-RNA.Harvest.Date_residuals-lmFit_SEX-Age-BMI_hierarchClust.pdf"))
doHierarchClust(expDat=resid.lm$resid,colLabels=T,colInfo=newinfo$State,outFile=paste(sep="",outPrefix,"_residuals-lmFit_Sequencing.Batch-RNA.Harvest.Date-SEX-Age-BMI_blockIndiv_hierarchClust.pdf"))
doHierarchClust(expDat=resid.lm.noblock$resid,colLabels=T,colInfo=newinfo$State,outFile=paste(sep="",outPrefix,"_residuals-lmFit_Sequencing.Batch-RNA.Harvest.Date-SEX-Age-BMI_hierarchClust.pdf"))
cat(file=logFile,append=T,"Done.\n\n")
    
# DE analysis
deAtorPairedFUN<-function(resid,info,outPrefix){
    resid_D<-resid[,match(info$SampleName[info$Atorvastatin=="D"],colnames(resid))]
    resid_A<-resid[,match(info$SampleName[info$Atorvastatin=="A"],colnames(resid))]

    resid_D_IR<-resid[,match(info$SampleName[info$Atorvastatin=="D" & info$State=="IR"],colnames(resid))]
    resid_A_IR<-resid[,match(info$SampleName[info$Atorvastatin=="A" & info$State=="IR"],colnames(resid))]

    resid_D_IS<-resid[,match(info$SampleName[info$Atorvastatin=="D" & info$State=="IS"],colnames(resid))]
    resid_A_IS<-resid[,match(info$SampleName[info$Atorvastatin=="A" & info$State=="IS"],colnames(resid))]

    de_pairwise<-data.frame(gene=rownames(resid),n_pairs=NA,mean_A=NA,mean_D=NA,meanPairedLog2FoldChange=NA,t_stat=NA,t_pval=NA,t_pval_adj=NA,median_A=NA,median_D=NA,medianPairedLog2FoldChange=NA,wilcoxon_pval=NA,wilcoxon_pval_adj=NA)
    de_pairwise_IR<-data.frame(gene=rownames(resid),n_pairs=NA,mean_A=NA,mean_D=NA,meanPairedLog2FoldChange=NA,t_stat=NA,t_pval=NA,t_pval_adj=NA,median_A=NA,median_D=NA,medianPairedLog2FoldChange=NA,wilcoxon_pval=NA,wilcoxon_pval_adj=NA)
    de_pairwise_IS<-data.frame(gene=rownames(resid),n_pairs=NA,mean_A=NA,mean_D=NA,meanPairedLog2FoldChange=NA,t_stat=NA,t_pval=NA,t_pval_adj=NA,median_A=NA,median_D=NA,medianPairedLog2FoldChange=NA,wilcoxon_pval=NA,wilcoxon_pval_adj=NA)

    computeDEmat_onerow<-function(deMat,i,resid_condA,resid_condB){
        t.tmp<-t.test(resid_condA[i,],resid_condB[i,],paired=T)
        w.tmp<-wilcox.test(resid_condA[i,],resid_condB[i,],paired=T,exact=T)
        n_pair<-t.tmp$parameter+1 # df + 1
        diff_mean<-t.tmp$estimate
        t_stat<-t.tmp$statistic
        t_pval<-t.tmp$p.value
        wilcoxon_pval<-w.tmp$p.value
        return(c(n_pair,diff_mean,t_stat,t_pval,wilcoxon_pval))
    }

    computeDEmat_statsFullMat<-function(deMat,resid_condA,resid_condB){
        tmpMat<-resid_condA-resid_condB
        deMat$medianPairedLog2FoldChange<-apply(FUN=median,MARGIN=1,X=tmpMat,na.rm=T) # median of the difference, not difference of the medians!
        deMat$mean_A<-rowMeans(ifelse(!is.na(resid_condA) & !is.na(resid_condB),resid_condA,NA))
        deMat$mean_D<-rowMeans(ifelse(!is.na(resid_condA) & !is.na(resid_condB),resid_condB,NA))
        deMat$median_A<-apply(FUN=median,MARGIN=1,X=ifelse(!is.na(resid_condA) & !is.na(resid_condB),resid_condA,NA))
        deMat$median_D<-apply(FUN=median,MARGIN=1,X=ifelse(!is.na(resid_condA) & !is.na(resid_condB),resid_condB,NA))
        return(deMat)
    }

    de_pairwise<-computeDEmat_statsFullMat(deMat=de_pairwise,resid_condA=resid_A,resid_condB=resid_D)
    de_pairwise_IR<-computeDEmat_statsFullMat(deMat=de_pairwise_IR,resid_condA=resid_A_IR,resid_condB=resid_D_IR)
    de_pairwise_IS<-computeDEmat_statsFullMat(deMat=de_pairwise_IS,resid_condA=resid_A_IS,resid_condB=resid_D_IS)

    idx.Col<-match(c("n_pairs","meanPairedLog2FoldChange","t_stat","t_pval","wilcoxon_pval"),colnames(de_pairwise))
    for(i in 1:nrow(resid)){
        de_pairwise[i,idx.Col]<-computeDEmat_onerow(deMat=de_pairwise,i=i,resid_condA=resid_A,resid_condB=resid_D)
        de_pairwise_IR[i,idx.Col]<-computeDEmat_onerow(deMat=de_pairwise_IR,i=i,resid_condA=resid_A_IR,resid_condB=resid_D_IR)
        de_pairwise_IS[i,idx.Col]<-computeDEmat_onerow(deMat=de_pairwise_IS,i=i,resid_condA=resid_A_IS,resid_condB=resid_D_IS)        
    }
    
    de_pairwise$t_pval_adj<-p.adjust(de_pairwise$t_pval,method="BH")
    de_pairwise$wilcoxon_pval_adj<-p.adjust(de_pairwise$wilcoxon_pval,method="BH")

    de_pairwise_IR$t_pval_adj<-p.adjust(de_pairwise_IR$t_pval,method="BH")
    de_pairwise_IR$wilcoxon_pval_adj<-p.adjust(de_pairwise_IR$wilcoxon_pval,method="BH")

    de_pairwise_IS$t_pval_adj<-p.adjust(de_pairwise_IS$t_pval,method="BH")
    de_pairwise_IS$wilcoxon_pval_adj<-p.adjust(de_pairwise_IS$wilcoxon_pval,method="BH")

    de_pairwise<-de_pairwise[order(de_pairwise$t_pval,decreasing=F),]
    de_pairwise_IR<-de_pairwise_IR[order(de_pairwise_IR$t_pval,decreasing=F),]
    de_pairwise_IS<-de_pairwise_IS[order(de_pairwise_IS$t_pval,decreasing=F),]
    
    write.table(de_pairwise,sep="\t",row.names=F,col.names=T,quote=F,file=paste(sep="",outPrefix,"_all.tab"))
    write.table(de_pairwise_IR,sep="\t",row.names=F,col.names=T,quote=F,file=paste(sep="",outPrefix,"_IR.tab"))
    write.table(de_pairwise_IS,sep="\t",row.names=F,col.names=T,quote=F,file=paste(sep="",outPrefix,"_IS.tab"))    
}

cat(file=logFile,append=T,"Do the DE analysis - pairwise tests.\n")
deAtorPairedFUN(resid=resid.vplm$resid,info=newinfo,outPrefix=paste(sep="",outPrefix,"_DE_pairwise_residuals-VP_Sequencing.Batch-RNA.Harvest.Date_residuals-lmFit_SEX-Age-BMI"))
deAtorPairedFUN(resid=resid.lm$resid,info=newinfo,outPrefix=paste(sep="",outPrefix,"_DE_pairwise_residuals-lmFit_Sequencing.Batch-RNA.Harvest.Date-SEX-Age-BMI_blockIndiv"))
deAtorPairedFUN(resid=resid.lm.noblock$resid,info=newinfo,outPrefix=paste(sep="",outPrefix,"_DE_pairwise_residuals-lmFit_Sequencing.Batch-RNA.Harvest.Date-SEX-Age-BMI"))

#topSetDE_all<-doDE("AtorvastatinA-AtorvastatinD",expDat=resid.vplm$resid,covDat=newinfo,form=~0+Atorvastatin,outPrefix=paste(sep="",outPrefix,"_DE_AtorvastatinAvsD"),blockVar=newinfo$Indiv,returnTopSet=TRUE,ENSG=F)
#idxIR<-which(newinfo$State=="IR")
#topSetDE_IR<-doDE("AtorvastatinA-AtorvastatinD",expDat=resid.vplm$resid[idxIR],covDat=newinfo[idxIR,],form=~0+Atorvastatin,outPrefix=paste(sep="",outPrefix,"_DE_AtorvastatinAvsD_IRonly"),blockVar=newinfo$Indiv[idxIR],returnTopSet=TRUE,ENSG=F)
#idxIS<-which(newinfo$State=="IS")
#topSetDE_IS<-doDE("AtorvastatinA-AtorvastatinD",expDat=resid.vplm$resid[idxIS],covDat=newinfo[idxIS,],form=~0+Atorvastatin,outPrefix=paste(sep="",outPrefix,"_DE_AtorvastatinAvsD_ISonly"),blockVar=newinfo$Indiv[idxIS],returnTopSet=TRUE,ENSG=F)
cat(file=logFile,append=T,"Done.\n\n")
}

# gene set enrichment analysis (since the results are not veyr much different for the 3 different ways of adjusting for covariates, we only focus on the first set of residuals here)
cat(file=logFile,append=T,"Doing gene set enrichment analyses.\n")

allDE<-read.table(sep="\t",header=T,stringsAsFactors=F,paste(sep="",outPrefix,"_DE_pairwise_residuals-VP_Sequencing.Batch-RNA.Harvest.Date_residuals-lmFit_SEX-Age-BMI_all.tab"))
irDE<-read.table(sep="\t",header=T,stringsAsFactors=F,paste(sep="",outPrefix,"_DE_pairwise_residuals-VP_Sequencing.Batch-RNA.Harvest.Date_residuals-lmFit_SEX-Age-BMI_IR.tab"))
isDE<-read.table(sep="\t",header=T,stringsAsFactors=F,paste(sep="",outPrefix,"_DE_pairwise_residuals-VP_Sequencing.Batch-RNA.Harvest.Date_residuals-lmFit_SEX-Age-BMI_IS.tab"))

sigAll<-allDE$gene[allDE$t_pval_adj<0.05]
sigIR<-allDE$gene[irDE$t_pval_adj<0.05]
sigIS<-allDE$gene[isDE$t_pval_adj<0.05]

diffIrIs<-setdiff(sigIR,allDE$gene[isDE$t_pval_adj<0.1])
diffIsIr<-setdiff(sigIS,allDE$gene[irDE$t_pval_adj<0.1])

write.table(sigAll,sep="\t",row.names=F,col.names=T,quote=F,file=paste(sep="",outPrefix,"_DE_pairwise_residuals-VP_Sequencing.Batch-RNA.Harvest.Date_residuals-lmFit_SEX-Age-BMI_all_Padj0.05.tab"))
write.table(sigIR,sep="\t",row.names=F,col.names=T,quote=F,file=paste(sep="",outPrefix,"_DE_pairwise_residuals-VP_Sequencing.Batch-RNA.Harvest.Date_residuals-lmFit_SEX-Age-BMI_IR_Padj0.05.tab"))
write.table(sigIS,sep="\t",row.names=F,col.names=T,quote=F,file=paste(sep="",outPrefix,"_DE_pairwise_residuals-VP_Sequencing.Batch-RNA.Harvest.Date_residuals-lmFit_SEX-Age-BMI_IS_Padj0.05.tab"))
write.table(diffIrIs,sep="\t",row.names=F,col.names=T,quote=F,file=paste(sep="",outPrefix,"_DE_pairwise_residuals-VP_Sequencing.Batch-RNA.Harvest.Date_residuals-lmFit_SEX-Age-BMI_Diff_IRPadj0.05_ISPadj0.1.tab"))
write.table(diffIsIr,sep="\t",row.names=F,col.names=T,quote=F,file=paste(sep="",outPrefix,"_DE_pairwise_residuals-VP_Sequencing.Batch-RNA.Harvest.Date_residuals-lmFit_SEX-Age-BMI_Diff_ISPadj0.05_IRPadj0.1.tab"))    
    
bgList<-allDE$gene # same for all three analyses

if(length(sigAll)>1){
    msigdbAll<-gsea_enrichment(geneList=sigAll,backgroundList=bgList,list.gs=msigdb_C5BP_File)
    write.table(sep="\t",row.names=F,col.names=T,quote=F,msigdbAll,file=paste(sep="",outPrefix,"_DE_pairwise_residuals-VP_Sequencing.Batch-RNA.Harvest.Date_residuals-lmFit_SEX-Age-BMI_all_Padj0.05_msigDB_C5BP.tab"))
    write.table(sep="\t",row.names=F,col.names=T,quote=F,msigdbAll[as.numeric(as.character(msigdbAll$Adjusted.Pvalue))<0.05,],file=paste(sep="",outPrefix,"_DE_pairwise_residuals-VP_Sequencing.Batch-RNA.Harvest.Date_residuals-lmFit_SEX-Age-BMI_all_Padj0.05_msigDB_C5BP_Padj0.05.tab"))    
    msigdbAll<-gsea_enrichment(geneList=sigAll,backgroundList=bgList,list.gs=msigdb_C5ALL_File)
    write.table(sep="\t",row.names=F,col.names=T,quote=F,msigdbAll,file=paste(sep="",outPrefix,"_DE_pairwise_residuals-VP_Sequencing.Batch-RNA.Harvest.Date_residuals-lmFit_SEX-Age-BMI_all_Padj0.05_msigDB_C5ALL.tab"))
    write.table(sep="\t",row.names=F,col.names=T,quote=F,msigdbAll[as.numeric(as.character(msigdbAll$Adjusted.Pvalue))<0.05,],file=paste(sep="",outPrefix,"_DE_pairwise_residuals-VP_Sequencing.Batch-RNA.Harvest.Date_residuals-lmFit_SEX-Age-BMI_all_Padj0.05_msigDB_C5ALL_Padj0.05.tab"))
}
if(length(sigIR)>1){
    msigdbIR<-gsea_enrichment(geneList=sigIR,backgroundList=bgList,list.gs=msigdb_C5BP_File)
    write.table(sep="\t",row.names=F,col.names=T,quote=F,msigdbIR,file=paste(sep="",outPrefix,"_DE_pairwise_residuals-VP_Sequencing.Batch-RNA.Harvest.Date_residuals-lmFit_SEX-Age-BMI_IR_Padj0.05_msigDB_C5BP.tab"))
    write.table(sep="\t",row.names=F,col.names=T,quote=F,msigdbIR[as.numeric(as.character(msigdbIR$Adjusted.Pvalue))<0.05,],file=paste(sep="",outPrefix,"_DE_pairwise_residuals-VP_Sequencing.Batch-RNA.Harvest.Date_residuals-lmFit_SEX-Age-BMI_IR_Padj0.05_msigDB_C5BP_Padj0.05.tab"))
    msigdbIR<-gsea_enrichment(geneList=sigIR,backgroundList=bgList,list.gs=msigdb_C5ALL_File)
    write.table(sep="\t",row.names=F,col.names=T,quote=F,msigdbIR,file=paste(sep="",outPrefix,"_DE_pairwise_residuals-VP_Sequencing.Batch-RNA.Harvest.Date_residuals-lmFit_SEX-Age-BMI_IR_Padj0.05_msigDB_C5ALL.tab"))
    write.table(sep="\t",row.names=F,col.names=T,quote=F,msigdbIR[as.numeric(as.character(msigdbIR$Adjusted.Pvalue))<0.05,],file=paste(sep="",outPrefix,"_DE_pairwise_residuals-VP_Sequencing.Batch-RNA.Harvest.Date_residuals-lmFit_SEX-Age-BMI_IR_Padj0.05_msigDB_C5ALL_Padj0.05.tab"))
}
if(length(sigIS)>1){
    msigdbIS<-gsea_enrichment(geneList=sigIS,backgroundList=bgList,list.gs=msigdb_C5BP_File)
    write.table(sep="\t",row.names=F,col.names=T,quote=F,msigdbIS,file=paste(sep="",outPrefix,"_DE_pairwise_residuals-VP_Sequencing.Batch-RNA.Harvest.Date_residuals-lmFit_SEX-Age-BMI_IS_Padj0.05_msigDB_C5BP.tab"))
    write.table(sep="\t",row.names=F,col.names=T,quote=F,msigdbIS[as.numeric(as.character(msigdbIS$Adjusted.Pvalue))<0.05,],file=paste(sep="",outPrefix,"_DE_pairwise_residuals-VP_Sequencing.Batch-RNA.Harvest.Date_residuals-lmFit_SEX-Age-BMI_IS_Padj0.05_msigDB_C5BP_Padj0.05.tab"))
    msigdbIS<-gsea_enrichment(geneList=sigIS,backgroundList=bgList,list.gs=msigdb_C5ALL_File)
    write.table(sep="\t",row.names=F,col.names=T,quote=F,msigdbIS,file=paste(sep="",outPrefix,"_DE_pairwise_residuals-VP_Sequencing.Batch-RNA.Harvest.Date_residuals-lmFit_SEX-Age-BMI_IS_Padj0.05_msigDB_C5ALL.tab"))
    write.table(sep="\t",row.names=F,col.names=T,quote=F,msigdbIS[as.numeric(as.character(msigdbIS$Adjusted.Pvalue))<0.05,],file=paste(sep="",outPrefix,"_DE_pairwise_residuals-VP_Sequencing.Batch-RNA.Harvest.Date_residuals-lmFit_SEX-Age-BMI_IS_Padj0.05_msigDB_C5ALL_Padj0.05.tab"))
}
if(length(diffIrIs)>1){
    msigdbDiffIrIs<-gsea_enrichment(geneList=diffIrIs,backgroundList=bgList,list.gs=msigdb_C5BP_File)
    write.table(sep="\t",row.names=F,col.names=T,quote=F,msigdbDiffIrIs,file=paste(sep="",outPrefix,"_DE_pairwise_residuals-VP_Sequencing.Batch-RNA.Harvest.Date_residuals-lmFit_SEX-Age-BMI_Diff_IRPadj0.05_ISPadj0.1_msigDB_C5BP.tab"))
    msigdbDiffIrIs<-gsea_enrichment(geneList=diffIrIs,backgroundList=bgList,list.gs=msigdb_C5ALL_File)
    write.table(sep="\t",row.names=F,col.names=T,quote=F,msigdbDiffIrIs,file=paste(sep="",outPrefix,"_DE_pairwise_residuals-VP_Sequencing.Batch-RNA.Harvest.Date_residuals-lmFit_SEX-Age-BMI_Diff_IRPadj0.05_ISPadj0.1_msigDB_C5ALL.tab"))
}
if(length(diffIsIr)>1){
    msigdbDiffIsIr<-gsea_enrichment(geneList=diffIsIr,backgroundList=bgList,list.gs=msigdb_C5BP_File)
    write.table(sep="\t",row.names=F,col.names=T,quote=F,msigdbDiffIsIr,file=paste(sep="",outPrefix,"_DE_pairwise_residuals-VP_Sequencing.Batch-RNA.Harvest.Date_residuals-lmFit_SEX-Age-BMI_Diff_ISPadj0.05_IRPadj0.1_msigDB_C5BP.tab"))
    msigdbDiffIsIr<-gsea_enrichment(geneList=diffIsIr,backgroundList=bgList,list.gs=msigdb_C5ALL_File)
    write.table(sep="\t",row.names=F,col.names=T,quote=F,msigdbDiffIsIr,file=paste(sep="",outPrefix,"_DE_pairwise_residuals-VP_Sequencing.Batch-RNA.Harvest.Date_residuals-lmFit_SEX-Age-BMI_Diff_ISPadj0.05_IRPadj0.1_msigDB_C5ALL.tab"))    
}

cat(file=logFile,append=T,"Done.\n\n")

# exit
cat(file=logFile,append=T,"This is the end.\n")	
