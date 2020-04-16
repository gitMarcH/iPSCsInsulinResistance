cat("Starting.\n\n")

# load required libraries and scripts
library(limma)
library(edgeR)
library(WGCNA)
library(coexpp)
library(RColorBrewer)
library(doParallel)

registerDoParallel(12)

coexppScript<-"/hpc/users/henrim02/VRC/scripts/R_scripts/coexppMinervaSpecifiableParameters.R"
source("/hpc/users/henrim02/VRC/scripts/R_scripts/rnaSeqLimmaPlotPCAFunction_general.R")

handoverDirCoexp<-"/sc/orga/projects/changr04a/marc_handOver/GENESIPS/IR_IS/coexpression"

args<-commandArgs(TRUE)
ver<-unlist(strsplit(split=",",args[1])) # v2 or v3 or v2,v3 or v3,v2
iris<-unlist(strsplit(split=",",args[2])) # ir or is or ir,is or is,ir

########
## V2 ##
########

if(is.element(el="v2",set=ver)){
cat("\n\n## V2 ##\n\n")

# set up (v2)
inputDirV2<-"/sc/orga/projects/changr04a/marc_VRC/IR-IS_Biomarker/v6.317.88_June23"

expFile<-paste(sep="",inputDirV2,"/output/Res.KNOWN.VP--BATCH_RNA.method_lmFit--Reprogramming.Source.Cell_SEX_RACE_Age_BMI.v6.317.88.cutoff_cpm1_sample0.3.RData")
infoFile<-paste(sep="",inputDirV2,"/combined_withnewdata_317samples_covariates_IRvsIS_IVAN.txt")

outDir<-paste(sep="/",handoverDirCoexp,"v2_AS/redo_v2_coexp")
#if(!dir.exists(outDir)){cat(paste(sep="","Creating ",outDir,".\n")); dir.create(outDir,recursive=T)}

# read data
cat("Reading data\n")
load(expFile) # object Res.KNOWN.AllKNOWN
newinfo<-read.table(infoFile,sep="\t",header=T,stringsAsFactor=F)

# match covariates and expression residual data, then extract IS, resp. IR samples only
cat("Matching samples and extracting IS, resp. IR samples only.\n")
commonSamples<-intersect(newinfo$Sample.Name,colnames(Res.KNOWN.AllKNOWN))
newinfo<-newinfo[match(commonSamples,newinfo$Sample.Name),]
Res.KNOWN.AllKNOWN<-Res.KNOWN.AllKNOWN[,match(commonSamples,colnames(Res.KNOWN.AllKNOWN))]

isSamples<-newinfo$Sample.Name[newinfo$State=="Insulin sensitive"]
isInfo<-newinfo[match(isSamples,newinfo$Sample.Name),]
isRes<-Res.KNOWN.AllKNOWN[,match(isSamples,colnames(Res.KNOWN.AllKNOWN))]
save(isRes,file=paste(sep="",outDir,"/v2_AS_IS_residuals.RData"))

irSamples<-newinfo$Sample.Name[newinfo$State=="Insulin resistant"]
irInfo<-newinfo[match(irSamples,newinfo$Sample.Name),]
irRes<-Res.KNOWN.AllKNOWN[,match(irSamples,colnames(Res.KNOWN.AllKNOWN))]
save(irRes,file=paste(sep="",outDir,"/v2_AS_IR_residuals.RData"))

# build co-expression networks for several the beta values that were originally used

# v2 IS beta=9
if(is.element(el="is",set=iris)){
beta<-9
cat(paste(sep="","Doing the coexpression network for v2 IS with beta",beta,".\n"))
betaDir<-paste(sep="",outDir,"/coexp_v2_IS_outlierThr150_beta",beta,"_R3.0.3")
#if(!dir.exists(betaDir)){dir.create(recursive=T,betaDir)}
cat(paste(sep=" ","Rscript",coexppScript,betaDir,paste(sep="/",outDir,"v2_AS_IS_residuals.RData"),"GENESIPS_IRIS_v2_IS",paste(sep="","coexp_beta",beta),beta,TRUE,150,"&>",paste(sep="/",betaDir,paste(sep="","GENESIPS_IRIS_v2_IS_coexp_beta",beta,"_R.log\n"))))
system2(command="Rscript",args=c(coexppScript,betaDir,paste(sep="/",outDir,"v2_AS_IS_residuals.RData"),"GENESIPS_IRIS_v2_IS",paste(sep="","coexp_beta",beta),beta,TRUE,150,"&>",paste(sep="/",betaDir,paste(sep="","GENESIPS_IRIS_v2_IS_coexp_beta",beta,"_R.log"))))
}

# v2 IR beta=7.5
if(is.element(el="ir",set=iris)){
beta<-7.5
cat(paste(sep="","Doing the coexpression network for v2 IR with beta",beta,".\n"))
betaDir<-paste(sep="",outDir,"/coexp_v2_IR_outlierThr600_beta",beta,"_R3.0.3")
#if(!dir.exists(betaDir)){dir.create(recursive=T,betaDir)}
cat(paste(sep=" ","Rscript",coexppScript,betaDir,paste(sep="/",outDir,"v2_AS_IR_residuals.RData"),"GENESIPS_IRIS_v2_IR",paste(sep="","coexp_beta",beta),beta,TRUE,"&>",paste(sep="/",betaDir,paste(sep="","GENESIPS_IRIS_v2_IR_coexp_beta",beta,"_R.log\n"))))
system2(command="Rscript",args=c(coexppScript,betaDir,paste(sep="/",outDir,"v2_AS_IR_residuals.RData"),"GENESIPS_IRIS_v2_IR",paste(sep="","coexp_beta",beta),beta,TRUE,"&>",paste(sep="/",betaDir,paste(sep="","GENESIPS_IRIS_v2_IR_coexp_beta",beta,"_R.log"))))
}
}

########
## V3 ##
########

if(is.element(el="v3",set=ver)){
cat("\n\n## V3 ##\n\n")

# set up (v3)
inputDirV3<-"/sc/orga/projects/changr04a/marc_VRC/IR-IS_Biomarker/v6.317.88_June23"

expFile<-paste(sep="",inputDirV3,"/output/Res.KNOWN.VP--BATCH_RNA.method_lmFit--Reprogramming.Source.Cell_SEX_RACE_Age_BMI.v6.317.88.cutoff_cpm1_sample0.3.AvgResPerPatient.RData")
infoFile<-paste(sep="",inputDirV3,"/combined_withnewdata_317samples_covariates_IRvsIS_IVAN.txt")

outDir<-paste(sep="/",handoverDirCoexp,"v3_ApG/redo_v3_coexp")
#if(!dir.exists(outDir)){cat(paste(sep="","Creating ",outDir,".\n")); dir.create(outDir,recursive=T)}

# read data
cat("Reading data\n")
load(expFile) # object Res.KNOWN.AllKNOWN.AvgRes.PATIENT"
newinfo<-read.table(infoFile,sep="\t",header=T,stringsAsFactor=F)
newinfo2<-data.frame(Sample.Name=ifelse(unique(newinfo$PatientID)=="845A","845A",formatC(as.integer(unique(newinfo$PatientID)),width = 3, format = "d", flag = "0")),State=NA,PatientID=unique(newinfo$PatientID))
for(i in 1:nrow(newinfo2)){newinfo2$State[i]<-newinfo$State[newinfo$PatientID==newinfo2$PatientID[i]][1]}
newinfo<-newinfo2
rm(newinfo2)

# match covariates and expression residual data, then extract IS, resp. IR samples only
cat("Matching samples and extracting IS, resp. IR samples only.\n")
commonSamples<-intersect(newinfo$Sample.Name,colnames(Res.KNOWN.AllKNOWN.AvgRes.PATIENT))
newinfo<-newinfo[match(commonSamples,newinfo$Sample.Name),]
Res.KNOWN.AllKNOWN.AvgRes.PATIENT<-Res.KNOWN.AllKNOWN.AvgRes.PATIENT[,match(commonSamples,colnames(Res.KNOWN.AllKNOWN.AvgRes.PATIENT))]

isSamples<-newinfo$Sample.Name[newinfo$State=="Insulin sensitive"]
isInfo<-newinfo[match(isSamples,newinfo$Sample.Name),]
isRes<-Res.KNOWN.AllKNOWN.AvgRes.PATIENT[,match(isSamples,colnames(Res.KNOWN.AllKNOWN.AvgRes.PATIENT))]
save(isRes,file=paste(sep="",outDir,"/v3_AS_IS_residuals.RData"))

irSamples<-newinfo$Sample.Name[newinfo$State=="Insulin resistant"]
irInfo<-newinfo[match(irSamples,newinfo$Sample.Name),]
irRes<-Res.KNOWN.AllKNOWN.AvgRes.PATIENT[,match(irSamples,colnames(Res.KNOWN.AllKNOWN.AvgRes.PATIENT))]
save(irRes,file=paste(sep="",outDir,"/v3_AS_IR_residuals.RData"))

# build co-expression networks for several the beta values that were originally used

# v3 IS beta=9.5
if(is.element(el="is",set=iris)){
beta<-9.5
cat(paste(sep="","Doing the coexpression network for v3 IS with beta",beta,".\n"))
betaDir<-paste(sep="",outDir,"/coexp_v3_IS_outlierThr600_beta",beta,"_R3.0.3")
#if(!dir.exists(betaDir)){dir.create(recursive=T,betaDir)}
cat(paste(sep=" ","Rscript",coexppScript,betaDir,paste(sep="/",outDir,"v3_AS_IS_residuals.RData"),"GENESIPS_IRIS_v3_IS",paste(sep="","coexp_beta",beta),beta,TRUE,"&>",paste(sep="/",betaDir,paste(sep="","GENESIPS_IRIS_v3_IS_coexp_beta",beta,"_R.log\n"))))
system2(command="Rscript",args=c(coexppScript,betaDir,paste(sep="/",outDir,"v3_AS_IS_residuals.RData"),"GENESIPS_IRIS_v3_IS",paste(sep="","coexp_beta",beta),beta,TRUE,"&>",paste(sep="/",betaDir,paste(sep="","GENESIPS_IRIS_v3_IS_coexp_beta",beta,"_R.log"))))
}

# v3 IR beta=7.5
if(is.element(el="ir",set=iris)){
beta<-7.5
cat(paste(sep="","Doing the coexpression network for v3 IR with beta",beta,".\n"))
betaDir<-paste(sep="",outDir,"/coexp_v3_IR_outlierThr600_beta",beta,"_R3.0.3")
#if(!dir.exists(betaDir)){dir.create(recursive=T,betaDir)}
cat(paste(sep=" ","Rscript",coexppScript,betaDir,paste(sep="/",outDir,"v3_AS_IR_residuals.RData"),"GENESIPS_IRIS_v3_IR",paste(sep="","coexp_beta",beta),beta,TRUE,"&>",paste(sep="/",betaDir,paste(sep="","GENESIPS_IRIS_v3_IR_coexp_beta",beta,"_R.log\n"))))
system2(command="Rscript",args=c(coexppScript,betaDir,paste(sep="/",outDir,"v3_AS_IR_residuals.RData"),"GENESIPS_IRIS_v3_IR",paste(sep="","coexp_beta",beta),beta,TRUE,"&>",paste(sep="/",betaDir,paste(sep="","GENESIPS_IRIS_v3_IR_coexp_beta",beta,"_R.log"))))
}
}

# exit
cat("\nThis is the end.\n")
