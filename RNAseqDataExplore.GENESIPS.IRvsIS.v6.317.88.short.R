rm(list = ls())
library(limma)
library(edgeR)
library(ggplot2)
library(reshape2)
library(ape)

workDir<-commandArgs(TRUE)[1]
hpcTRUE<-as.logical(commandArgs(TRUE)[2])

if(!is.na(hpcTRUE) & hpcTRUE==TRUE){
	library(variancePartition,lib.loc="/hpc/users/henrim02/.Rlib")
}else{
	library(variancePartition)
}

cat(paste(sep="","Input prameters:\n\tworkDir = < ",workDir," >\n\thpcTRUE = < ",hpcTRUE," >\n\n"))

if(!is.na(workDir)){
	setwd(workDir)
}else{
	#setwd("/Users/henrim02/work/geneIPS/IR-IS_Biomarker/v6.317.88_June27")
	setwd("/Users/henrim02/work/geneIPS/IR-IS_Biomarker/v6.317.88_June23")
}
system("mkdir -p output")

#For GenesiPSC's looser cutoff, genes with > 0.1 CPM in at least 10% of experiments were retained for the liberal cutoff.
#For GenesiPSC's tight cutoff for Gabriel's analysis, 1 CPM in at least 30% of the experiments were retained 
    MIN_GENE_CPM=1#.1 #original value was 1 according to Menachem. 
    MIN_SAMPLE_PERCENT_WITH_MIN_GENE_CPM=0.3#0.1
    #MIN_GENE_CPM=0.1 #original value was 1 according to Menachem. 
    #MIN_SAMPLE_PERCENT_WITH_MIN_GENE_CPM=0.1
    low_expression_cutoff=paste("cutoff",paste("cpm=",MIN_GENE_CPM,sep=''),paste("sample=",MIN_SAMPLE_PERCENT_WITH_MIN_GENE_CPM,sep=''),sep='_')
    low_expression_cutoff<-gsub(pattern="=",replacement="",low_expression_cutoff)

################################################################################
# Functions:
################################################################################
getGeneFilteredGeneExprMatrix <- function(gene_expression_counts) {
       # if (!is.null(ONLY_USE_GENES)) {
       #     useGenes = colnames(gene_expression_counts)
       #     useGenes = useGenes[useGenes %in% ONLY_USE_GENES]
       #     gene_expression_counts = gene_expression_counts[, useGenes]
       #     writeLines(paste("\nLimiting expression data to ", length(useGenes), " genes specified by the ONLY_USE_GENES parameter.", sep=""))
       # }
        
        # Make edgeR object:
        MATRIX.ALL_GENES = DGEList(counts=gene_expression_counts, genes=rownames(gene_expression_counts))

        # Keep genes with at least MIN_GENE_CPM count-per-million reads (cpm) in at least (MIN_SAMPLE_PERCENT_WITH_MIN_GENE_CPM)% of the samples:
        #Gabriel: version-2: new low-expression gene cutoff: this is due to PCA test on different cutoff to see the batch effect converge
        fracSamplesWithMinCPM = rowSums(cpm(MATRIX.ALL_GENES) >MIN_GENE_CPM)
        isNonLowExpr = fracSamplesWithMinCPM >= MIN_SAMPLE_PERCENT_WITH_MIN_GENE_CPM*ncol(gene_expression_counts)
        
        #Gabriel: version-1.
        #fracSamplesWithMinCPM = rowSums(cpm(MATRIX.ALL_GENES) >1)
        #isNonLowExpr = fracSamplesWithMinCPM >= 50
        #ORG Menachem
        #fracSamplesWithMinCPM = rowMeans(cpm(MATRIX.ALL_GENES) >=MIN_GENE_CPM)
        #isNonLowExpr = fracSamplesWithMinCPM >= MIN_SAMPLE_PERCENT_WITH_MIN_GENE_CPM
        
        
        
        MATRIX.NON_LOW_GENES = MATRIX.ALL_GENES[isNonLowExpr, ]
        writeLines(paste("\nWill normalize expression counts for ", nrow(MATRIX.NON_LOW_GENES), " genes (those with a minimum of ", MIN_GENE_CPM, " CPM in at least ", sprintf("%.2f", 100 * MIN_SAMPLE_PERCENT_WITH_MIN_GENE_CPM), "% of the ", ncol(MATRIX.NON_LOW_GENES), " samples).", sep=""))


        # FRACTION_BIN_WIDTH = 0.02
        # plotFracSamplesWithMinCPM = data.frame(GeneFeature=names(fracSamplesWithMinCPM), fracSamplesWithMinCPM=as.numeric(fracSamplesWithMinCPM))
        # gRes = ggplot(plotFracSamplesWithMinCPM, aes(x=fracSamplesWithMinCPM))
        # gRes = gRes + geom_vline(xintercept=MIN_SAMPLE_PERCENT_WITH_MIN_GENE_CPM, linetype="solid", col="red")
        # gRes = gRes + geom_histogram(color="black", fill="white", binwidth=FRACTION_BIN_WIDTH) #+ scale_x_log10()
        # gRes = gRes + xlab(paste("Fraction of samples with at least ", MIN_GENE_CPM, " CPM", sep="")) + ylab("# of genes")

        return(list(filteredExprMatrix=MATRIX.NON_LOW_GENES, plotHist=NULL))
    }



calcResiduals <- function(geneBySampleValues, samplesByCovariates, varsToAddBackIn=NULL, sampleWeights=NULL) {
        #################################################################################
        # Use the calcResiduals() code after this section in a for loop:
        #################################################################################
        if (is.matrix(sampleWeights)) {
            residualizedMat = matrix(NA, nrow=nrow(geneBySampleValues), ncol=ncol(geneBySampleValues), dimnames=dimnames(geneBySampleValues))
            for (gInd in 1:nrow(geneBySampleValues)) {
                gRow = calcResiduals(geneBySampleValues[gInd, , drop=FALSE], samplesByCovariates, varsToAddBackIn, sampleWeights[gInd, ])
                residualizedMat[gInd, ] = gRow
            }
            return(residualizedMat)
        }
        #################################################################################
        
        #result.lm = lsfit(x=samplesByCovariates, y=t(geneBySampleValues), wt=sampleWeights, intercept=FALSE)
        
        # Formula of "y ~ 0 + x" means no intercept:
        result.lm = lm(t(geneBySampleValues) ~ 0 + samplesByCovariates, weights=sampleWeights)
        covarNames = colnames(samplesByCovariates)
        
        coef = result.lm$coefficients
        isMatrixForm = is.matrix(coef)
        if (isMatrixForm) {
            rownames(coef) = covarNames
        }
        else {
            names(coef) = covarNames
        }

        allVarsToAddBack = '(Intercept)' # not sure why this is hardcoded: the fitting formula does not fit an intercept; also even if it did, taking the intersection with covarNames would get rid of the intercept anyway
        if (!is.null(varsToAddBackIn)) {
            allVarsToAddBack = c(allVarsToAddBack, varsToAddBackIn)
        }
        allVarsToAddBack = intersect(allVarsToAddBack, covarNames)
        
        residualizedMat = result.lm$residuals
        for (v in allVarsToAddBack) {
            if (isMatrixForm) {
                multCoef = coef[v, , drop=FALSE]
            }
            else {
                multCoef = coef[v]
            }
            residualizedMat = residualizedMat + samplesByCovariates[, v, drop=FALSE] %*% multCoef
        }

        residualizedMat = t(residualizedMat)
        rownames(residualizedMat) = rownames(geneBySampleValues)
        colnames(residualizedMat) = colnames(geneBySampleValues)

        return(residualizedMat)
    }

#############################################################################################################
# Analysis: 3rd round: removed outlier samples (THIS IS AFTER OUTLIER DETECTION, USE THIS FOR THE PAPER!!!!)#
# For IR-IS paper, we use 317 samples samples (3rd round data-final freeze)
#############################################################################################################
##new data: (version-3: This data consist of 3rd batch data+299 samples=317 sample)
dataversion="v6.317.88"
rawcount<-read.table("GENESIPS.combined.withnewdata.317samples.txt", header = TRUE,sep="\t")
rawcountHoldout<-read.table("GENESIPS.IRIS.holdout.geneCounts.csv",header=T,sep=",")
#outlier removal:
#845_1, 845_6, 845_10, 864_2, 558_1, 558_2, 541_5
outlier<-cbind("864_2","845_1","845_10","845_6","558_1","558_2","541_5")

# 317 samples
genenames<-rawcount[,1]
sample0<-rawcount[,-1]
rownames(sample0)<-genenames
colnames(sample0) = gsub("^X", "", colnames(sample0))
#remove the outliers
idx = match(as.character(outlier),colnames(sample0))
sample<-sample0[,-idx]
#write.table(rownames(sample),quote=FALSE,sep="\t",file=paste("output/total.gene.ensemble",dataversion,"txt",sep="."),col.names=F,row.names=F)

# 88 holdout samples (85 after removing patient 102 who has also been included in the 317 samples)
genenamesHO<-rawcountHoldout[,1]
sample0HO<-rawcountHoldout[,-1]
rownames(sample0HO)<-genenamesHO
colnames(sample0HO) = gsub("^X", "", colnames(sample0HO))
sampleHoldout<-sample0HO[,gsub(pattern="_[0-9]*",replacement="",colnames(sample0HO))!="102"]
samples2Remove.Holdout<-c("148_01","148_03") # clones frm patient 148 appear in 2 batches; we keep only 1 batch of these
sampleHoldout<-sampleHoldout[,!is.element(el=colnames(sampleHoldout),set=samples2Remove.Holdout)]
rm(sample0HO)

# joint sample matrix
commonGenes<-intersect(rownames(sample),rownames(sampleHoldout))
sampleAll<-cbind(sample[match(commonGenes,rownames(sampleHoldout)),],sampleHoldout[match(commonGenes,rownames(sampleHoldout)),])

#************************************************************
#Remove low-expressed genes
#************************************************************
filteredMatPlot = getGeneFilteredGeneExprMatrix(sample)

TRUE.GENE_EXPRESSION_DGELIST_MAT = filteredMatPlot$filteredExprMatrix
write.table(rownames(TRUE.GENE_EXPRESSION_DGELIST_MAT),quote=FALSE,sep="\t",file=paste("output/expressed.gene.ensemble","CPMcut=",MIN_GENE_CPM,"Percent=",MIN_SAMPLE_PERCENT_WITH_MIN_GENE_CPM,dataversion,"txt",sep="."),col.names=F,row.names=F)

TRUE.GENE_EXPRESSION_MAT = TRUE.GENE_EXPRESSION_DGELIST_MAT$counts

#-- extra: average the counts for all clones per patient ID
patIDs<-gsub(pattern="_[0-9]*",replacement="",colnames(TRUE.GENE_EXPRESSION_MAT))
uniqPatIDs<-unique(patIDs)
TRUE.GENE_EXPRESSION_MAT.PATIENT<-matrix(nrow=nrow(TRUE.GENE_EXPRESSION_MAT),ncol=length(uniqPatIDs))
rownames(TRUE.GENE_EXPRESSION_MAT.PATIENT)<-rownames(TRUE.GENE_EXPRESSION_MAT)
colnames(TRUE.GENE_EXPRESSION_MAT.PATIENT)<-uniqPatIDs
for(patID in uniqPatIDs){
	if(sum(patIDs==patID)>1){
		TRUE.GENE_EXPRESSION_MAT.PATIENT[,colnames(TRUE.GENE_EXPRESSION_MAT.PATIENT)==patID]<-rowSums(TRUE.GENE_EXPRESSION_MAT[,patIDs==patID]) # originally we used rowMeans; but rowSums should be betetr for subsequent variance modelling by voom etc -- after all for a patient with 4 clones we will have many more counts (i.e. information) than for a patient with just 1 clone
	}else{
		TRUE.GENE_EXPRESSION_MAT.PATIENT[,colnames(TRUE.GENE_EXPRESSION_MAT.PATIENT)==patID]<-TRUE.GENE_EXPRESSION_MAT[,patIDs==patID]
	}
}
TRUE.GENE_EXPRESSION_DGELIST_MAT.PATIENT<-DGEList(counts=TRUE.GENE_EXPRESSION_MAT.PATIENT, genes=rownames(TRUE.GENE_EXPRESSION_MAT.PATIENT))

#------------NEW-Nov-2015: Just to double check the effect of gene length on gene expression variance------------#
# RPKM: normalize against gene length-----------#
#write.table(rownames(TRUE.GENE_EXPRESSION_MAT),quote=FALSE,sep="\n",file=paste("output/filtered.gene.ensemble",dataversion,"txt",sep="."),col.names=F,row.names=F)
#rpkm(TRUE.GENE_EXPRESSION_MAT,)

#------------# TMM normalization for effective library size-----------------#
TRUE.GENE_EXPRESSION_DGELIST_MAT.NORM <-calcNormFactors(TRUE.GENE_EXPRESSION_DGELIST_MAT,method='TMM')

#-- extra: do the same for the average count matrix
TRUE.GENE_EXPRESSION_DGELIST_MAT.PATIENT.NORM <-calcNormFactors(TRUE.GENE_EXPRESSION_DGELIST_MAT.PATIENT,method='TMM')

# same thing for holdout data
filteredMatPlotHO = getGeneFilteredGeneExprMatrix(sampleHoldout)
TRUE.GENE_EXPRESSION_DGELIST_MAT_Holdout = filteredMatPlotHO$filteredExprMatrix
TRUE.GENE_EXPRESSION_MAT_Holdout = TRUE.GENE_EXPRESSION_DGELIST_MAT_Holdout$counts
TRUE.GENE_EXPRESSION_DGELIST_MAT_Holdout.NORM <-calcNormFactors(TRUE.GENE_EXPRESSION_DGELIST_MAT_Holdout,method='TMM')

patIDHOs<-gsub(pattern="_[0-9]*",replacement="",colnames(TRUE.GENE_EXPRESSION_MAT_Holdout))
uniqPatIDHOs<-unique(patIDHOs)
TRUE.GENE_EXPRESSION_MAT_Holdout.PATIENT<-matrix(nrow=nrow(TRUE.GENE_EXPRESSION_MAT_Holdout),ncol=length(uniqPatIDHOs))
rownames(TRUE.GENE_EXPRESSION_MAT_Holdout.PATIENT)<-rownames(TRUE.GENE_EXPRESSION_MAT_Holdout)
colnames(TRUE.GENE_EXPRESSION_MAT_Holdout.PATIENT)<-uniqPatIDHOs
for(patIDHO in uniqPatIDHOs){
	if(sum(patIDHOs==patIDHO)>1){
		TRUE.GENE_EXPRESSION_MAT_Holdout.PATIENT[,colnames(TRUE.GENE_EXPRESSION_MAT_Holdout.PATIENT)==patIDHO]<-rowMeans(TRUE.GENE_EXPRESSION_MAT_Holdout[,patIDHOs==patIDHO])
	}else{
		TRUE.GENE_EXPRESSION_MAT_Holdout.PATIENT[,colnames(TRUE.GENE_EXPRESSION_MAT_Holdout.PATIENT)==patIDHO]<-TRUE.GENE_EXPRESSION_MAT_Holdout[,patIDHOs==patIDHO]
	}
}
TRUE.GENE_EXPRESSION_DGELIST_MAT_Holdout.PATIENT<-DGEList(counts=TRUE.GENE_EXPRESSION_MAT_Holdout.PATIENT, genes=rownames(TRUE.GENE_EXPRESSION_MAT_Holdout.PATIENT))
TRUE.GENE_EXPRESSION_DGELIST_MAT_Holdout.PATIENT.NORM <-calcNormFactors(TRUE.GENE_EXPRESSION_DGELIST_MAT_Holdout.PATIENT,method='TMM')

# same thing for combined data
filteredMatPlotAll = getGeneFilteredGeneExprMatrix(sampleAll)
TRUE.GENE_EXPRESSION_DGELIST_MAT_All = filteredMatPlotAll$filteredExprMatrix
TRUE.GENE_EXPRESSION_MAT_All = TRUE.GENE_EXPRESSION_DGELIST_MAT_All$counts
TRUE.GENE_EXPRESSION_DGELIST_MAT_All.NORM <-calcNormFactors(TRUE.GENE_EXPRESSION_DGELIST_MAT_All,method='TMM')

patIDAlls<-gsub(pattern="_[0-9]*",replacement="",colnames(TRUE.GENE_EXPRESSION_MAT_All))
uniqPatIDAlls<-unique(patIDAlls)
TRUE.GENE_EXPRESSION_MAT_All.PATIENT<-matrix(nrow=nrow(TRUE.GENE_EXPRESSION_MAT_All),ncol=length(uniqPatIDAlls))
rownames(TRUE.GENE_EXPRESSION_MAT_All.PATIENT)<-rownames(TRUE.GENE_EXPRESSION_MAT_All)
colnames(TRUE.GENE_EXPRESSION_MAT_All.PATIENT)<-uniqPatIDAlls
for(patIDAll in uniqPatIDAlls){
	if(sum(patIDAlls==patIDAll)>1){
		TRUE.GENE_EXPRESSION_MAT_All.PATIENT[,colnames(TRUE.GENE_EXPRESSION_MAT_All.PATIENT)==patIDAll]<-rowMeans(TRUE.GENE_EXPRESSION_MAT_All[,patIDAlls==patIDAll])
	}else{
		TRUE.GENE_EXPRESSION_MAT_All.PATIENT[,colnames(TRUE.GENE_EXPRESSION_MAT_All.PATIENT)==patIDAll]<-TRUE.GENE_EXPRESSION_MAT_All[,patIDAlls==patIDAll]
	}
}
TRUE.GENE_EXPRESSION_DGELIST_MAT_All.PATIENT<-DGEList(counts=TRUE.GENE_EXPRESSION_MAT_All.PATIENT, genes=rownames(TRUE.GENE_EXPRESSION_MAT_All.PATIENT))
TRUE.GENE_EXPRESSION_DGELIST_MAT_All.PATIENT.NORM <-calcNormFactors(TRUE.GENE_EXPRESSION_DGELIST_MAT_All.PATIENT,method='TMM')

#****************Calculate log-CPM by voom****************
# voom normalization
#*********************************************************
pdf(paste("output/voom.OutlierRemove",dataversion,"pdf",sep="."))
vobj = voom(TRUE.GENE_EXPRESSION_DGELIST_MAT.NORM, plot=T)
VOOM_RAW_LOG_EXPRESSION_MAT = vobj$E
dev.off()
write.table(VOOM_RAW_LOG_EXPRESSION_MAT,quote=FALSE,sep="\t",col.names=TRUE,file="output/Voom_Log2CPM_noCov_OutlierRemoval.txt")

#-- extra: do the same for the average count matrix
vobj.PATIENT = voom(TRUE.GENE_EXPRESSION_DGELIST_MAT.PATIENT.NORM, plot=F)
VOOM_RAW_LOG_EXPRESSION_MAT.PATIENT = vobj.PATIENT$E

# Plot distribtuion of CPM values
pdf(paste("output/distribution.OutlierRemoval.cpm",dataversion,"pdf",sep="."))
plot(density(vobj$E[,1]))
for(i in 2:ncol(vobj)){
	lines(density(vobj$E[,i]))
}
dev.off()

# do the same for holdout data
pdf(paste("output/voom.OutlierRemove",dataversion,"Holdout.pdf",sep="."))
vobjHO = voom(TRUE.GENE_EXPRESSION_DGELIST_MAT_Holdout.NORM, plot=T)
VOOM_RAW_LOG_EXPRESSION_MAT_Holdout = vobjHO$E
dev.off()
write.table(VOOM_RAW_LOG_EXPRESSION_MAT,quote=FALSE,sep="\t",col.names=TRUE,file="output/Voom_Log2CPM_noCov_OutlierRemoval.Holdout.txt")
# avg count data (holdout0)
vobjHO.PATIENT = voom(TRUE.GENE_EXPRESSION_DGELIST_MAT_Holdout.PATIENT.NORM, plot=F)
VOOM_RAW_LOG_EXPRESSION_MAT_Holdout.PATIENT = vobjHO.PATIENT$E

# do the same for all data
pdf(paste("output/voom.OutlierRemove",dataversion,"All.pdf",sep="."))
vobjAll = voom(TRUE.GENE_EXPRESSION_DGELIST_MAT_All.NORM, plot=T)
VOOM_RAW_LOG_EXPRESSION_MAT_All = vobjAll$E
dev.off()
write.table(VOOM_RAW_LOG_EXPRESSION_MAT,quote=FALSE,sep="\t",col.names=TRUE,file="output/Voom_Log2CPM_noCov_OutlierRemoval.All.txt")
# avg count data (combined)
vobjAll.PATIENT = voom(TRUE.GENE_EXPRESSION_DGELIST_MAT_All.PATIENT.NORM, plot=F)
VOOM_RAW_LOG_EXPRESSION_MAT_All.PATIENT = vobjAll.PATIENT$E

# # same thing for all (combined) data
# NB voom & library size normalisations are equivalent whether matrices are combined before or after the voom command (up to numerical rounding errors and very small changes due to the slightly different normalisation factors of the gene lists used for the indivdual and combined DGELists objects)
# commonGenes.Tmp<-intersect(rownames(TRUE.GENE_EXPRESSION_DGELIST_MAT$counts),rownames(TRUE.GENE_EXPRESSION_DGELIST_MAT_Holdout$counts))
# TRUE.GENE_EXPRESSION_DGELIST_MAT_All.NORM.Alt<-list(counts=cbind(TRUE.GENE_EXPRESSION_DGELIST_MAT$counts[match(commonGenes.Tmp,rownames(TRUE.GENE_EXPRESSION_DGELIST_MAT$counts)),],TRUE.GENE_EXPRESSION_DGELIST_MAT_Holdout.NORM$counts[match(commonGenes.Tmp,rownames(TRUE.GENE_EXPRESSION_DGELIST_MAT_Holdout.NORM$counts)),]),samples=data.frame(group=factor(rep(1,ncol(TRUE.GENE_EXPRESSION_DGELIST_MAT.NORM$counts)+ncol(TRUE.GENE_EXPRESSION_DGELIST_MAT_Holdout.NORM$counts))),lib.size=colSums(cbind(TRUE.GENE_EXPRESSION_DGELIST_MAT$counts[match(commonGenes.Tmp,rownames(TRUE.GENE_EXPRESSION_DGELIST_MAT$counts)),],TRUE.GENE_EXPRESSION_DGELIST_MAT_Holdout.NORM$counts[match(commonGenes.Tmp,rownames(TRUE.GENE_EXPRESSION_DGELIST_MAT_Holdout.NORM$counts)),])),norm.factors=c(TRUE.GENE_EXPRESSION_DGELIST_MAT.NORM$samples$norm.factors,TRUE.GENE_EXPRESSION_DGELIST_MAT_Holdout.NORM$samples$norm.factors)),genes=list(commonGenes.Tmp))
# TRUE.GENE_EXPRESSION_DGELIST_MAT_All.NORM.Alt<-as(TRUE.GENE_EXPRESSION_DGELIST_MAT_All.NORM.Alt,Class="DGEList")
# vobjAll.Alt = voom(TRUE.GENE_EXPRESSION_DGELIST_MAT_All.NORM.Alt, plot=T)
# VOOM_RAW_LOG_EXPRESSION_MAT_All.Alt = vobjAll.Alt$E

#*****************************************************************************#
#Before any adjustment, we use PCA to determine which covariates are relavent
#*****************************************************************************#
##new data: (version-3: This data consist of 3rd batch data+299 samples=317 sample)
#Load in RDS object; USE ONLY ONCE
#info3.rds=readRDS(file="info3.RDS")
#write.table(info3.rds,quote=FALSE,sep="\t",col.names=TRUE,file="output/combined_withnewdata_317samples_covariates_IRvsIS.txt.txt")
info = read.table("combined_withnewdata_317samples_covariates_IRvsIS_IVAN.txt", header=TRUE, sep="\t")
info = info[match(colnames(sample), as.character(info$Sample.Name)),]
tmpPatIDs<-as.character(info$PatientID)
for(i in 1:length(tmpPatIDs)){if(length(unlist(strsplit(split="",tmpPatIDs[i])))==2){tmpPatIDs[i]<-paste(sep="","0",tmpPatIDs[i])}}
info$PatientID<-factor(tmpPatIDs)
info$RACE[info$RACE=="Black"]<-"African.American" # consolidate 'Black' & 'African.American' race labels
info$RACE[is.element(el=info$PatientID,set=c("054","912","212","285","729"))]<-"East.Asian" # correct IDs of patient recorded as simply Asian"
info$RACE<-factor(info$RACE)

info.Holdout<-read.table("newSamples_iPSC_personalAndTechnicalCovariates.tab",sep="\t",stringsAsFactors=F,header=T)
info.Holdout<-info.Holdout[match(colnames(sampleHoldout),info.Holdout$Sample.Name),]
tmpPatIDs<-as.character(info.Holdout$PatientID)
for(i in 1:length(tmpPatIDs)){if(length(unlist(strsplit(split="",tmpPatIDs[i])))==2){tmpPatIDs[i]<-paste(sep="","0",tmpPatIDs[i])}}
info.Holdout$PatientID<-factor(tmpPatIDs)
info.Holdout$RACE<-as.character(info.Holdout$RACE)
info.Holdout$RACE[info.Holdout$RACE=="Black" | info.Holdout$RACE=="Black or African American"]<-"African.American" # consolidate redundant labels
info.Holdout$RACE[info.Holdout$RACE=="East Asian"]<-"East.Asian" # use same labels as the info for the 317 samples
info.Holdout$RACE[info.Holdout$RACE=="South Asian"]<-"South.Asian" # use same labels as the info for the 317 samples
info.Holdout$RACE[info.Holdout$RACE=="White Hispanic"]<-"White.Hispanic" # use same labels as the info for the 317 samples
info.Holdout$RACE[is.element(el=info.Holdout$PatientID,set=c("054","912","212","285","729"))]<-"East.Asian" # correct IDs of patient recorded as simply Asian"
info.Holdout$RACE<-factor(info.Holdout$RACE)

commonCovs<-intersect(colnames(info),colnames(info.Holdout))
info.All<-rbind(info[,match(commonCovs,colnames(info))],info.Holdout[,match(commonCovs,colnames(info.Holdout))])
datOrigin<-c(rep("original317",nrow(info)),rep("holdout85",nrow(info.Holdout)))
info.All<-data.frame(info.All,datOrigin=datOrigin)
info.All<-info.All[match(colnames(sampleAll),info.All$Sample.Name),]

##new data: (version-2: This data consist of new batch data+299 samples=317 sample)
espluripotentmarker<-read.table("ES_PLUR_MARKER_new_tier1_ensemble.txt", header = FALSE,sep="\t")
esdiffmarker<-read.table("ES_PLUR_MARKER_new_tier2_ensemble.txt", header = FALSE,sep="\t")

ES.marker.matrix<-as.matrix(espluripotentmarker)
DIF.marker.matrix<-as.matrix(esdiffmarker)

ES_marker_index.voom<-vector()
DIFF_marker_index.voom<-vector()
for (index in 1:length(as.matrix(espluripotentmarker))){
	ES.marker.name<-ES.marker.matrix[index,]
	#tmp.index<-pmatch(ES.marker.name,rownames(VOOM_WEIGHTED_RESIDUALIZED_WITH_PRIMARY_MAT),nomatch=-1)
	tmp.index<-pmatch(ES.marker.name,rownames(VOOM_RAW_LOG_EXPRESSION_MAT),nomatch=-1)
	if(tmp.index>-1){
		ES_marker_index.voom<-c(ES_marker_index.voom, tmp.index)
	}
}


for (index in 1:length(as.matrix(esdiffmarker))){
	DIFF.marker.name<-DIF.marker.matrix[index,]
	#tmp.index<-pmatch(DIFF.marker.name,rownames(VOOM_WEIGHTED_RESIDUALIZED_WITH_PRIMARY_MAT),nomatch=-1)
	tmp.index<-pmatch(DIFF.marker.name,rownames(VOOM_RAW_LOG_EXPRESSION_MAT),nomatch=-1)

	if(tmp.index>-1){
		DIFF_marker_index.voom<-c(DIFF_marker_index.voom, tmp.index)	
	}
}

ES.marker.sample.voom<-VOOM_RAW_LOG_EXPRESSION_MAT[ES_marker_index.voom,]
DIFF.marker.sample.voom<-VOOM_RAW_LOG_EXPRESSION_MAT[DIFF_marker_index.voom,]
marker.sample.voom<- rbind(ES.marker.sample.voom, DIFF.marker.sample.voom)
#***************************************#
#Plot PCA with ES marker outlier removal#
#***************************************#
pdf(paste("output/PCA_ES.marker_explore_outlierremoval",dataversion,"pdf",sep="."))
SampleByVariable=t(marker.sample.voom)
clonename<-rownames(SampleByVariable)
pca <- prcomp(SampleByVariable, scale=T)
#=======pca-1 vs pca-2=======
ggplot(data.frame(pca$x), aes(x= pca$x[,1], y= pca$x[,2], colour=factor(info$Reprogramming.Source.Cell), label=clonename))+
  scale_colour_discrete(guide =FALSE) +geom_point() +geom_text(aes(label=clonename),hjust=0, vjust=0,size=2)+labs(title="Reprogramming.Source.Cell")

ggplot(data.frame(pca$x), aes(x= pca$x[,1], y= pca$x[,2], colour=factor(info$BATCH), label=clonename))+
  scale_colour_discrete(guide = FALSE) +geom_point() +geom_text(aes(label=clonename),hjust=0, vjust=0,size=2)+labs(title="BATCH")

ggplot(data.frame(pca$x), aes(x= pca$x[,1], y= pca$x[,2], colour=factor(info$Sample.Name), label=clonename))+
  scale_colour_discrete(guide = FALSE) +geom_point() +geom_text(aes(label=clonename),hjust=0, vjust=0,size=2)+labs(title="Sample.Name")
  
ggplot(data.frame(pca$x), aes(x= pca$x[,1], y= pca$x[,2], colour=factor(info$Reprogramming.Batch), label=clonename))+
  scale_colour_discrete(guide = FALSE) +geom_point() +geom_text(aes(label=clonename),hjust=0, vjust=0,size=2) +labs(title="Reprogramming.Batch")
  
ggplot(data.frame(pca$x), aes(x= pca$x[,1], y= pca$x[,2], colour=factor(info$Sendai.Virus.Lot), label=clonename))+
  scale_colour_discrete(guide = FALSE) +geom_point() +geom_text(aes(label=clonename),hjust=0, vjust=0,size=2) +labs(title="Sendai.Virus.Lot")  

ggplot(data.frame(pca$x), aes(x= pca$x[,1], y= pca$x[,2], colour=factor(info$SEX), label=clonename))+
  scale_colour_discrete(guide = FALSE) +geom_point() +geom_text(aes(label=clonename),hjust=0, vjust=0,size=2) +labs(title="SEX") 
  
  
   ggplot(data.frame(pca$x), aes(x= pca$x[,1], y= pca$x[,2], colour=factor(info$RACE), label=clonename))+
  scale_colour_discrete(guide = FALSE) +geom_point() +geom_text(aes(label=clonename),hjust=0, vjust=0,size=2)+labs(title="RACE") 
  
     ggplot(data.frame(pca$x), aes(x= pca$x[,1], y= pca$x[,2], colour=factor(info$State), label=clonename))+scale_colour_discrete(guide = FALSE) +geom_point() +geom_text(aes(label=clonename),hjust=0, vjust=0,size=2) +labs(title="State") 
  
  ggplot(data.frame(pca$x), aes(x= pca$x[,1], y= pca$x[,2], colour=factor(info$RNA.method), label=clonename))+
  scale_colour_discrete(guide = FALSE) +geom_point() +geom_text(aes(label=clonename),hjust=0, vjust=0,size=2) +labs(title="RNA.method") 
  
    ggplot(data.frame(pca$x), aes(x= pca$x[,1], y= pca$x[,2], colour=factor(info$Age), label=clonename))+
  scale_colour_discrete(guide = FALSE) +geom_point() +geom_text(aes(label=clonename),hjust=0, vjust=0,size=2)+labs(title="Age")  
  #=========pca-1 vs pca-3======#
  ggplot(data.frame(pca$x), aes(x= pca$x[,1], y= pca$x[,3], colour=factor(info$Reprogramming.Source.Cell), label=clonename))+
  scale_colour_discrete(guide =FALSE) +geom_point() +geom_text(aes(label=clonename),hjust=0, vjust=0,size=2)+labs(title="Reprogramming.Source.Cell")

ggplot(data.frame(pca$x), aes(x= pca$x[,1], y= pca$x[,3], colour=factor(info$BATCH), label=clonename))+
  scale_colour_discrete(guide = FALSE) +geom_point() +geom_text(aes(label=clonename),hjust=0, vjust=0,size=2)+labs(title="BATCH")

ggplot(data.frame(pca$x), aes(x= pca$x[,1], y= pca$x[,3], colour=factor(info$Sample.Name), label=clonename))+
  scale_colour_discrete(guide = FALSE) +geom_point() +geom_text(aes(label=clonename),hjust=0, vjust=0,size=2)+labs(title="Sample.Name")
  
ggplot(data.frame(pca$x), aes(x= pca$x[,1], y= pca$x[,3], colour=factor(info$Reprogramming.Batch), label=clonename))+
  scale_colour_discrete(guide = FALSE) +geom_point() +geom_text(aes(label=clonename),hjust=0, vjust=0,size=2) +labs(title="Reprogramming.Batch")
  
ggplot(data.frame(pca$x), aes(x= pca$x[,1], y= pca$x[,3], colour=factor(info$Sendai.Virus.Lot), label=clonename))+
  scale_colour_discrete(guide = FALSE) +geom_point() +geom_text(aes(label=clonename),hjust=0, vjust=0,size=2) +labs(title="Sendai.Virus.Lot")  

ggplot(data.frame(pca$x), aes(x= pca$x[,1], y= pca$x[,3], colour=factor(info$SEX), label=clonename))+
  scale_colour_discrete(guide = FALSE) +geom_point() +geom_text(aes(label=clonename),hjust=0, vjust=0,size=2) +labs(title="SEX") 
  
  
   ggplot(data.frame(pca$x), aes(x= pca$x[,1], y= pca$x[,3], colour=factor(info$RACE), label=clonename))+
  scale_colour_discrete(guide = FALSE) +geom_point() +geom_text(aes(label=clonename),hjust=0, vjust=0,size=2)+labs(title="RACE") 
  
     ggplot(data.frame(pca$x), aes(x= pca$x[,1], y= pca$x[,3], colour=factor(info$State), label=clonename))+scale_colour_discrete(guide = FALSE) +geom_point() +geom_text(aes(label=clonename),hjust=0, vjust=0,size=2) +labs(title="State") 
  
  ggplot(data.frame(pca$x), aes(x= pca$x[,1], y= pca$x[,3], colour=factor(info$RNA.method), label=clonename))+
  scale_colour_discrete(guide = FALSE) +geom_point() +geom_text(aes(label=clonename),hjust=0, vjust=0,size=2) +labs(title="RNA.method") 
  
    ggplot(data.frame(pca$x), aes(x= pca$x[,1], y= pca$x[,3], colour=factor(info$Age), label=clonename))+
  scale_colour_discrete(guide = FALSE) +geom_point() +geom_text(aes(label=clonename),hjust=0, vjust=0,size=2)+labs(title="Age") 
  
  #=======pca-2 vs pca-3===========#
   ggplot(data.frame(pca$x), aes(x= pca$x[,2], y= pca$x[,3], colour=factor(info$Reprogramming.Source.Cell), label=clonename))+
  scale_colour_discrete(guide =FALSE) +geom_point() +geom_text(aes(label=clonename),hjust=0, vjust=0,size=2)+labs(title="Reprogramming.Source.Cell")

ggplot(data.frame(pca$x), aes(x= pca$x[,2], y= pca$x[,3], colour=factor(info$BATCH), label=clonename))+
  scale_colour_discrete(guide = FALSE) +geom_point() +geom_text(aes(label=clonename),hjust=0, vjust=0,size=2)+labs(title="BATCH")

ggplot(data.frame(pca$x), aes(x= pca$x[,2], y= pca$x[,3], colour=factor(info$Sample.Name), label=clonename))+
  scale_colour_discrete(guide = FALSE) +geom_point() +geom_text(aes(label=clonename),hjust=0, vjust=0,size=2)+labs(title="Sample.Name")
  
ggplot(data.frame(pca$x), aes(x= pca$x[,2], y= pca$x[,3], colour=factor(info$Reprogramming.Batch), label=clonename))+
  scale_colour_discrete(guide = FALSE) +geom_point() +geom_text(aes(label=clonename),hjust=0, vjust=0,size=2) +labs(title="Reprogramming.Batch")
  
ggplot(data.frame(pca$x), aes(x= pca$x[,2], y= pca$x[,3], colour=factor(info$Sendai.Virus.Lot), label=clonename))+
  scale_colour_discrete(guide = FALSE) +geom_point() +geom_text(aes(label=clonename),hjust=0, vjust=0,size=2) +labs(title="Sendai.Virus.Lot")  

ggplot(data.frame(pca$x), aes(x= pca$x[,2], y= pca$x[,3], colour=factor(info$SEX), label=clonename))+
  scale_colour_discrete(guide = FALSE) +geom_point() +geom_text(aes(label=clonename),hjust=0, vjust=0,size=2) +labs(title="SEX") 
  
  
   ggplot(data.frame(pca$x), aes(x= pca$x[,2], y= pca$x[,3], colour=factor(info$RACE), label=clonename))+
  scale_colour_discrete(guide = FALSE) +geom_point() +geom_text(aes(label=clonename),hjust=0, vjust=0,size=2)+labs(title="RACE") 
  
     ggplot(data.frame(pca$x), aes(x= pca$x[,2], y= pca$x[,3], colour=factor(info$State), label=clonename))+scale_colour_discrete(guide = FALSE) +geom_point() +geom_text(aes(label=clonename),hjust=0, vjust=0,size=2) +labs(title="State") 
  
  ggplot(data.frame(pca$x), aes(x= pca$x[,2], y= pca$x[,3], colour=factor(info$RNA.method), label=clonename))+
  scale_colour_discrete(guide = FALSE) +geom_point() +geom_text(aes(label=clonename),hjust=0, vjust=0,size=2) +labs(title="RNA.method") 
  
    ggplot(data.frame(pca$x), aes(x= pca$x[,2], y= pca$x[,3], colour=factor(info$Age), label=clonename))+
  scale_colour_discrete(guide = FALSE) +geom_point() +geom_text(aes(label=clonename),hjust=0, vjust=0,size=2)+labs(title="Age") 
  
dev.off()


for(datSet in c("train","test","combined")){
	if(datSet=="train"){
		pdf(paste("output/PCA_All.gene_explore_outlierremoval",dataversion,"pdf",sep="."))
		SampleByVariable=t(VOOM_RAW_LOG_EXPRESSION_MAT)
		infoTmp<-info
	}else if(datSet=="test"){
		pdf(paste("output/PCA_All.gene_explore_outlierremoval",dataversion,"Holdout.pdf",sep="."))
		SampleByVariable=t(VOOM_RAW_LOG_EXPRESSION_MAT_Holdout)
		infoTmp<-info.Holdout
	}else if(datSet=="combined"){
		pdf(paste("output/PCA_All.gene_explore_outlierremoval",dataversion,"All.pdf",sep="."))
		SampleByVariable=t(VOOM_RAW_LOG_EXPRESSION_MAT_All)
		infoTmp<-info.All
	}
	print(datSet)
	clonename<-rownames(SampleByVariable)
	pca <- prcomp(SampleByVariable, scale=T)
	
	#=======pca-1 vs pca-2=======
	print(ggplot(data.frame(pca$x), aes(x= pca$x[,1], y= pca$x[,2], colour=factor(infoTmp$Reprogramming.Source.Cell), label=clonename))+
		  scale_colour_discrete(guide =FALSE) +geom_point() +geom_text(aes(label=clonename),hjust=0, vjust=0,size=2)+labs(title="Reprogramming.Source.Cell")+ theme(legend.position=c(0,0),legend.justification=c(0,0),legend.box.just="top",legend.background=element_rect(fill=alpha("white",0.4))) )
		
	print(ggplot(data.frame(pca$x), aes(x= pca$x[,1], y= pca$x[,2], colour=factor(infoTmp$BATCH), label=clonename))+
		  scale_colour_discrete(guide = FALSE) +geom_point() +geom_text(aes(label=clonename),hjust=0, vjust=0,size=2)+labs(title="BATCH")+ theme(legend.position=c(0,1),legend.justification=c(0,1),legend.box.just="top",legend.background=element_rect(fill=alpha("white",0.4))) )
		
	print(ggplot(data.frame(pca$x), aes(x= pca$x[,1], y= pca$x[,2], colour=factor(infoTmp$Sample.Name), label=clonename))+
		  scale_colour_discrete(guide = FALSE) +geom_point() +geom_text(aes(label=clonename),hjust=0, vjust=0,size=2)+labs(title="Sample.Name")+ theme(legend.position=c(0,1),legend.justification=c(0,1),legend.box.just="top",legend.background=element_rect(fill=alpha("white",0.4))) )
		  
	print(ggplot(data.frame(pca$x), aes(x= pca$x[,1], y= pca$x[,2], colour=factor(infoTmp$Reprogramming.Batch), label=clonename))+
		  scale_colour_discrete(guide = FALSE) +geom_point() +geom_text(aes(label=clonename),hjust=0, vjust=0,size=2) +labs(title="Reprogramming.Batch")+ theme(legend.position=c(0,1),legend.justification=c(0,1),legend.box.just="top",legend.background=element_rect(fill=alpha("white",0.4))) )
		  
	if(datSet=="train"){print(ggplot(data.frame(pca$x), aes(x= pca$x[,1], y= pca$x[,2], colour=factor(infoTmp$Sendai.Virus.Lot), label=clonename))+
		  scale_colour_discrete(guide = FALSE) +geom_point() +geom_text(aes(label=clonename),hjust=0, vjust=0,size=2) +labs(title="Sendai.Virus.Lot") + theme(legend.position=c(0,1),legend.justification=c(0,1),legend.box.just="top",legend.background=element_rect(fill=alpha("white",0.4)))  )}
		
	print(ggplot(data.frame(pca$x), aes(x= pca$x[,1], y= pca$x[,2], colour=factor(infoTmp$SEX), label=clonename))+
		  scale_colour_discrete(guide = FALSE) +geom_point() +geom_text(aes(label=clonename),hjust=0, vjust=0,size=2) +labs(title="SEX") + theme(legend.position=c(0,1),legend.justification=c(0,1),legend.box.just="top",legend.background=element_rect(fill=alpha("white",0.4))) )
		  
	print(ggplot(data.frame(pca$x), aes(x= pca$x[,1], y= pca$x[,2], colour=factor(infoTmp$RACE), label=clonename))+
		  scale_colour_discrete(guide = FALSE) +geom_point() +geom_text(aes(label=clonename),hjust=0, vjust=0,size=2)+labs(title="RACE")+ theme(legend.position=c(0,1),legend.justification=c(0,1),legend.box.just="top",legend.background=element_rect(fill=alpha("white",0.4)))  )
		  
	print(ggplot(data.frame(pca$x), aes(x= pca$x[,1], y= pca$x[,2], colour=factor(infoTmp$State), label=clonename))+
		  scale_colour_discrete(guide = FALSE) +geom_point() +geom_text(aes(label=clonename),hjust=0, vjust=0,size=2) +labs(title="State")+ theme(legend.position=c(0,1),legend.justification=c(0,1),legend.box.just="top",legend.background=element_rect(fill=alpha("white",0.4)))  )
		  
	print(ggplot(data.frame(pca$x), aes(x= pca$x[,1], y= pca$x[,2], colour=factor(infoTmp$RNA.method), label=clonename))+
		  scale_colour_discrete(guide = FALSE) +geom_point() +geom_text(aes(label=clonename),hjust=0, vjust=0,size=2) +labs(title="RNA.method")+ theme(legend.position=c(0,1),legend.justification=c(0,1),legend.box.just="top",legend.background=element_rect(fill=alpha("white",0.4)))  )
		  
	print(ggplot(data.frame(pca$x), aes(x= pca$x[,1], y= pca$x[,2], colour=factor(infoTmp$Age), label=clonename))+
		  scale_colour_discrete(guide = FALSE) +geom_point() +geom_text(aes(label=clonename),hjust=0, vjust=0,size=2)+labs(title="Age")+ theme(legend.position=c(0,1),legend.justification=c(0,1),legend.box.just="top",legend.background=element_rect(fill=alpha("white",0.4)))  ) 
		  
	if(datSet=="combined"){print(ggplot(data.frame(pca$x), aes(x= pca$x[,1], y= pca$x[,2], colour=factor(infoTmp$datOrigin), label=clonename))+
		  scale_colour_discrete(guide =FALSE) +geom_point() +geom_text(aes(label=clonename),hjust=0, vjust=0,size=2)+labs(title="datOrigin")+ theme(legend.position=c(0,0),legend.justification=c(0,0),legend.box.just="top",legend.background=element_rect(fill=alpha("white",0.4))) )}
	
	#=========pca-1 vs pca-3======#
	print(ggplot(data.frame(pca$x), aes(x= pca$x[,1], y= pca$x[,3], colour=factor(infoTmp$Reprogramming.Source.Cell), label=clonename))+
		  scale_colour_discrete(guide =FALSE) +geom_point() +geom_text(aes(label=clonename),hjust=0, vjust=0,size=2)+labs(title="Reprogramming.Source.Cell")+ theme(legend.position=c(0,1),legend.justification=c(0,1),legend.box.just="top",legend.background=element_rect(fill=alpha("white",0.4))) )
		
	print(ggplot(data.frame(pca$x), aes(x= pca$x[,1], y= pca$x[,3], colour=factor(infoTmp$BATCH), label=clonename))+
		  scale_colour_discrete(guide = FALSE) +geom_point() +geom_text(aes(label=clonename),hjust=0, vjust=0,size=2)+labs(title="BATCH")+ theme(legend.position=c(0,1),legend.justification=c(0,1),legend.box.just="top",legend.background=element_rect(fill=alpha("white",0.4))) )
		
	print(ggplot(data.frame(pca$x), aes(x= pca$x[,1], y= pca$x[,3], colour=factor(infoTmp$Sample.Name), label=clonename))+
		  scale_colour_discrete(guide = FALSE) +geom_point() +geom_text(aes(label=clonename),hjust=0, vjust=0,size=2)+labs(title="Sample.Name")+ theme(legend.position=c(0,1),legend.justification=c(0,1),legend.box.just="top",legend.background=element_rect(fill=alpha("white",0.4))) )
		  
	print(ggplot(data.frame(pca$x), aes(x= pca$x[,1], y= pca$x[,3], colour=factor(infoTmp$Reprogramming.Batch), label=clonename))+
		  scale_colour_discrete(guide = FALSE) +geom_point() +geom_text(aes(label=clonename),hjust=0, vjust=0,size=2) +labs(title="Reprogramming.Batch")+ theme(legend.position=c(0,1),legend.justification=c(0,1),legend.box.just="top",legend.background=element_rect(fill=alpha("white",0.4))) )
		  
	if(datSet=="train"){print(ggplot(data.frame(pca$x), aes(x= pca$x[,1], y= pca$x[,3], colour=factor(infoTmp$Sendai.Virus.Lot), label=clonename))+
		  scale_colour_discrete(guide = FALSE) +geom_point() +geom_text(aes(label=clonename),hjust=0, vjust=0,size=2) +labs(title="Sendai.Virus.Lot") + theme(legend.position=c(0,1),legend.justification=c(0,1),legend.box.just="top",legend.background=element_rect(fill=alpha("white",0.4)))  )}
		
	print(ggplot(data.frame(pca$x), aes(x= pca$x[,1], y= pca$x[,3], colour=factor(infoTmp$SEX), label=clonename))+
		  scale_colour_discrete(guide = FALSE) +geom_point() +geom_text(aes(label=clonename),hjust=0, vjust=0,size=2) +labs(title="SEX") + theme(legend.position=c(0,1),legend.justification=c(0,1),legend.box.just="top",legend.background=element_rect(fill=alpha("white",0.4))) )
		  		  
	print(ggplot(data.frame(pca$x), aes(x= pca$x[,1], y= pca$x[,3], colour=factor(infoTmp$RACE), label=clonename))+
		  scale_colour_discrete(guide = FALSE) +geom_point() +geom_text(aes(label=clonename),hjust=0, vjust=0,size=2)+labs(title="RACE") + theme(legend.position=c(0,1),legend.justification=c(0,1),legend.box.just="top",legend.background=element_rect(fill=alpha("white",0.4))) )
		  
	print(ggplot(data.frame(pca$x), aes(x= pca$x[,1], y= pca$x[,3], colour=factor(infoTmp$State), label=clonename))+
		  scale_colour_discrete(guide = FALSE) +geom_point() +geom_text(aes(label=clonename),hjust=0, vjust=0,size=2) +labs(title="State") + theme(legend.position=c(0,1),legend.justification=c(0,1),legend.box.just="top",legend.background=element_rect(fill=alpha("white",0.4))) )
		  
	print(ggplot(data.frame(pca$x), aes(x= pca$x[,1], y= pca$x[,3], colour=factor(infoTmp$RNA.method), label=clonename))+
		  scale_colour_discrete(guide = FALSE) +geom_point() +geom_text(aes(label=clonename),hjust=0, vjust=0,size=2) +labs(title="RNA.method") + theme(legend.position=c(0,1),legend.justification=c(0,1),legend.box.just="top",legend.background=element_rect(fill=alpha("white",0.4))) )
		  
	print(ggplot(data.frame(pca$x), aes(x= pca$x[,1], y= pca$x[,3], colour=factor(infoTmp$Age), label=clonename))+
		  scale_colour_discrete(guide = FALSE) +geom_point() +geom_text(aes(label=clonename),hjust=0, vjust=0,size=2)+labs(title="Age") + theme(legend.position=c(0,1),legend.justification=c(0,1),legend.box.just="top",legend.background=element_rect(fill=alpha("white",0.4))) )
		  
	if(datSet=="combined"){print(ggplot(data.frame(pca$x), aes(x= pca$x[,1], y= pca$x[,3], colour=factor(infoTmp$datOrigin), label=clonename))+
		  scale_colour_discrete(guide =FALSE) +geom_point() +geom_text(aes(label=clonename),hjust=0, vjust=0,size=2)+labs(title="datOrigin")+ theme(legend.position=c(0,0),legend.justification=c(0,0),legend.box.just="top",legend.background=element_rect(fill=alpha("white",0.4))) )}
		  
	#=======pca-2 vs pca-3===========#
	print(ggplot(data.frame(pca$x), aes(x= pca$x[,2], y= pca$x[,3], colour=factor(infoTmp$Reprogramming.Source.Cell), label=clonename))+
		  scale_colour_discrete(guide =FALSE) +geom_point() +geom_text(aes(label=clonename),hjust=0, vjust=0,size=2)+labs(title="Reprogramming.Source.Cell")+ theme(legend.position=c(0,1),legend.justification=c(0,1),legend.box.just="top",legend.background=element_rect(fill=alpha("white",0.4))) )
		
	print(ggplot(data.frame(pca$x), aes(x= pca$x[,2], y= pca$x[,3], colour=factor(infoTmp$BATCH), label=clonename))+
		  scale_colour_discrete(guide = FALSE) +geom_point() +geom_text(aes(label=clonename),hjust=0, vjust=0,size=2)+labs(title="BATCH")+ theme(legend.position=c(0,1),legend.justification=c(0,1),legend.box.just="top",legend.background=element_rect(fill=alpha("white",0.4))) )
		
	print(ggplot(data.frame(pca$x), aes(x= pca$x[,2], y= pca$x[,3], colour=factor(infoTmp$Sample.Name), label=clonename))+
		  scale_colour_discrete(guide = FALSE) +geom_point() +geom_text(aes(label=clonename),hjust=0, vjust=0,size=2)+labs(title="Sample.Name")+ theme(legend.position=c(0,1),legend.justification=c(0,1),legend.box.just="top",legend.background=element_rect(fill=alpha("white",0.4))) )
		  
	print(ggplot(data.frame(pca$x), aes(x= pca$x[,2], y= pca$x[,3], colour=factor(infoTmp$Reprogramming.Batch), label=clonename))+
		  scale_colour_discrete(guide = FALSE) +geom_point() +geom_text(aes(label=clonename),hjust=0, vjust=0,size=2) +labs(title="Reprogramming.Batch")+ theme(legend.position=c(0,1),legend.justification=c(0,1),legend.box.just="top",legend.background=element_rect(fill=alpha("white",0.4))) )
		  
	if(datSet=="train"){print(ggplot(data.frame(pca$x), aes(x= pca$x[,2], y= pca$x[,3], colour=factor(infoTmp$Sendai.Virus.Lot), label=clonename))+
		  scale_colour_discrete(guide = FALSE) +geom_point() +geom_text(aes(label=clonename),hjust=0, vjust=0,size=2) +labs(title="Sendai.Virus.Lot")+ theme(legend.position=c(0,1),legend.justification=c(0,1),legend.box.just="top",legend.background=element_rect(fill=alpha("white",0.4))) )  }
		
	print(ggplot(data.frame(pca$x), aes(x= pca$x[,2], y= pca$x[,3], colour=factor(infoTmp$SEX), label=clonename))+
		  scale_colour_discrete(guide = FALSE) +geom_point() +geom_text(aes(label=clonename),hjust=0, vjust=0,size=2) +labs(title="SEX") + theme(legend.position=c(0,1),legend.justification=c(0,1),legend.box.just="top",legend.background=element_rect(fill=alpha("white",0.4))) )
		  
	print(ggplot(data.frame(pca$x), aes(x= pca$x[,2], y= pca$x[,3], colour=factor(infoTmp$RACE), label=clonename))+
		  scale_colour_discrete(guide = FALSE) +geom_point() +geom_text(aes(label=clonename),hjust=0, vjust=0,size=2)+labs(title="RACE") + theme(legend.position=c(0,1),legend.justification=c(0,1),legend.box.just="top",legend.background=element_rect(fill=alpha("white",0.4))) )
		  
	print(ggplot(data.frame(pca$x), aes(x= pca$x[,2], y= pca$x[,3], colour=factor(infoTmp$State), label=clonename))+
		  scale_colour_discrete(guide = FALSE) +geom_point() +geom_text(aes(label=clonename),hjust=0, vjust=0,size=2) +labs(title="State")+ theme(legend.position=c(0,1),legend.justification=c(0,1),legend.box.just="top",legend.background=element_rect(fill=alpha("white",0.4)))  )
		  
	print(ggplot(data.frame(pca$x), aes(x= pca$x[,2], y= pca$x[,3], colour=factor(infoTmp$RNA.method), label=clonename))+
		  scale_colour_discrete(guide = FALSE) +geom_point() +geom_text(aes(label=clonename),hjust=0, vjust=0,size=2) +labs(title="RNA.method") + theme(legend.position=c(0,1),legend.justification=c(0,1),legend.box.just="top",legend.background=element_rect(fill=alpha("white",0.4))) )
		  
	print(ggplot(data.frame(pca$x), aes(x= pca$x[,2], y= pca$x[,3], colour=factor(infoTmp$Age), label=clonename))+
		  scale_colour_discrete(guide = FALSE) +geom_point() +geom_text(aes(label=clonename),hjust=0, vjust=0,size=2)+labs(title="Age") + theme(legend.position=c(0,1),legend.justification=c(0,1),legend.box.just="top",legend.background=element_rect(fill=alpha("white",0.4))) )

	if(datSet=="combined"){print(ggplot(data.frame(pca$x), aes(x= pca$x[,2], y= pca$x[,3], colour=factor(infoTmp$datOrigin), label=clonename))+
		  scale_colour_discrete(guide =FALSE) +geom_point() +geom_text(aes(label=clonename),hjust=0, vjust=0,size=2)+labs(title="datOrigin")+ theme(legend.position=c(0,0),legend.justification=c(0,0),legend.box.just="top",legend.background=element_rect(fill=alpha("white",0.4))) )}

	dev.off()
}



#***************************************************************************************************************************************************
#--------Residuals by removing BATCH+SEX+Reprogramming.Source.Cell+SendaiVirusLot+Reprogramming.Batch+Reprogramming.Source.Cell+Dataset-------------
#Remove Batch Effect from RNA-seq data and normalized by voom
#No hidden factors are used. This correction only KNOWN COVARIATES
#BATCH+SEX+Reprogramming.Source.Cell+SendaiVirusLot+Reprogramming.Batch+Reprogramming.Source.Cell+Dataset
#****************************************************************
#---------------------
#  data version-3
#---------------------
info = read.table("combined_withnewdata_317samples_covariates_IRvsIS_IVAN.txt", header=TRUE, sep="\t")
info = info[match(colnames(sample), as.character(info$Sample.Name)),]
tmpPatIDs<-as.character(info$PatientID)
for(i in 1:length(tmpPatIDs)){if(length(unlist(strsplit(split="",tmpPatIDs[i])))==2){tmpPatIDs[i]<-paste(sep="","0",tmpPatIDs[i])}}
info$PatientID<-factor(tmpPatIDs)
info$RACE[info$RACE=="Black"]<-"African.American" # consolidate 'Black' & 'African.American' race labels
info$RACE[is.element(el=info$PatientID,set=c("054","912","212","285","729"))]<-"East.Asian" # correct IDs of patient recorded as simply Asian"
info$RACE<-factor(info$RACE)

#========Marc's design matrix=======#
# the idea is to normalise & adjust the IR sample separately from the IS samples while adjusting for patient ID to remove the intra-patient variability without automatically also loosing the IR/IS variability
# NB all covariates (sex, age, bmi, race, reprogramming.source.cell, rna.method) are co-linear with patient ID, so just correct for that one
idxIR<-which(info$State=="Insulin resistant")
irSampIDs<-info$Sample.Name[idxIR]
TRUE.GENE_EXPRESSION_DGELIST_MAT.NORM.IR<-TRUE.GENE_EXPRESSION_DGELIST_MAT.NORM[,is.element(el=colnames(TRUE.GENE_EXPRESSION_DGELIST_MAT.NORM),set=irSampIDs)]
TRUE.GENE_EXPRESSION_DGELIST_MAT.NORM.IS<-TRUE.GENE_EXPRESSION_DGELIST_MAT.NORM[,!is.element(el=colnames(TRUE.GENE_EXPRESSION_DGELIST_MAT.NORM),set=irSampIDs)]
info.IR<-info[idxIR,]
info.IR$RACE<-factor(info.IR$RACE)
info.IR$PatientID<-factor(info.IR$PatientID)
info.IS<-info[-idxIR,]
info.IS$PatientID<-factor(info.IS$PatientID)
info.IS$RACE<-factor(info.IS$RACE)

design.Common<-model.matrix(~Reprogramming.Source.Cell+SEX+RACE+Age+BMI,info) # even though each of these is colinear with PatientID, we want to remove these effects first as we want to estimate a patient-only effect, without contributions for sex, age, race or technical effects etc; this is important, otherwise e.g. the intercept for male patients for Y-chr genes will be unecessarily high
design.Common.withState<-model.matrix(~Reprogramming.Source.Cell+SEX+RACE+Age+BMI+State,info)
adjusted_covariates="VP--BATCH_RNA.method_lmFit--Reprogramming.Source.Cell_SEX_RACE_Age_BMI" # later in the code we run variancePartition for BATCH and RNA.method first and then lmFit for the above design matrix
adjusted_covariates.withState="VP--BATCH_RNA.method_lmFit--Reprogramming.Source.Cell_SEX_RACE_Age_BMI_State" # later in the code we run variancePartition for BATCH and RNA.method first and then lmFit for the above design matrix

#-- create info matrix for avg count matrix (can also be used for avg count post normalisation)
cols2Remove<-which(is.element(el=colnames(info),set=c("Sample.Name.old","Sample.Name.Original","Sample.Name","Reprogramming.Batch","Sendai.Virus.Lot","RIN","Passage","totalReads","sendai","sendaiCPM")))
idx4Info<-integer(0)
for(patID in uniqPatIDs){
	idxTmp<-which(as.character(info$PatientID)==patID)
	if(length(idxTmp)>1){
		tmpCheck<-character(length(idxTmp))
		for(j in setdiff(1:ncol(info),cols2Remove)){
			tmpCheck<-paste(sep="::",tmpCheck,info[idxTmp,j])
		}
		if(sum(tmpCheck==tmpCheck[1])!=length(tmpCheck)){
			stop(paste(sep="","ERROR: not all covariates are identical for different clones from patient ",patID,"!"))
		}
	}
	idx4Info<-c(idx4Info,idxTmp[1])
}
info.PATIENT<-info[idx4Info,-cols2Remove]
info.PATIENT$PatientID2<-uniqPatIDs
design.Common.PATIENT<-model.matrix(~SEX+Reprogramming.Source.Cell+RACE+Age+BMI,info.PATIENT)

# repeat all this for the holdout samples
#info.Holdout<-read.table("iPSC_metadata_v2.csv",sep=",",stringsAsFactors=F,header=T)
info.Holdout<-read.table("newSamples_iPSC_personalAndTechnicalCovariates.tab",sep="\t",stringsAsFactors=F,header=T)
## remove samples from patients 108 and 946 (each has their own, unique batch) and the clones from patient 148 that are in batch B786
#remove.Holdout.IDs<-info.Holdout$Sample.Name[info.Holdout$PatientID=="148" & info.Holdout$BATCH=="B786"]
#sampleHoldout<-sampleHoldout[,!is.element(el=colnames(sampleHoldout),set=remove.Holdout.IDs)]

info.Holdout<-info.Holdout[match(colnames(sampleHoldout),info.Holdout$Sample.Name),]
tmpPatIDs<-as.character(info.Holdout$PatientID)
for(i in 1:length(tmpPatIDs)){if(length(unlist(strsplit(split="",tmpPatIDs[i])))==2){tmpPatIDs[i]<-paste(sep="","0",tmpPatIDs[i])}}
info.Holdout$PatientID<-factor(tmpPatIDs)
info.Holdout$RACE<-as.character(info.Holdout$RACE)
info.Holdout$RACE[info.Holdout$RACE=="Black" | info.Holdout$RACE=="Black or African American"]<-"African.American" # consolidate redundant labels
info.Holdout$RACE[info.Holdout$RACE=="East Asian"]<-"East.Asian" # use same labels as the info for the 317 samples
info.Holdout$RACE[info.Holdout$RACE=="South Asian"]<-"South.Asian" # use same labels as the info for the 317 samples
info.Holdout$RACE[info.Holdout$RACE=="White Hispanic"]<-"White.Hispanic" # use same labels as the info for the 317 samples
info.Holdout$RACE[is.element(el=info.Holdout$PatientID,set=c("054","912","212","285","729"))]<-"East.Asian" # correct IDs of patient recorded as simply Asian"
info.Holdout$RACE<-factor(info.Holdout$RACE)

design.Common.Holdout<-model.matrix(~SEX+RACE+Age+BMI,info.Holdout) # no need to adjust for Reprogramming.source.cell or RNA.method -- just one cell type for these 85 samples
#design.Common.Holdout<-design.Common.Holdout[,colnames(design.Common.Holdout)!="RACEEast Asian"] # east asian samples all from one patient whose samples were processed in just one batch

cols2Remove.Holdout<-which(is.element(el=colnames(info.Holdout),set=c("Sample.Name","REASON","Reprogramming.Batch","Sendai.Virus.Lot","RIN","Passage","totalReads","sendai","sendaiCPM","REAL.NAME")))
idx4Info.Holdout<-integer(0)
for(patID in uniqPatIDHOs){
	idxTmp.HO<-which(as.character(info.Holdout$PatientID)==patID)
	if(length(idxTmp.HO)>1){
		tmpCheck<-character(length(idxTmp.HO))
		for(j in setdiff(1:ncol(info.Holdout),cols2Remove.Holdout)){
			tmpCheck<-paste(sep="::",tmpCheck,info.Holdout[idxTmp.HO,j])
		}
		if(sum(tmpCheck==tmpCheck[1])!=length(tmpCheck)){
			stop(paste(sep="","ERROR: not all covariates are identical for different clones from patient ",patID,"!"))
		}
	}
	idx4Info.Holdout<-c(idx4Info.Holdout,idxTmp.HO[1])
}
info.Holdout.PATIENT<-info.Holdout[idx4Info.Holdout,-cols2Remove.Holdout]
design.Common.Holdout.PATIENT<-model.matrix(~SEX+RACE+Age+BMI,info.Holdout.PATIENT) # cannot adjust for BATCH since different clones for same patients may have different BATCH values -- it really is just patient 148 that has 2 different batches for a total of 5 clones; NB patients 108 and 946 have their own unique batches when the holdout samples are not combined with the other 317/310 samples (? remove the 2 clones for patient 148 in batch B786)

# create common info matrix
commonCovs<-intersect(colnames(info),colnames(info.Holdout))
info.All<-rbind(info[,match(commonCovs,colnames(info))],info.Holdout[,match(commonCovs,colnames(info.Holdout))])
datOrigin<-c(rep("original317",nrow(info)),rep("holdout85",nrow(info.Holdout)))
info.All<-data.frame(info.All,datOrigin=datOrigin)
info.All<-info.All[match(colnames(sampleAll),info.All$Sample.Name),]
design.Common.All<-model.matrix(~Reprogramming.Source.Cell+SEX+RACE+Age+BMI,info.All)

# create avg count info matrix for all (combined) data
cols2Remove.All<-which(is.element(el=colnames(info.All),set=c("Sample.Name","Sample.Name.Original","Sample.Name","REASON","Reprogramming.Batch","Sendai.Virus.Lot","RIN","Passage","totalReads","sendai","sendaiCPM","REAL.NAME")))
idx4Info.All<-integer(0)
tmpPatIDs<-as.character(info.All$PatientID)
for(i in 1:length(tmpPatIDs)){if(length(unlist(strsplit(split="",tmpPatIDs[i])))==2){tmpPatIDs[i]<-paste(sep="","0",tmpPatIDs[i])}}
info.All$PatientID<-factor(tmpPatIDs)
for(patID in uniqPatIDAlls){
	idxTmp.All<-which(as.character(info.All$PatientID)==patID)
	if(length(idxTmp.All)>1){
		tmpCheck<-character(length(idxTmp.All))
		for(j in setdiff(1:ncol(info.All),cols2Remove.All)){
			tmpCheck<-paste(sep="::",tmpCheck,info.All[idxTmp.All,j])
		}
		if(sum(tmpCheck==tmpCheck[1])!=length(tmpCheck)){
			stop(paste(sep="","ERROR: not all covariates are identical for different clones from patient ",patID,"!"))
		}
	}
	idx4Info.All<-c(idx4Info.All,idxTmp.All[1])
}
info.All.PATIENT<-info.All[idx4Info.All,-cols2Remove.All]
design.Common.All.PATIENT<-model.matrix(~Reprogramming.Source.Cell+SEX+RACE+Age+BMI,info.All.PATIENT)
write.table(info.All.PATIENT,sep="\t",row.names=F,col.names=T,quote=F,file="combined_317samples-Train_plus_85samples-Holdout_covariates_IRvsIS_IVAN_PATIENT-LEVEL.txt")

#===================================#

#Apply voom
PRED.GENE_EXPRESSION_MAT<- voom(TRUE.GENE_EXPRESSION_DGELIST_MAT.NORM,design=design.Common,plot=TRUE) 

#Get the numeric matrix of normalized expression values on the log2 scale
VOOM_NORMALIZED_LOG_EXPRESSION_MAT = PRED.GENE_EXPRESSION_MAT$E

#Get the numeric matrix of inverse variance weights
VOOM_WEIGHTS_MAT = PRED.GENE_EXPRESSION_MAT$weights
colnames(VOOM_WEIGHTS_MAT)<-colnames(PRED.GENE_EXPRESSION_MAT$E)
VOOM_WEIGHTS_MAT.IR<-VOOM_WEIGHTS_MAT[,match(info.IR$Sample.Name,colnames(VOOM_WEIGHTS_MAT))]
VOOM_WEIGHTS_MAT.IS<-VOOM_WEIGHTS_MAT[,match(info.IS$Sample.Name,colnames(VOOM_WEIGHTS_MAT))]

# compute residuals -- variancePartition for BATCH, RNA.method first, then lmFit for SEX, RACE, Age, BMI
fit.Common.VarPart <- fitVarPartModel( exprObj=VOOM_NORMALIZED_LOG_EXPRESSION_MAT, ~ (1|BATCH) + (1|RNA.method), data=info , useWeights=TRUE, weightsMatrix=VOOM_WEIGHTS_MAT)
Res.KNOWN.AllKNOWN.VarPart <- residuals( fit.Common.VarPart )
colnames(Res.KNOWN.AllKNOWN.VarPart)<-colnames(VOOM_NORMALIZED_LOG_EXPRESSION_MAT)
fit.Common=lmFit(Res.KNOWN.AllKNOWN.VarPart, design.Common)
Res.KNOWN.AllKNOWN = residuals(fit.Common, Res.KNOWN.AllKNOWN.VarPart) # no intercept added back in!
save(Res.KNOWN.AllKNOWN,file=paste("output/Res.KNOWN",adjusted_covariates,dataversion,low_expression_cutoff,"RData",sep="."))
write.table(Res.KNOWN.AllKNOWN,quote=FALSE,sep="\t",col.names=TRUE,row.names=TRUE,file=paste("output/Res.KNOWN",adjusted_covariates,dataversion,low_expression_cutoff,"txt",sep="."))
Res.KNOWN.AllKNOWN.scale=t(scale(t(Res.KNOWN.AllKNOWN),scale = FALSE, center = TRUE))
#save(Res.KNOWN.AllKNOWN.scale,file=paste("output/Res.KNOWN",adjusted_covariates,dataversion,low_expression_cutoff,"centered","RData",sep="."))
#write.table(Res.KNOWN.AllKNOWN.scale,quote=FALSE,sep="\t",col.names=TRUE,row.names=TRUE,file=paste("output/Res.KNOWN",adjusted_covariates,dataversion,low_expression_cutoff,"centered","txt",sep="."))

fit.Common.withState=lmFit(Res.KNOWN.AllKNOWN.VarPart, design.Common.withState)
Res.KNOWN.AllKNOWN.withState = residuals(fit.Common.withState, Res.KNOWN.AllKNOWN.VarPart) # no intercept added back in!
save(Res.KNOWN.AllKNOWN.withState,file=paste("output/Res.KNOWN",adjusted_covariates.withState,dataversion,low_expression_cutoff,"RData",sep="."))
write.table(Res.KNOWN.AllKNOWN.withState,quote=FALSE,sep="\t",col.names=TRUE,row.names=TRUE,file=paste("output/Res.KNOWN",adjusted_covariates.withState,dataversion,low_expression_cutoff,"txt",sep="."))


# get average residuals
patIDs<-gsub(pattern="_[0-9]*",replacement="",colnames(Res.KNOWN.AllKNOWN))
uniqPatIDs<-unique(patIDs)
Res.KNOWN.AllKNOWN.AvgRes.PATIENT<-matrix(nrow=nrow(Res.KNOWN.AllKNOWN),ncol=length(uniqPatIDs))
rownames(Res.KNOWN.AllKNOWN.AvgRes.PATIENT)<-rownames(Res.KNOWN.AllKNOWN)
colnames(Res.KNOWN.AllKNOWN.AvgRes.PATIENT)<-uniqPatIDs
for(patID in uniqPatIDs){
	if(sum(patIDs==patID)>1){
		Res.KNOWN.AllKNOWN.AvgRes.PATIENT[,colnames(Res.KNOWN.AllKNOWN.AvgRes.PATIENT)==patID]<-rowMeans(Res.KNOWN.AllKNOWN[,patIDs==patID])
	}else{
		Res.KNOWN.AllKNOWN.AvgRes.PATIENT[,colnames(Res.KNOWN.AllKNOWN.AvgRes.PATIENT)==patID]<-Res.KNOWN.AllKNOWN[,patIDs==patID]
	}
}
save(Res.KNOWN.AllKNOWN.AvgRes.PATIENT,file=paste("output/Res.KNOWN",adjusted_covariates,dataversion,low_expression_cutoff,"AvgResPerPatient","RData",sep="."))
write.table(Res.KNOWN.AllKNOWN.AvgRes.PATIENT,quote=FALSE,sep="\t",col.names=TRUE,row.names=TRUE,file=paste("output/Res.KNOWN",adjusted_covariates,dataversion,low_expression_cutoff,"AvgResPerPatient","txt",sep="."))

# repeat the above for the avg. counts matrix
PRED.GENE_EXPRESSION_MAT.PATIENT<- voom(TRUE.GENE_EXPRESSION_DGELIST_MAT.PATIENT.NORM,design=design.Common.PATIENT,plot=TRUE) 
VOOM_NORMALIZED_LOG_EXPRESSION_MAT.PATIENT = PRED.GENE_EXPRESSION_MAT.PATIENT$E
VOOM_WEIGHTS_MAT.PATIENT = PRED.GENE_EXPRESSION_MAT.PATIENT$weights
colnames(VOOM_WEIGHTS_MAT.PATIENT)<-colnames(PRED.GENE_EXPRESSION_MAT.PATIENT$E)
fit.Common.PATIENT.VarPart <- fitVarPartModel( exprObj=VOOM_NORMALIZED_LOG_EXPRESSION_MAT.PATIENT, ~ (1|BATCH) + (1|RNA.method), data=info.PATIENT , useWeights=TRUE, weightsMatrix=VOOM_WEIGHTS_MAT.PATIENT)
Res.KNOWN.AllKNOWN.PATIENT.VarPart <- residuals( fit.Common.PATIENT.VarPart )
colnames(Res.KNOWN.AllKNOWN.PATIENT.VarPart)<-colnames(VOOM_NORMALIZED_LOG_EXPRESSION_MAT.PATIENT)
fit.Common.PATIENT=lmFit(Res.KNOWN.AllKNOWN.PATIENT.VarPart, design.Common.PATIENT)
Res.KNOWN.AllKNOWN.PATIENT = residuals(fit.Common.PATIENT, Res.KNOWN.AllKNOWN.PATIENT.VarPart) # no intercept added back in!
save(Res.KNOWN.AllKNOWN.PATIENT,file=paste("output/Res.KNOWN",adjusted_covariates,dataversion,low_expression_cutoff,"AvgCountPerPatient","RData",sep="."))
write.table(Res.KNOWN.AllKNOWN.PATIENT,quote=FALSE,sep="\t",col.names=TRUE,row.names=TRUE,file=paste("output/Res.KNOWN",adjusted_covariates,dataversion,low_expression_cutoff,"AvgCountPerPatient","txt",sep="."))
Res.KNOWN.AllKNOWN.PATIENT.scale=t(scale(t(Res.KNOWN.AllKNOWN.PATIENT),scale = FALSE, center = TRUE))
#save(Res.KNOWN.AllKNOWN.PATIENT.scale,file=paste("output/Res.KNOWN",adjusted_covariates,dataversion,low_expression_cutoff,"centered","RData",sep="."))
#write.table(Res.KNOWN.AllKNOWN.PATIENT.scale,quote=FALSE,sep="\t",col.names=TRUE,row.names=TRUE,file=paste("output/Res.KNOWN",adjusted_covariates,dataversion,low_expression_cutoff,"centered","txt",sep="."))

# repeat the above for the avg. counts post normalisation matrix
# new: create library size & voom normalised matrix of training data by combining the library size & voomed data from the individual training samples
# i.e. compute avg values par patient from vobj$E -- the non-adjusted but library-size and voom normalised data (NB the E matrices are actually the same for vobj and PRED.GENE_EXPRESSION_MAT)
# NB we also need to compute and average the weight matrix -- for this we will use PRED.GENE_EXPRESSION_MAT since its weights are computed using the coefficients for the covariates we will adjust for
if(sum(rownames(vobj$E)!=rownames(PRED.GENE_EXPRESSION_MAT$E))>0){stop("vobj$E and PRED.GENE_EXPRESSION_MAT$E have to have samples in the same order!")}
if(sum(colnames(vobj$E)!=colnames(PRED.GENE_EXPRESSION_MAT$E))>0){stop("vobj$E and PRED.GENE_EXPRESSION_MAT$E have to have genes in the same order!")}
patIDs.PostLibSizeAdj<-gsub(pattern="_[0-9]*",replacement="",colnames(vobj$E))
uniqPatIDs.PostLibSizeAdj<-unique(patIDs.PostLibSizeAdj)
VOOM_RAW_LOG_EXPRESSION_MAT.PATIENT.PostLibSizeAdj<-matrix(nrow=nrow(vobj$E),ncol=length(uniqPatIDs.PostLibSizeAdj))
VOOMWEIGHTS_LOG_EXPRESSION_MAT.PATIENT.PostLibSizeAdj<-matrix(nrow=nrow(vobj$weights),ncol=length(uniqPatIDs.PostLibSizeAdj))
rownames(VOOM_RAW_LOG_EXPRESSION_MAT.PATIENT.PostLibSizeAdj)<-rownames(vobj$E)
colnames(VOOM_RAW_LOG_EXPRESSION_MAT.PATIENT.PostLibSizeAdj)<-uniqPatIDs.PostLibSizeAdj
rownames(VOOMWEIGHTS_LOG_EXPRESSION_MAT.PATIENT.PostLibSizeAdj)<-rownames(VOOM_RAW_LOG_EXPRESSION_MAT.PATIENT.PostLibSizeAdj)
colnames(VOOMWEIGHTS_LOG_EXPRESSION_MAT.PATIENT.PostLibSizeAdj)<-uniqPatIDs.PostLibSizeAdj
for(patID in uniqPatIDs.PostLibSizeAdj){
	if(sum(patIDs.PostLibSizeAdj==patID)>1){
		VOOM_RAW_LOG_EXPRESSION_MAT.PATIENT.PostLibSizeAdj[,colnames(VOOM_RAW_LOG_EXPRESSION_MAT.PATIENT.PostLibSizeAdj)==patID]<-rowMeans(vobj$E[,patIDs.PostLibSizeAdj==patID])
		VOOMWEIGHTS_LOG_EXPRESSION_MAT.PATIENT.PostLibSizeAdj[,colnames(VOOMWEIGHTS_LOG_EXPRESSION_MAT.PATIENT.PostLibSizeAdj)==patID]<-rowMeans(PRED.GENE_EXPRESSION_MAT$weights[,patIDs.PostLibSizeAdj==patID])
	}else{
		VOOM_RAW_LOG_EXPRESSION_MAT.PATIENT.PostLibSizeAdj[,colnames(VOOM_RAW_LOG_EXPRESSION_MAT.PATIENT.PostLibSizeAdj)==patID]<-VOOM_RAW_LOG_EXPRESSION_MAT[,patIDs.PostLibSizeAdj==patID]
		VOOMWEIGHTS_LOG_EXPRESSION_MAT.PATIENT.PostLibSizeAdj[,colnames(VOOMWEIGHTS_LOG_EXPRESSION_MAT.PATIENT.PostLibSizeAdj)==patID]<-PRED.GENE_EXPRESSION_MAT$weights[,patIDs.PostLibSizeAdj==patID]
	}
}
# create vobj.PATIENT.PostLibSizeAdj
targets.Tmp<-data.frame(group=rep(1,ncol(VOOM_RAW_LOG_EXPRESSION_MAT.PATIENT.PostLibSizeAdj)),lib.size=rep(mean(vobj$targets$lib.size),ncol(VOOM_RAW_LOG_EXPRESSION_MAT.PATIENT.PostLibSizeAdj)),norm.factors=rep(1,ncol(VOOM_RAW_LOG_EXPRESSION_MAT.PATIENT.PostLibSizeAdj))) # since this is after normalisation we assume all library sizes and normalisation factors to be equal
rownames(targets.Tmp)<-colnames(VOOM_RAW_LOG_EXPRESSION_MAT.PATIENT.PostLibSizeAdj)
vobj.PATIENT.PostLibSizeAdj<-list(genes=list(rownames(VOOM_RAW_LOG_EXPRESSION_MAT.PATIENT.PostLibSizeAdj)),targets=targets.Tmp,E=VOOM_RAW_LOG_EXPRESSION_MAT.PATIENT.PostLibSizeAdj,weights=VOOMWEIGHTS_LOG_EXPRESSION_MAT.PATIENT.PostLibSizeAdj,design=NULL)
vobj.PATIENT.PostLibSizeAdj<-as(vobj.PATIENT.PostLibSizeAdj,Class="EList")

VOOM_NORMALIZED_LOG_EXPRESSION_MAT.PATIENT.PostLibSizeAdj<-vobj.PATIENT.PostLibSizeAdj$E
VOOM_WEIGHTS_MAT.PATIENT.PostLibSizeAdj = vobj.PATIENT.PostLibSizeAdj$weights
colnames(VOOM_WEIGHTS_MAT.PATIENT.PostLibSizeAdj)<-colnames(vobj.PATIENT.PostLibSizeAdj$E)
fit.Common.PATIENT.PostLibSizeAdj.VarPart <- fitVarPartModel( exprObj=VOOM_NORMALIZED_LOG_EXPRESSION_MAT.PATIENT.PostLibSizeAdj, ~ (1|BATCH) + (1|RNA.method), data=info.PATIENT , useWeights=TRUE, weightsMatrix=VOOM_WEIGHTS_MAT.PATIENT.PostLibSizeAdj)
Res.KNOWN.AllKNOWN.PATIENT.PostLibSizeAdj.VarPart <- residuals( fit.Common.PATIENT.PostLibSizeAdj.VarPart )
colnames(Res.KNOWN.AllKNOWN.PATIENT.PostLibSizeAdj.VarPart)<-colnames(VOOM_NORMALIZED_LOG_EXPRESSION_MAT.PATIENT.PostLibSizeAdj)
fit.Common.PATIENT.PostLibSizeAdj=lmFit(Res.KNOWN.AllKNOWN.PATIENT.PostLibSizeAdj.VarPart, design.Common.PATIENT, weights=VOOM_WEIGHTS_MAT.PATIENT.PostLibSizeAdj)
Res.KNOWN.AllKNOWN.PATIENT.PostLibSizeAdj = residuals(fit.Common.PATIENT.PostLibSizeAdj, Res.KNOWN.AllKNOWN.PATIENT.PostLibSizeAdj.VarPart) # no intercept added back in!
save(Res.KNOWN.AllKNOWN.PATIENT.PostLibSizeAdj,file=paste("output/Res.KNOWN",adjusted_covariates,dataversion,low_expression_cutoff,"AvgCountPerPatient.PostLibSizeAdj","RData",sep="."))
write.table(Res.KNOWN.AllKNOWN.PATIENT.PostLibSizeAdj,quote=FALSE,sep="\t",col.names=TRUE,row.names=TRUE,file=paste("output/Res.KNOWN",adjusted_covariates,dataversion,low_expression_cutoff,"AvgCountPerPatient.PostLibSizeAdj","txt",sep="."))
Res.KNOWN.AllKNOWN.PATIENT.PostLibSizeAdj.scale=t(scale(t(Res.KNOWN.AllKNOWN.PATIENT.PostLibSizeAdj),scale = FALSE, center = TRUE))

# repeat the above for the holdout samples
PRED.GENE_EXPRESSION_MAT.Holdout<- voom(TRUE.GENE_EXPRESSION_DGELIST_MAT_Holdout.NORM,design=design.Common.Holdout,plot=TRUE) 
VOOM_NORMALIZED_LOG_EXPRESSION_MAT.Holdout = PRED.GENE_EXPRESSION_MAT.Holdout$E
VOOM_WEIGHTS_MAT.Holdout = PRED.GENE_EXPRESSION_MAT.Holdout$weights
colnames(VOOM_WEIGHTS_MAT.Holdout)<-colnames(PRED.GENE_EXPRESSION_MAT.Holdout$E)
adjusted_covariatesHoldout="VP--BATCH_lmFit--SEX_RACE_Age_BMI"
fit.Common.Holdout.VarPart <- fitVarPartModel( exprObj=VOOM_NORMALIZED_LOG_EXPRESSION_MAT.Holdout, ~ (1|BATCH), data=info.Holdout , useWeights=TRUE, weightsMatrix=VOOM_WEIGHTS_MAT.Holdout) # no need to adjust for Reprogramming.Source.Cell or RNA.method -- only 1 level for all samples
Res.KNOWN.AllKNOWN.Holdout.VarPart <- residuals( fit.Common.Holdout.VarPart )
colnames(Res.KNOWN.AllKNOWN.Holdout.VarPart)<-colnames(VOOM_NORMALIZED_LOG_EXPRESSION_MAT.Holdout) 
fit.Common.Holdout=lmFit(Res.KNOWN.AllKNOWN.Holdout.VarPart, design.Common.Holdout, weights=VOOM_WEIGHTS_MAT.Holdout)
Res.KNOWN.AllKNOWN.Holdout = residuals(fit.Common.Holdout, Res.KNOWN.AllKNOWN.Holdout.VarPart) # no intercept added back in!
save(Res.KNOWN.AllKNOWN.Holdout,file=paste("output/Res.KNOWN",adjusted_covariatesHoldout,dataversion,low_expression_cutoff,"Holdout","RData",sep="."))
write.table(Res.KNOWN.AllKNOWN.Holdout,quote=FALSE,sep="\t",col.names=TRUE,row.names=TRUE,file=paste("output/Res.KNOWN",adjusted_covariatesHoldout,dataversion,low_expression_cutoff,"Holdout","txt",sep="."))

# get average residuals
patIDs.Holdout<-gsub(pattern="_[0-9]*",replacement="",colnames(Res.KNOWN.AllKNOWN.Holdout))
uniqPatIDs.Holdout<-unique(patIDs.Holdout)
Res.KNOWN.AllKNOWN.Holdout.AvgRes.PATIENT<-matrix(nrow=nrow(Res.KNOWN.AllKNOWN.Holdout),ncol=length(uniqPatIDs.Holdout))
rownames(Res.KNOWN.AllKNOWN.Holdout.AvgRes.PATIENT)<-rownames(Res.KNOWN.AllKNOWN.Holdout)
colnames(Res.KNOWN.AllKNOWN.Holdout.AvgRes.PATIENT)<-uniqPatIDs.Holdout
for(patID in uniqPatIDs.Holdout){
	if(sum(patIDs.Holdout==patID)>1){
		Res.KNOWN.AllKNOWN.Holdout.AvgRes.PATIENT[,colnames(Res.KNOWN.AllKNOWN.Holdout.AvgRes.PATIENT)==patID]<-rowMeans(Res.KNOWN.AllKNOWN.Holdout[,patIDs.Holdout==patID])
	}else{
		Res.KNOWN.AllKNOWN.Holdout.AvgRes.PATIENT[,colnames(Res.KNOWN.AllKNOWN.Holdout.AvgRes.PATIENT)==patID]<-Res.KNOWN.AllKNOWN.Holdout[,patIDs.Holdout==patID]
	}
}
save(Res.KNOWN.AllKNOWN.Holdout.AvgRes.PATIENT,file=paste("output/Res.KNOWN",adjusted_covariatesHoldout,dataversion,low_expression_cutoff,"Holdout","AvgResPerPatient","RData",sep="."))
write.table(Res.KNOWN.AllKNOWN.Holdout.AvgRes.PATIENT,quote=FALSE,sep="\t",col.names=TRUE,row.names=TRUE,file=paste("output/Res.KNOWN",adjusted_covariatesHoldout,dataversion,low_expression_cutoff,"Holdout","AvgResPerPatient","txt",sep="."))

# repeat the above for the holdout avg. counts matrix
PRED.GENE_EXPRESSION_MAT.Holdout.PATIENT<- voom(TRUE.GENE_EXPRESSION_DGELIST_MAT_Holdout.PATIENT.NORM,design=design.Common.Holdout.PATIENT,plot=TRUE) 
VOOM_NORMALIZED_LOG_EXPRESSION_MAT.Holdout.PATIENT = PRED.GENE_EXPRESSION_MAT.Holdout.PATIENT$E
VOOM_WEIGHTS_MAT.Holdout.PATIENT = PRED.GENE_EXPRESSION_MAT.Holdout.PATIENT$weights
colnames(VOOM_WEIGHTS_MAT.Holdout.PATIENT)<-colnames(PRED.GENE_EXPRESSION_MAT.Holdout.PATIENT$E)
fit.Common.Holdout.PATIENT.VarPart <- fitVarPartModel( exprObj=VOOM_NORMALIZED_LOG_EXPRESSION_MAT.Holdout.PATIENT, ~ (1|BATCH), data=info.Holdout.PATIENT , useWeights=TRUE, weightsMatrix=VOOM_WEIGHTS_MAT.Holdout.PATIENT)
Res.KNOWN.AllKNOWN.Holdout.PATIENT.VarPart <- residuals( fit.Common.Holdout.PATIENT.VarPart )
colnames(Res.KNOWN.AllKNOWN.Holdout.PATIENT.VarPart)<-colnames(VOOM_NORMALIZED_LOG_EXPRESSION_MAT.Holdout.PATIENT)
fit.Common.Holdout.PATIENT=lmFit(Res.KNOWN.AllKNOWN.Holdout.PATIENT.VarPart, design.Common.Holdout.PATIENT, weights=VOOM_WEIGHTS_MAT.Holdout.PATIENT)
Res.KNOWN.AllKNOWN.Holdout.PATIENT = residuals(fit.Common.Holdout.PATIENT, Res.KNOWN.AllKNOWN.Holdout.PATIENT.VarPart) # no intercept added back in!
save(Res.KNOWN.AllKNOWN.Holdout.PATIENT,file=paste("output/Res.KNOWN",adjusted_covariatesHoldout,dataversion,low_expression_cutoff,"Holdout","AvgCountPerPatient","RData",sep="."))
write.table(Res.KNOWN.AllKNOWN.Holdout.PATIENT,quote=FALSE,sep="\t",col.names=TRUE,row.names=TRUE,file=paste("output/Res.KNOWN",adjusted_covariatesHoldout,dataversion,low_expression_cutoff,"Holdout","AvgCountPerPatient","txt",sep="."))

# repeat the above for the holdout avg. counts post normalisation matrix
# new: create library size & voom normalised matrix of holdout data by combining the library size & voomed data from the individual holdout samples
# i.e. compute avg values par patient from vobjHO$E -- the non-adjusted but library-size and voom normalised data (NB the E matrices are actually the same for vobjHO and PRED.GENE_EXPRESSION_MAT.Holdout)
# NB we also need to compute and average the weight matrix -- for this we will use PRED.GENE_EXPRESSION_MAT.Holdout since its weights are computed using the coefficients for the covariates we will adjust for
if(sum(rownames(vobjHO$E)!=rownames(PRED.GENE_EXPRESSION_MAT.Holdout$E))>0){stop("vobjHO$E and PRED.GENE_EXPRESSION_MAT.Holdout$E have to have samples in the same order!")}
if(sum(colnames(vobjHO$E)!=colnames(PRED.GENE_EXPRESSION_MAT.Holdout$E))>0){stop("vobjHO$E and PRED.GENE_EXPRESSION_MAT.Holdout$E have to have genes in the same order!")}
patIDs.Holdout.PostLibSizeAdj<-gsub(pattern="_[0-9]*",replacement="",colnames(vobjHO$E))
uniqPatIDs.Holdout.PostLibSizeAdj<-unique(patIDs.Holdout.PostLibSizeAdj)
VOOM_RAW_LOG_EXPRESSION_MAT_Holdout.PATIENT.PostLibSizeAdj<-matrix(nrow=nrow(vobjHO$E),ncol=length(uniqPatIDs.Holdout.PostLibSizeAdj))
VOOMWEIGHTS_LOG_EXPRESSION_MAT_Holdout.PATIENT.PostLibSizeAdj<-matrix(nrow=nrow(vobjHO$weights),ncol=length(uniqPatIDs.Holdout.PostLibSizeAdj))
rownames(VOOM_RAW_LOG_EXPRESSION_MAT_Holdout.PATIENT.PostLibSizeAdj)<-rownames(vobjHO$E)
colnames(VOOM_RAW_LOG_EXPRESSION_MAT_Holdout.PATIENT.PostLibSizeAdj)<-uniqPatIDs.Holdout.PostLibSizeAdj
rownames(VOOMWEIGHTS_LOG_EXPRESSION_MAT_Holdout.PATIENT.PostLibSizeAdj)<-rownames(VOOM_RAW_LOG_EXPRESSION_MAT_Holdout.PATIENT.PostLibSizeAdj)
colnames(VOOMWEIGHTS_LOG_EXPRESSION_MAT_Holdout.PATIENT.PostLibSizeAdj)<-uniqPatIDs.Holdout.PostLibSizeAdj
for(patID in uniqPatIDs.Holdout.PostLibSizeAdj){
	if(sum(patIDs.Holdout.PostLibSizeAdj==patID)>1){
		VOOM_RAW_LOG_EXPRESSION_MAT_Holdout.PATIENT.PostLibSizeAdj[,colnames(VOOM_RAW_LOG_EXPRESSION_MAT_Holdout.PATIENT.PostLibSizeAdj)==patID]<-rowMeans(vobjHO$E[,patIDs.Holdout.PostLibSizeAdj==patID])
		VOOMWEIGHTS_LOG_EXPRESSION_MAT_Holdout.PATIENT.PostLibSizeAdj[,colnames(VOOMWEIGHTS_LOG_EXPRESSION_MAT_Holdout.PATIENT.PostLibSizeAdj)==patID]<-rowMeans(PRED.GENE_EXPRESSION_MAT.Holdout$weights[,patIDs.Holdout.PostLibSizeAdj==patID])
	}else{
		VOOM_RAW_LOG_EXPRESSION_MAT_Holdout.PATIENT.PostLibSizeAdj[,colnames(VOOM_RAW_LOG_EXPRESSION_MAT_Holdout.PATIENT.PostLibSizeAdj)==patID]<-vobjHO$E[,patIDs.Holdout.PostLibSizeAdj==patID]
		VOOMWEIGHTS_LOG_EXPRESSION_MAT_Holdout.PATIENT.PostLibSizeAdj[,colnames(VOOMWEIGHTS_LOG_EXPRESSION_MAT_Holdout.PATIENT.PostLibSizeAdj)==patID]<-PRED.GENE_EXPRESSION_MAT.Holdout$weights[,patIDs.Holdout.PostLibSizeAdj==patID]
	}
}
# create vobjHO.PATIENT.PostLibSizeAdj
targets.Tmp<-data.frame(group=rep(1,ncol(VOOM_RAW_LOG_EXPRESSION_MAT_Holdout.PATIENT.PostLibSizeAdj)),lib.size=rep(mean(vobjHO$targets$lib.size),ncol(VOOM_RAW_LOG_EXPRESSION_MAT_Holdout.PATIENT.PostLibSizeAdj)),norm.factors=rep(1,ncol(VOOM_RAW_LOG_EXPRESSION_MAT_Holdout.PATIENT.PostLibSizeAdj))) # since this is after normalisation we assume all library sizes and normalisation factors to be equal
rownames(targets.Tmp)<-colnames(VOOM_RAW_LOG_EXPRESSION_MAT_Holdout.PATIENT.PostLibSizeAdj)
vobjHO.PATIENT.PostLibSizeAdj<-list(genes=list(rownames(VOOM_RAW_LOG_EXPRESSION_MAT_Holdout.PATIENT.PostLibSizeAdj)),targets=targets.Tmp,E=VOOM_RAW_LOG_EXPRESSION_MAT_Holdout.PATIENT.PostLibSizeAdj,weights=VOOMWEIGHTS_LOG_EXPRESSION_MAT_Holdout.PATIENT.PostLibSizeAdj,design=NULL)
vobjHO.PATIENT.PostLibSizeAdj<-as(vobjHO.PATIENT.PostLibSizeAdj,Class="EList")

VOOM_NORMALIZED_LOG_EXPRESSION_MAT.Holdout.PATIENT.PostLibSizeAdj<-vobjHO.PATIENT.PostLibSizeAdj$E
VOOM_WEIGHTS_MAT.Holdout.PATIENT.PostLibSizeAdj = vobjHO.PATIENT.PostLibSizeAdj$weights
colnames(VOOM_WEIGHTS_MAT.Holdout.PATIENT.PostLibSizeAdj)<-colnames(vobjHO.PATIENT.PostLibSizeAdj$E)
fit.Common.Holdout.PATIENT.PostLibSizeAdj.VarPart <- fitVarPartModel( exprObj=VOOM_NORMALIZED_LOG_EXPRESSION_MAT.Holdout.PATIENT.PostLibSizeAdj, ~ (1|BATCH), data=info.Holdout.PATIENT , useWeights=TRUE, weightsMatrix=VOOM_WEIGHTS_MAT.Holdout.PATIENT.PostLibSizeAdj)
Res.KNOWN.AllKNOWN.Holdout.PATIENT.PostLibSizeAdj.VarPart <- residuals( fit.Common.Holdout.PATIENT.PostLibSizeAdj.VarPart )
colnames(Res.KNOWN.AllKNOWN.Holdout.PATIENT.PostLibSizeAdj.VarPart)<-colnames(VOOM_NORMALIZED_LOG_EXPRESSION_MAT.Holdout.PATIENT.PostLibSizeAdj)
fit.Common.Holdout.PATIENT.PostLibSizeAdj=lmFit(Res.KNOWN.AllKNOWN.Holdout.PATIENT.PostLibSizeAdj.VarPart, design.Common.Holdout.PATIENT, weights=VOOM_WEIGHTS_MAT.Holdout.PATIENT.PostLibSizeAdj)
Res.KNOWN.AllKNOWN.Holdout.PATIENT.PostLibSizeAdj = residuals(fit.Common.Holdout.PATIENT.PostLibSizeAdj, Res.KNOWN.AllKNOWN.Holdout.PATIENT.PostLibSizeAdj.VarPart) # no intercept added back in!
save(Res.KNOWN.AllKNOWN.Holdout.PATIENT.PostLibSizeAdj,file=paste("output/Res.KNOWN",adjusted_covariatesHoldout,dataversion,low_expression_cutoff,"Holdout.AvgCountPerPatient.PostLibSizeAdj","RData",sep="."))
write.table(Res.KNOWN.AllKNOWN.Holdout.PATIENT.PostLibSizeAdj,quote=FALSE,sep="\t",col.names=TRUE,row.names=TRUE,file=paste("output/Res.KNOWN",adjusted_covariatesHoldout,dataversion,low_expression_cutoff,"Holdout.AvgCountPerPatient.PostLibSizeAdj","txt",sep="."))
Res.KNOWN.AllKNOWN.Holdout.PATIENT.PostLibSizeAdj.scale=t(scale(t(Res.KNOWN.AllKNOWN.Holdout.PATIENT.PostLibSizeAdj),scale = FALSE, center = TRUE))

# repeat for combined data
PRED.GENE_EXPRESSION_MAT.All<- voom(TRUE.GENE_EXPRESSION_DGELIST_MAT_All.NORM,design=design.Common.All,plot=TRUE) 
VOOM_NORMALIZED_LOG_EXPRESSION_MAT.All = PRED.GENE_EXPRESSION_MAT.All$E
VOOM_WEIGHTS_MAT.All = PRED.GENE_EXPRESSION_MAT.All$weights
colnames(VOOM_WEIGHTS_MAT.All)<-colnames(PRED.GENE_EXPRESSION_MAT.All$E)
fit.Common.All.VarPart <- fitVarPartModel( exprObj=VOOM_NORMALIZED_LOG_EXPRESSION_MAT.All, ~ (1|BATCH) + (1|RNA.method), data=info.All , useWeights=TRUE, weightsMatrix=VOOM_WEIGHTS_MAT.All)
Res.KNOWN.AllKNOWN.All.VarPart <- residuals( fit.Common.All.VarPart )
colnames(Res.KNOWN.AllKNOWN.All.VarPart)<-colnames(VOOM_NORMALIZED_LOG_EXPRESSION_MAT.All)
fit.Common.All=lmFit(Res.KNOWN.AllKNOWN.All.VarPart, design.Common.All, weights=VOOM_WEIGHTS_MAT.All)
Res.KNOWN.AllKNOWN.All = residuals(fit.Common.All, Res.KNOWN.AllKNOWN.All.VarPart) # no intercept added back in!
save(Res.KNOWN.AllKNOWN.All,file=paste("output/Res.KNOWN",adjusted_covariates,dataversion,low_expression_cutoff,"All","RData",sep="."))
write.table(Res.KNOWN.AllKNOWN.All,quote=FALSE,sep="\t",col.names=TRUE,row.names=TRUE,file=paste("output/Res.KNOWN",adjusted_covariates,dataversion,low_expression_cutoff,"All","txt",sep="."))

# get average residuals
patIDs.All<-gsub(pattern="_[0-9]*",replacement="",colnames(Res.KNOWN.AllKNOWN.All))
uniqPatIDs.All<-unique(patIDs.All)
Res.KNOWN.AllKNOWN.All.AvgRes.PATIENT<-matrix(nrow=nrow(Res.KNOWN.AllKNOWN.All),ncol=length(uniqPatIDs.All))
rownames(Res.KNOWN.AllKNOWN.All.AvgRes.PATIENT)<-rownames(Res.KNOWN.AllKNOWN.All)
colnames(Res.KNOWN.AllKNOWN.All.AvgRes.PATIENT)<-uniqPatIDs.All
for(patID in uniqPatIDs.All){
	if(sum(patIDs.All==patID)>1){
		Res.KNOWN.AllKNOWN.All.AvgRes.PATIENT[,colnames(Res.KNOWN.AllKNOWN.All.AvgRes.PATIENT)==patID]<-rowMeans(Res.KNOWN.AllKNOWN.All[,patIDs.All==patID])
	}else{
		Res.KNOWN.AllKNOWN.All.AvgRes.PATIENT[,colnames(Res.KNOWN.AllKNOWN.All.AvgRes.PATIENT)==patID]<-Res.KNOWN.AllKNOWN.All[,patIDs.All==patID]
	}
}
save(Res.KNOWN.AllKNOWN.All.AvgRes.PATIENT,file=paste("output/Res.KNOWN",adjusted_covariates,dataversion,low_expression_cutoff,"All","AvgResPerPatient","RData",sep="."))
write.table(Res.KNOWN.AllKNOWN.All.AvgRes.PATIENT,quote=FALSE,sep="\t",col.names=TRUE,row.names=TRUE,file=paste("output/Res.KNOWN",adjusted_covariates,dataversion,low_expression_cutoff,"All","AvgResPerPatient","txt",sep="."))

# repeat for the combined avg count data
PRED.GENE_EXPRESSION_MAT.All.PATIENT<- voom(TRUE.GENE_EXPRESSION_DGELIST_MAT_All.PATIENT.NORM,design=design.Common.All.PATIENT,plot=TRUE) 
VOOM_NORMALIZED_LOG_EXPRESSION_MAT.All.PATIENT = PRED.GENE_EXPRESSION_MAT.All.PATIENT$E
VOOM_WEIGHTS_MAT.All.PATIENT = PRED.GENE_EXPRESSION_MAT.All.PATIENT$weights
colnames(VOOM_WEIGHTS_MAT.All.PATIENT)<-colnames(PRED.GENE_EXPRESSION_MAT.All.PATIENT$E)
fit.Common.All.PATIENT.VarPart <- fitVarPartModel( exprObj=VOOM_NORMALIZED_LOG_EXPRESSION_MAT.All.PATIENT, ~ (1|BATCH) + (1|RNA.method), data=info.All.PATIENT , useWeights=TRUE, weightsMatrix=VOOM_WEIGHTS_MAT.All.PATIENT)
Res.KNOWN.AllKNOWN.All.PATIENT.VarPart <- residuals( fit.Common.All.PATIENT.VarPart )
colnames(Res.KNOWN.AllKNOWN.All.PATIENT.VarPart)<-colnames(VOOM_NORMALIZED_LOG_EXPRESSION_MAT.All.PATIENT)
fit.Common.All.PATIENT=lmFit(Res.KNOWN.AllKNOWN.All.PATIENT.VarPart, design.Common.All.PATIENT, weights=VOOM_WEIGHTS_MAT.All.PATIENT)
Res.KNOWN.AllKNOWN.All.PATIENT = residuals(fit.Common.All.PATIENT, Res.KNOWN.AllKNOWN.All.PATIENT.VarPart) # no intercept added back in!
save(Res.KNOWN.AllKNOWN.All.PATIENT,file=paste("output/Res.KNOWN",adjusted_covariates,dataversion,low_expression_cutoff,"All","AvgCountPerPatient","RData",sep="."))
write.table(Res.KNOWN.AllKNOWN.All.PATIENT,quote=FALSE,sep="\t",col.names=TRUE,row.names=TRUE,file=paste("output/Res.KNOWN",adjusted_covariates,dataversion,low_expression_cutoff,"All","AvgCountPerPatient","txt",sep="."))

# repeat for the combined avg count post normalisation matrix
# new: create library size & voom normalised matrix of combined data by combining the library size & voomed data from the individual train & holdout samples
# i.e. compute avg values par patient from VOOM_RAW_LOG_EXPRESSION_MAT_All
# NB we also need to compute and average the weight matrix
# i.e. compute avg values par patient from vobjAll$E -- the non-adjusted but library-size and voom normalised data (NB the E matrices are actually the same for vobjAll and PRED.GENE_EXPRESSION_MAT.All)
# NB we also need to compute and average the weight matrix -- for this we will use PRED.GENE_EXPRESSION_MAT.All since its weights are computed using the coefficients for the covariates we will adjust for
if(sum(rownames(vobjAll$E)!=rownames(PRED.GENE_EXPRESSION_MAT.All$E))>0){stop("vobjAll$E and PRED.GENE_EXPRESSION_MAT.All$E have to have samples in the same order!")}
if(sum(colnames(vobjAll$E)!=colnames(PRED.GENE_EXPRESSION_MAT.All$E))>0){stop("vobjAll$E and PRED.GENE_EXPRESSION_MAT.All$E have to have genes in the same order!")}
patIDs.All.PostLibSizeAdj<-gsub(pattern="_[0-9]*",replacement="",colnames(vobjAll$E))
uniqPatIDs.All.PostLibSizeAdj<-unique(patIDs.All.PostLibSizeAdj)
VOOM_RAW_LOG_EXPRESSION_MAT_All.PATIENT.PostLibSizeAdj<-matrix(nrow=nrow(vobjAll$E),ncol=length(uniqPatIDs.All.PostLibSizeAdj))
VOOMWEIGHTS_LOG_EXPRESSION_MAT_All.PATIENT.PostLibSizeAdj<-matrix(nrow=nrow(vobjAll$weights),ncol=length(uniqPatIDs.All.PostLibSizeAdj))
rownames(VOOM_RAW_LOG_EXPRESSION_MAT_All.PATIENT.PostLibSizeAdj)<-rownames(vobjAll$E)
colnames(VOOM_RAW_LOG_EXPRESSION_MAT_All.PATIENT.PostLibSizeAdj)<-uniqPatIDs.All.PostLibSizeAdj
rownames(VOOMWEIGHTS_LOG_EXPRESSION_MAT_All.PATIENT.PostLibSizeAdj)<-rownames(VOOM_RAW_LOG_EXPRESSION_MAT_All.PATIENT.PostLibSizeAdj)
colnames(VOOMWEIGHTS_LOG_EXPRESSION_MAT_All.PATIENT.PostLibSizeAdj)<-uniqPatIDs.All.PostLibSizeAdj
for(patID in uniqPatIDs.All.PostLibSizeAdj){
	if(sum(patIDs.All.PostLibSizeAdj==patID)>1){
		VOOM_RAW_LOG_EXPRESSION_MAT_All.PATIENT.PostLibSizeAdj[,colnames(VOOM_RAW_LOG_EXPRESSION_MAT_All.PATIENT.PostLibSizeAdj)==patID]<-rowMeans(vobjAll$E[,patIDs.All.PostLibSizeAdj==patID])
		VOOMWEIGHTS_LOG_EXPRESSION_MAT_All.PATIENT.PostLibSizeAdj[,colnames(VOOMWEIGHTS_LOG_EXPRESSION_MAT_All.PATIENT.PostLibSizeAdj)==patID]<-rowMeans(vobjAll$weights[,patIDs.All.PostLibSizeAdj==patID])
	}else{
		VOOM_RAW_LOG_EXPRESSION_MAT_All.PATIENT.PostLibSizeAdj[,colnames(VOOM_RAW_LOG_EXPRESSION_MAT_All.PATIENT.PostLibSizeAdj)==patID]<-vobjAll$E[,patIDs.All.PostLibSizeAdj==patID]
		VOOMWEIGHTS_LOG_EXPRESSION_MAT_All.PATIENT.PostLibSizeAdj[,colnames(VOOMWEIGHTS_LOG_EXPRESSION_MAT_All.PATIENT.PostLibSizeAdj)==patID]<-vobjAll$weights[,patIDs.All.PostLibSizeAdj==patID]
	}
}
# create vobjAll.PATIENT.PostLibSizeAdj
targets.Tmp<-data.frame(group=rep(1,ncol(VOOM_RAW_LOG_EXPRESSION_MAT_All.PATIENT.PostLibSizeAdj)),lib.size=rep(mean(vobjAll$targets$lib.size),ncol(VOOM_RAW_LOG_EXPRESSION_MAT_All.PATIENT.PostLibSizeAdj)),norm.factors=rep(1,ncol(VOOM_RAW_LOG_EXPRESSION_MAT_All.PATIENT.PostLibSizeAdj))) # since this is after normalisation we assuem all library sizes and normalisation factors to be equal
rownames(targets.Tmp)<-colnames(VOOM_RAW_LOG_EXPRESSION_MAT_All.PATIENT.PostLibSizeAdj)
vobjAll.PATIENT.PostLibSizeAdj<-list(genes=list(rownames(VOOM_RAW_LOG_EXPRESSION_MAT_All.PATIENT.PostLibSizeAdj)),targets=targets.Tmp,E=VOOM_RAW_LOG_EXPRESSION_MAT_All.PATIENT.PostLibSizeAdj,weights=VOOMWEIGHTS_LOG_EXPRESSION_MAT_All.PATIENT.PostLibSizeAdj,design=NULL)
vobjAll.PATIENT.PostLibSizeAdj<-as(vobjAll.PATIENT.PostLibSizeAdj,Class="EList")

VOOM_NORMALIZED_LOG_EXPRESSION_MAT.All.PATIENT.PostLibSizeAdj<-vobjAll.PATIENT.PostLibSizeAdj$E
VOOM_WEIGHTS_MAT.All.PATIENT.PostLibSizeAdj = vobjAll.PATIENT.PostLibSizeAdj$weights
colnames(VOOM_WEIGHTS_MAT.All.PATIENT.PostLibSizeAdj)<-colnames(vobjAll.PATIENT.PostLibSizeAdj$E)
fit.Common.All.PATIENT.PostLibSizeAdj.VarPart <- fitVarPartModel( exprObj=VOOM_NORMALIZED_LOG_EXPRESSION_MAT.All.PATIENT.PostLibSizeAdj, ~ (1|BATCH) + (1|RNA.method), data=info.All.PATIENT , useWeights=TRUE, weightsMatrix=VOOM_WEIGHTS_MAT.All.PATIENT.PostLibSizeAdj)
Res.KNOWN.AllKNOWN.All.PATIENT.PostLibSizeAdj.VarPart <- residuals( fit.Common.All.PATIENT.PostLibSizeAdj.VarPart )
colnames(Res.KNOWN.AllKNOWN.All.PATIENT.PostLibSizeAdj.VarPart)<-colnames(VOOM_NORMALIZED_LOG_EXPRESSION_MAT.All.PATIENT.PostLibSizeAdj)
fit.Common.All.PATIENT.PostLibSizeAdj=lmFit(Res.KNOWN.AllKNOWN.All.PATIENT.PostLibSizeAdj.VarPart, design.Common.All.PATIENT, weights=VOOM_WEIGHTS_MAT.All.PATIENT.PostLibSizeAdj)
Res.KNOWN.AllKNOWN.All.PATIENT.PostLibSizeAdj = residuals(fit.Common.All.PATIENT.PostLibSizeAdj, Res.KNOWN.AllKNOWN.All.PATIENT.PostLibSizeAdj.VarPart) # no intercept added back in!
save(Res.KNOWN.AllKNOWN.All.PATIENT.PostLibSizeAdj,file=paste("output/Res.KNOWN",adjusted_covariates,dataversion,low_expression_cutoff,"All.AvgCountPerPatient.PostLibSizeAdj","RData",sep="."))
write.table(Res.KNOWN.AllKNOWN.All.PATIENT.PostLibSizeAdj,quote=FALSE,sep="\t",col.names=TRUE,row.names=TRUE,file=paste("output/Res.KNOWN",adjusted_covariates,dataversion,low_expression_cutoff,"All.AvgCountPerPatient.PostLibSizeAdj","txt",sep="."))
Res.KNOWN.AllKNOWN.All.PATIENT.PostLibSizeAdj.scale=t(scale(t(Res.KNOWN.AllKNOWN.All.PATIENT.PostLibSizeAdj),scale = FALSE, center = TRUE))

# just to check how different the separately processed datasets would be
commonGenes.PostAdj<-intersect(rownames(Res.KNOWN.AllKNOWN),rownames(Res.KNOWN.AllKNOWN.Holdout))
Res.KNOWN.AllKNOWN.All.PostAdj<-cbind(Res.KNOWN.AllKNOWN[match(commonGenes.PostAdj,rownames(Res.KNOWN.AllKNOWN)),],Res.KNOWN.AllKNOWN.Holdout[match(commonGenes.PostAdj,rownames(Res.KNOWN.AllKNOWN.Holdout)),])
info.All.PostAdj<-info.All[match(colnames(Res.KNOWN.AllKNOWN.All.PostAdj),info.All$Sample.Name),]
write.table(info.All.PostAdj,sep="\t",row.names=F,col.names=T,quote=F,file="combined_317samples-Train_plus_85samples-Holdout_covariates_IRvsIS_IVAN.txt")

#---------Check PCA again on residuals after adjustment for Batch.RNAmethod.Sex.RACE.RepSrcCell.BMI--------#
for(datSet in c("train","test","combined","combined.postadj")){
	if(datSet=="train"){
		pdf(paste("output/PCA_All.residual",adjusted_covariates,"outlierremoval",dataversion,"pdf",sep="."))
		SampleByVariable=t(Res.KNOWN.AllKNOWN)
		infoTmp<-info
	}else if(datSet=="test"){
		pdf(paste("output/PCA_All.residual",adjusted_covariates,"outlierremoval",dataversion,"Holdout.pdf",sep="."))
		SampleByVariable=t(Res.KNOWN.AllKNOWN.Holdout)
		infoTmp<-info.Holdout
	}else if(datSet=="combined"){
		pdf(paste("output/PCA_All.residual",adjusted_covariates,"outlierremoval",dataversion,"All.pdf",sep="."))
		SampleByVariable=t(Res.KNOWN.AllKNOWN.All)
		infoTmp<-info.All
	}else if(datSet=="combined.postadj"){
		pdf(paste("output/PCA_All.residual",adjusted_covariates,"outlierremoval",dataversion,"All.PostAdj.pdf",sep="."))
		SampleByVariable=t(Res.KNOWN.AllKNOWN.All.PostAdj)
		infoTmp<-info.All.PostAdj
	}
	print(datSet)
	clonename<-rownames(SampleByVariable)
	pca <- prcomp(SampleByVariable, scale=T)
	
	#=======pca-1 vs pca-2=======
	print(ggplot(data.frame(pca$x), aes(x= pca$x[,1], y= pca$x[,2], colour=factor(infoTmp$Reprogramming.Source.Cell), label=clonename))+
		  scale_colour_discrete(guide =FALSE) +geom_point() +geom_text(aes(label=clonename),hjust=0, vjust=0,size=2)+labs(title="Reprogramming.Source.Cell")+ theme(legend.position=c(0,0),legend.justification=c(0,0),legend.box.just="top",legend.background=element_rect(fill=alpha("white",0.4))) )
		
	print(ggplot(data.frame(pca$x), aes(x= pca$x[,1], y= pca$x[,2], colour=factor(infoTmp$BATCH), label=clonename))+
		  scale_colour_discrete(guide = FALSE) +geom_point() +geom_text(aes(label=clonename),hjust=0, vjust=0,size=2)+labs(title="BATCH")+ theme(legend.position=c(0,1),legend.justification=c(0,1),legend.box.just="top",legend.background=element_rect(fill=alpha("white",0.4))) )
		
	print(ggplot(data.frame(pca$x), aes(x= pca$x[,1], y= pca$x[,2], colour=factor(infoTmp$Sample.Name), label=clonename))+
		  scale_colour_discrete(guide = FALSE) +geom_point() +geom_text(aes(label=clonename),hjust=0, vjust=0,size=2)+labs(title="Sample.Name")+ theme(legend.position=c(0,1),legend.justification=c(0,1),legend.box.just="top",legend.background=element_rect(fill=alpha("white",0.4))) )
		  
	print(ggplot(data.frame(pca$x), aes(x= pca$x[,1], y= pca$x[,2], colour=factor(infoTmp$Reprogramming.Batch), label=clonename))+
		  scale_colour_discrete(guide = FALSE) +geom_point() +geom_text(aes(label=clonename),hjust=0, vjust=0,size=2) +labs(title="Reprogramming.Batch")+ theme(legend.position=c(0,1),legend.justification=c(0,1),legend.box.just="top",legend.background=element_rect(fill=alpha("white",0.4))) )
		  
	if(datSet=="train"){print(ggplot(data.frame(pca$x), aes(x= pca$x[,1], y= pca$x[,2], colour=factor(infoTmp$Sendai.Virus.Lot), label=clonename))+
		  scale_colour_discrete(guide = FALSE) +geom_point() +geom_text(aes(label=clonename),hjust=0, vjust=0,size=2) +labs(title="Sendai.Virus.Lot") + theme(legend.position=c(0,1),legend.justification=c(0,1),legend.box.just="top",legend.background=element_rect(fill=alpha("white",0.4)))  )}
		
	print(ggplot(data.frame(pca$x), aes(x= pca$x[,1], y= pca$x[,2], colour=factor(infoTmp$SEX), label=clonename))+
		  scale_colour_discrete(guide = FALSE) +geom_point() +geom_text(aes(label=clonename),hjust=0, vjust=0,size=2) +labs(title="SEX") + theme(legend.position=c(0,1),legend.justification=c(0,1),legend.box.just="top",legend.background=element_rect(fill=alpha("white",0.4))) )
		  
	print(ggplot(data.frame(pca$x), aes(x= pca$x[,1], y= pca$x[,2], colour=factor(infoTmp$RACE), label=clonename))+
		  scale_colour_discrete(guide = FALSE) +geom_point() +geom_text(aes(label=clonename),hjust=0, vjust=0,size=2)+labs(title="RACE")+ theme(legend.position=c(0,1),legend.justification=c(0,1),legend.box.just="top",legend.background=element_rect(fill=alpha("white",0.4)))  )
		  
	print(ggplot(data.frame(pca$x), aes(x= pca$x[,1], y= pca$x[,2], colour=factor(infoTmp$State), label=clonename))+
		  scale_colour_discrete(guide = FALSE) +geom_point() +geom_text(aes(label=clonename),hjust=0, vjust=0,size=2) +labs(title="State")+ theme(legend.position=c(0,1),legend.justification=c(0,1),legend.box.just="top",legend.background=element_rect(fill=alpha("white",0.4)))  )
		  
	print(ggplot(data.frame(pca$x), aes(x= pca$x[,1], y= pca$x[,2], colour=factor(infoTmp$RNA.method), label=clonename))+
		  scale_colour_discrete(guide = FALSE) +geom_point() +geom_text(aes(label=clonename),hjust=0, vjust=0,size=2) +labs(title="RNA.method")+ theme(legend.position=c(0,1),legend.justification=c(0,1),legend.box.just="top",legend.background=element_rect(fill=alpha("white",0.4)))  )
		  
	print(ggplot(data.frame(pca$x), aes(x= pca$x[,1], y= pca$x[,2], colour=factor(infoTmp$Age), label=clonename))+
		  scale_colour_discrete(guide = FALSE) +geom_point() +geom_text(aes(label=clonename),hjust=0, vjust=0,size=2)+labs(title="Age")+ theme(legend.position=c(0,1),legend.justification=c(0,1),legend.box.just="top",legend.background=element_rect(fill=alpha("white",0.4)))  ) 
		  
	if(datSet=="combined" | datSet=="combined.postadj"){print(ggplot(data.frame(pca$x), aes(x= pca$x[,1], y= pca$x[,2], colour=factor(infoTmp$datOrigin), label=clonename))+
		  scale_colour_discrete(guide =FALSE) +geom_point() +geom_text(aes(label=clonename),hjust=0, vjust=0,size=2)+labs(title="datOrigin")+ theme(legend.position=c(0,0),legend.justification=c(0,0),legend.box.just="top",legend.background=element_rect(fill=alpha("white",0.4))) )}
	
	#=========pca-1 vs pca-3======#
	print(ggplot(data.frame(pca$x), aes(x= pca$x[,1], y= pca$x[,3], colour=factor(infoTmp$Reprogramming.Source.Cell), label=clonename))+
		  scale_colour_discrete(guide =FALSE) +geom_point() +geom_text(aes(label=clonename),hjust=0, vjust=0,size=2)+labs(title="Reprogramming.Source.Cell")+ theme(legend.position=c(0,1),legend.justification=c(0,1),legend.box.just="top",legend.background=element_rect(fill=alpha("white",0.4))) )
		
	print(ggplot(data.frame(pca$x), aes(x= pca$x[,1], y= pca$x[,3], colour=factor(infoTmp$BATCH), label=clonename))+
		  scale_colour_discrete(guide = FALSE) +geom_point() +geom_text(aes(label=clonename),hjust=0, vjust=0,size=2)+labs(title="BATCH")+ theme(legend.position=c(0,1),legend.justification=c(0,1),legend.box.just="top",legend.background=element_rect(fill=alpha("white",0.4))) )
		
	print(ggplot(data.frame(pca$x), aes(x= pca$x[,1], y= pca$x[,3], colour=factor(infoTmp$Sample.Name), label=clonename))+
		  scale_colour_discrete(guide = FALSE) +geom_point() +geom_text(aes(label=clonename),hjust=0, vjust=0,size=2)+labs(title="Sample.Name")+ theme(legend.position=c(0,1),legend.justification=c(0,1),legend.box.just="top",legend.background=element_rect(fill=alpha("white",0.4))) )
		  
	print(ggplot(data.frame(pca$x), aes(x= pca$x[,1], y= pca$x[,3], colour=factor(infoTmp$Reprogramming.Batch), label=clonename))+
		  scale_colour_discrete(guide = FALSE) +geom_point() +geom_text(aes(label=clonename),hjust=0, vjust=0,size=2) +labs(title="Reprogramming.Batch")+ theme(legend.position=c(0,1),legend.justification=c(0,1),legend.box.just="top",legend.background=element_rect(fill=alpha("white",0.4))) )
		  
	if(datSet=="train"){print(ggplot(data.frame(pca$x), aes(x= pca$x[,1], y= pca$x[,3], colour=factor(infoTmp$Sendai.Virus.Lot), label=clonename))+
		  scale_colour_discrete(guide = FALSE) +geom_point() +geom_text(aes(label=clonename),hjust=0, vjust=0,size=2) +labs(title="Sendai.Virus.Lot") + theme(legend.position=c(0,1),legend.justification=c(0,1),legend.box.just="top",legend.background=element_rect(fill=alpha("white",0.4)))  )}
		
	print(ggplot(data.frame(pca$x), aes(x= pca$x[,1], y= pca$x[,3], colour=factor(infoTmp$SEX), label=clonename))+
		  scale_colour_discrete(guide = FALSE) +geom_point() +geom_text(aes(label=clonename),hjust=0, vjust=0,size=2) +labs(title="SEX") + theme(legend.position=c(0,1),legend.justification=c(0,1),legend.box.just="top",legend.background=element_rect(fill=alpha("white",0.4))) )
		  		  
	print(ggplot(data.frame(pca$x), aes(x= pca$x[,1], y= pca$x[,3], colour=factor(infoTmp$RACE), label=clonename))+
		  scale_colour_discrete(guide = FALSE) +geom_point() +geom_text(aes(label=clonename),hjust=0, vjust=0,size=2)+labs(title="RACE") + theme(legend.position=c(0,1),legend.justification=c(0,1),legend.box.just="top",legend.background=element_rect(fill=alpha("white",0.4))) )
		  
	print(ggplot(data.frame(pca$x), aes(x= pca$x[,1], y= pca$x[,3], colour=factor(infoTmp$State), label=clonename))+
		  scale_colour_discrete(guide = FALSE) +geom_point() +geom_text(aes(label=clonename),hjust=0, vjust=0,size=2) +labs(title="State") + theme(legend.position=c(0,1),legend.justification=c(0,1),legend.box.just="top",legend.background=element_rect(fill=alpha("white",0.4))) )
		  
	print(ggplot(data.frame(pca$x), aes(x= pca$x[,1], y= pca$x[,3], colour=factor(infoTmp$RNA.method), label=clonename))+
		  scale_colour_discrete(guide = FALSE) +geom_point() +geom_text(aes(label=clonename),hjust=0, vjust=0,size=2) +labs(title="RNA.method") + theme(legend.position=c(0,1),legend.justification=c(0,1),legend.box.just="top",legend.background=element_rect(fill=alpha("white",0.4))) )
		  
	print(ggplot(data.frame(pca$x), aes(x= pca$x[,1], y= pca$x[,3], colour=factor(infoTmp$Age), label=clonename))+
		  scale_colour_discrete(guide = FALSE) +geom_point() +geom_text(aes(label=clonename),hjust=0, vjust=0,size=2)+labs(title="Age") + theme(legend.position=c(0,1),legend.justification=c(0,1),legend.box.just="top",legend.background=element_rect(fill=alpha("white",0.4))) )
		  
	if(datSet=="combined" | datSet=="combined.postadj"){print(ggplot(data.frame(pca$x), aes(x= pca$x[,1], y= pca$x[,3], colour=factor(infoTmp$datOrigin), label=clonename))+
		  scale_colour_discrete(guide =FALSE) +geom_point() +geom_text(aes(label=clonename),hjust=0, vjust=0,size=2)+labs(title="datOrigin")+ theme(legend.position=c(0,0),legend.justification=c(0,0),legend.box.just="top",legend.background=element_rect(fill=alpha("white",0.4))) )}
		  
	#=======pca-2 vs pca-3===========#
	print(ggplot(data.frame(pca$x), aes(x= pca$x[,2], y= pca$x[,3], colour=factor(infoTmp$Reprogramming.Source.Cell), label=clonename))+
		  scale_colour_discrete(guide =FALSE) +geom_point() +geom_text(aes(label=clonename),hjust=0, vjust=0,size=2)+labs(title="Reprogramming.Source.Cell")+ theme(legend.position=c(0,1),legend.justification=c(0,1),legend.box.just="top",legend.background=element_rect(fill=alpha("white",0.4))) )
		
	print(ggplot(data.frame(pca$x), aes(x= pca$x[,2], y= pca$x[,3], colour=factor(infoTmp$BATCH), label=clonename))+
		  scale_colour_discrete(guide = FALSE) +geom_point() +geom_text(aes(label=clonename),hjust=0, vjust=0,size=2)+labs(title="BATCH")+ theme(legend.position=c(0,1),legend.justification=c(0,1),legend.box.just="top",legend.background=element_rect(fill=alpha("white",0.4))) )
		
	print(ggplot(data.frame(pca$x), aes(x= pca$x[,2], y= pca$x[,3], colour=factor(infoTmp$Sample.Name), label=clonename))+
		  scale_colour_discrete(guide = FALSE) +geom_point() +geom_text(aes(label=clonename),hjust=0, vjust=0,size=2)+labs(title="Sample.Name")+ theme(legend.position=c(0,1),legend.justification=c(0,1),legend.box.just="top",legend.background=element_rect(fill=alpha("white",0.4))) )
		  
	print(ggplot(data.frame(pca$x), aes(x= pca$x[,2], y= pca$x[,3], colour=factor(infoTmp$Reprogramming.Batch), label=clonename))+
		  scale_colour_discrete(guide = FALSE) +geom_point() +geom_text(aes(label=clonename),hjust=0, vjust=0,size=2) +labs(title="Reprogramming.Batch")+ theme(legend.position=c(0,1),legend.justification=c(0,1),legend.box.just="top",legend.background=element_rect(fill=alpha("white",0.4))) )
		  
	if(datSet=="train"){print(ggplot(data.frame(pca$x), aes(x= pca$x[,2], y= pca$x[,3], colour=factor(infoTmp$Sendai.Virus.Lot), label=clonename))+
		  scale_colour_discrete(guide = FALSE) +geom_point() +geom_text(aes(label=clonename),hjust=0, vjust=0,size=2) +labs(title="Sendai.Virus.Lot")+ theme(legend.position=c(0,1),legend.justification=c(0,1),legend.box.just="top",legend.background=element_rect(fill=alpha("white",0.4))) )  }
		
	print(ggplot(data.frame(pca$x), aes(x= pca$x[,2], y= pca$x[,3], colour=factor(infoTmp$SEX), label=clonename))+
		  scale_colour_discrete(guide = FALSE) +geom_point() +geom_text(aes(label=clonename),hjust=0, vjust=0,size=2) +labs(title="SEX") + theme(legend.position=c(0,1),legend.justification=c(0,1),legend.box.just="top",legend.background=element_rect(fill=alpha("white",0.4))) )
		  
	print(ggplot(data.frame(pca$x), aes(x= pca$x[,2], y= pca$x[,3], colour=factor(infoTmp$RACE), label=clonename))+
		  scale_colour_discrete(guide = FALSE) +geom_point() +geom_text(aes(label=clonename),hjust=0, vjust=0,size=2)+labs(title="RACE") + theme(legend.position=c(0,1),legend.justification=c(0,1),legend.box.just="top",legend.background=element_rect(fill=alpha("white",0.4))) )
		  
	print(ggplot(data.frame(pca$x), aes(x= pca$x[,2], y= pca$x[,3], colour=factor(infoTmp$State), label=clonename))+
		  scale_colour_discrete(guide = FALSE) +geom_point() +geom_text(aes(label=clonename),hjust=0, vjust=0,size=2) +labs(title="State")+ theme(legend.position=c(0,1),legend.justification=c(0,1),legend.box.just="top",legend.background=element_rect(fill=alpha("white",0.4)))  )
		  
	print(ggplot(data.frame(pca$x), aes(x= pca$x[,2], y= pca$x[,3], colour=factor(infoTmp$RNA.method), label=clonename))+
		  scale_colour_discrete(guide = FALSE) +geom_point() +geom_text(aes(label=clonename),hjust=0, vjust=0,size=2) +labs(title="RNA.method") + theme(legend.position=c(0,1),legend.justification=c(0,1),legend.box.just="top",legend.background=element_rect(fill=alpha("white",0.4))) )
		  
	print(ggplot(data.frame(pca$x), aes(x= pca$x[,2], y= pca$x[,3], colour=factor(infoTmp$Age), label=clonename))+
		  scale_colour_discrete(guide = FALSE) +geom_point() +geom_text(aes(label=clonename),hjust=0, vjust=0,size=2)+labs(title="Age") + theme(legend.position=c(0,1),legend.justification=c(0,1),legend.box.just="top",legend.background=element_rect(fill=alpha("white",0.4))) )

	if(datSet=="combined" | datSet=="combined.postadj"){print(ggplot(data.frame(pca$x), aes(x= pca$x[,2], y= pca$x[,3], colour=factor(infoTmp$datOrigin), label=clonename))+
		  scale_colour_discrete(guide =FALSE) +geom_point() +geom_text(aes(label=clonename),hjust=0, vjust=0,size=2)+labs(title="datOrigin")+ theme(legend.position=c(0,0),legend.justification=c(0,0),legend.box.just="top",legend.background=element_rect(fill=alpha("white",0.4))) )}

	dev.off()
}


#=============================
# do the DE analysis
#=============================

design.DE<-model.matrix(~State,info)

#rankedResiduals<-apply(Res.KNOWN.AllKNOWN,2,rank)
#zscoreResiduals<-apply(rankedResiduals,2,function(x){qnorm(x/(length(x)+1),0,1)})

fit=lmFit(Res.KNOWN.AllKNOWN, design.DE)
#fit=lmFit(zscoreResiduals, design.DE)
fit <- eBayes(fit)
topSet = topTable(fit, "StateInsulin sensitive", number=nrow(fit))
save(topSet,file=paste("output/Res.KNOWN",adjusted_covariates,dataversion,low_expression_cutoff,"DE_topSet","RData",sep="."))

# # one clone per sample analysis (doing it 100 times to average out some of the randomness)
# DE.OCP.results<-list()
# for(i in 1:100){
	# patIDs<-gsub(pattern="_[0-9]*",replacement="",colnames(Res.KNOWN.AllKNOWN))
	# idxOCP<-integer(0)
	# for(id in unique(patIDs)){
		# if(sum(patIDs==id)>1){
			# idxOCP<-c(idxOCP,sample(size=1,x=which(patIDs==id)))
		# }else{
			# idxOCP<-c(idxOCP,which(patIDs==id))
		# }
	# }
	# Res.KNOWN.AllKNOWN.OCP<-Res.KNOWN.AllKNOWN[,idxOCP]
	# info.OCP<-info[idxOCP,]
	# design.DE.OCP<-model.matrix(~State,info.OCP)
	# fit.OCP=lmFit(Res.KNOWN.AllKNOWN.OCP, design.DE.OCP)
	# fit.OCP <- eBayes(fit.OCP)
	# topSet.OCP = topTable(fit.OCP, "StateInsulin sensitive", number=nrow(fit.OCP))
	# DE.OCP.results[[i]]<-topSet.OCP
# }
# save(DE.OCP.results,file=paste("output/Res.KNOWN",adjusted_covariates,dataversion,low_expression_cutoff,"DE_topSet_OCP_100times","RData",sep="."))
# allPNonAdj0.05DEGenes<-character(0)
# for(i in 1:100){
	# allPNonAdj0.05DEGenes<-c(allPNonAdj0.05DEGenes,rownames(DE.OCP.results[[i]])[DE.OCP.results[[i]]$P.Value<0.05])
# }
# allPNonAdj0.05DEGenes<-sort(decreasing=T,table(allPNonAdj0.05DEGenes))
# consistentDEGenes<-names(allPNonAdj0.05DEGenes)[allPNonAdj0.05DEGenes>=50]
# write.table(consistentDEGenes,row.names=F,col.names=F,quote=F,file=paste("output/Res.KNOWN",adjusted_covariates,dataversion,low_expression_cutoff,"DE_topSet_OCP_100times_consistentDEGenes","txt",sep="."))

# # do permutation based multiple testing correction
# nRandom<-1000
# allPVals<-matrix(nrow=nrow(topSet),ncol=0)
# rownames(allPVals)<-rownames(topSet)
# for(i in 1:nRandom){
	# print(i)
	# if(i %% 10 == 1 & i>1){
		# rownames(allPVals)<-rownames(topSet)
		# colnames(allPVals)<-paste(sep="","perm",1:(i-1))
		# write.table(allPVals,sep="\t",row.names=T,col.names=T,quote=F,file=paste("output/Res.KNOWN",adjusted_covariates,dataversion,low_expression_cutoff,"DE_allPvalues_matrix","txt",sep="."))
	# }
	
	# idxRand<-sample(1:ncol(Res.KNOWN.AllKNOWN),size=ncol(Res.KNOWN.AllKNOWN)/2,replace=F)
	
	# info.rand<-info
	# stateNew<-rep(levels(info$State)[2],nrow(info.rand))
	# stateNew[idxRand]<-levels(info$State)[1]
	# info.rand$State<-factor(stateNew)	
	# info.IRrand<-info.rand[idxRand,]
	# #info.IRrand$RACE<-factor(info.IRrand$RACE) # won't use RACE here
	# info.IRrand$PatientID<-factor(info.IRrand$PatientID)
	# info.ISrand<-info.rand[-idxRand,]
	# #info.ISrand$RACE<-factor(info.ISrand$RACE) # won't use RACE here
	# info.ISrand$PatientID<-factor(info.ISrand$PatientID)
	
	# # now adjust the above residuals for patient ID separately, get the new residuals and add back the intercept
	# Res.KNOWN.AllKNOWN.IRrand<-Res.KNOWN.AllKNOWN[,match(info.IRrand$Sample.Name,colnames(Res.KNOWN.AllKNOWN))]
	# Res.KNOWN.AllKNOWN.ISrand<-Res.KNOWN.AllKNOWN[,match(info.ISrand$Sample.Name,colnames(Res.KNOWN.AllKNOWN))]

	# design.IRrand<-model.matrix(~PatientID,info.IRrand)
	# design.ISrand<-model.matrix(~PatientID,info.ISrand)

	# VOOM_WEIGHTS_MAT.IRrand<-VOOM_WEIGHTS_MAT[,match(info.IRrand$Sample.Name,colnames(VOOM_WEIGHTS_MAT))]
	# VOOM_WEIGHTS_MAT.ISrand<-VOOM_WEIGHTS_MAT[,match(info.ISrand$Sample.Name,colnames(VOOM_WEIGHTS_MAT))]

	# fit.IRrand=lmFit(Res.KNOWN.AllKNOWN.IRrand, design.IRrand, weights=VOOM_WEIGHTS_MAT.IRrand)
	# Res.KNOWN.AllKNOWN.IRrand = residuals(fit.IRrand, Res.KNOWN.AllKNOWN.IRrand) + fit.IRrand$coefficients[,1] # intercept added back in!
	# fit.ISrand=lmFit(Res.KNOWN.AllKNOWN.ISrand, design.ISrand, weights=VOOM_WEIGHTS_MAT.ISrand)
	# Res.KNOWN.AllKNOWN.ISrand = residuals(fit.ISrand, Res.KNOWN.AllKNOWN.ISrand) + fit.ISrand$coefficients[,1] # intercept added back in!

	# Res.KNOWN.AllKNOWN.rand<-cbind(Res.KNOWN.AllKNOWN.IRrand,Res.KNOWN.AllKNOWN.ISrand)
	# Res.KNOWN.AllKNOWN.rand<-Res.KNOWN.AllKNOWN.rand[,match(info.rand$Sample.Name,colnames(Res.KNOWN.AllKNOWN.rand))]
	
	# design.DE.rand<-model.matrix(~State,info.rand)

	# fit.rand=lmFit(Res.KNOWN.AllKNOWN.rand, design.DE.rand)
	# fit.rand <- eBayes(fit.rand)
	# pValsTmp = topTable(fit.rand, "StateInsulin sensitive", number=nrow(fit.rand))
	# pValsTmp<-pValsTmp[match(rownames(topSet),rownames(pValsTmp)),]
	# allPVals<-cbind(allPVals,pValsTmp$P.Value)
# }

# rownames(allPVals)<-rownames(topSet)
# colnames(allPVals)<-paste(sep="","perm",1:nRandom)
# write.table(allPVals,sep="\t",row.names=T,col.names=T,quote=F,file=paste("output/Res.KNOWN",adjusted_covariates,dataversion,low_expression_cutoff,"DE_allPvalues_matrix","txt",sep="."))

# countMoreExtreme<-rep(0,nrow(topSet))
# for(i in 1:nRandom){
	# idxMoreExtreme<-which(allPVals[,i]<=topSet$P.Value)
	# countMoreExtreme[idxMoreExtreme]<-countMoreExtreme[idxMoreExtreme]+1
# }
# adj.permute.Pvalue<-countMoreExtreme/nRandom

# # alternative adjustment
# adj.permute.Pvalue<-numeric(nrow(topSet))
# allPValsVector<-unlist(as.vector(allPVals))
# for(i in 1:nrow(topSet)){
	# adj.permute.Pvalue[i]<-sum(allPValsVector<=topSet$P.Value[i])/length(allPValsVector)
# }

#==============================
# Load HUGO gene mapping information
#==============================

hugo = read.table("hgnc_complete_set.txt", header=TRUE, sep="\t",fill=TRUE,strip.white = FALSE,stringsAsFactors=FALSE, quote='')
hugo.ENSG=hugo$ensembl_gene_id
hugo.SYMBOL=hugo$symbol
hugo.HGNC=hugo$hgnc_id
hugo.ENTREZ=hugo$entrez_id

bgList<-rownames(topSet)
bgList_symbol<-hugo.SYMBOL[match(bgList,hugo.ENSG)]; bgList_symbol<-bgList_symbol[!is.na(bgList_symbol)]
bgList_entrez<-hugo.ENTREZ[match(bgList,hugo.ENSG)]; bgList_entrez<-bgList_entrez[!is.na(bgList_entrez)]
write.table(bgList,row.names=F,col.names=F,quote=F,file=paste(sep="","output/backgroundGeneList.ensembl.",dataversion,".",low_expression_cutoff,".txt"))
write.table(bgList_symbol,row.names=F,col.names=F,quote=F,file=paste(sep="","output/backgroundGeneList.symbol.",dataversion,".",low_expression_cutoff,".txt"))
write.table(bgList_entrez,row.names=F,col.names=F,quote=F,file=paste(sep="","output/backgroundGeneList.entrez.",dataversion,".",low_expression_cutoff,".txt"))

# idxSig<-which(adj.permute.Pvalue<0.05)
idxSig<-which(topSet$adj.P.Val<0.05)
sigDE<-rownames(topSet)[idxSig]
sigDE_symbol<-hugo.SYMBOL[match(sigDE,hugo.ENSG)]; sigDE_symbol<-sigDE_symbol[!is.na(sigDE_symbol)]
sigDE_entrez<-hugo.ENTREZ[match(sigDE,hugo.ENSG)]; sigDE_entrez<-sigDE_entrez[!is.na(sigDE_entrez)]
write.table(sigDE,row.names=F,col.names=F,quote=F,file=paste("output/Res.KNOWN",adjusted_covariates,dataversion,low_expression_cutoff,"DE_sigGenes_Padj0.05_ensembl","txt",sep="."))
write.table(sigDE_symbol,row.names=F,col.names=F,quote=F,file=paste("output/Res.KNOWN",adjusted_covariates,dataversion,low_expression_cutoff,"DE_sigGenes_Padj0.05_symbol","txt",sep="."))
write.table(sigDE_entrez,row.names=F,col.names=F,quote=F,file=paste("output/Res.KNOWN",adjusted_covariates,dataversion,low_expression_cutoff,"DE_sigGenes_Padj0.05_entrez","txt",sep="."))


# idxSig<-which(adj.permute.Pvalue<0.1)
# sigDE<-rownames(topSet)[idxSig]
# sigDE_symbol<-hugo.SYMBOL[match(sigDE,hugo.ENSG)]; sigDE_symbol<-sigDE_symbol[!is.na(sigDE_symbol)]
# sigDE_entrez<-hugo.ENTREZ[match(sigDE,hugo.ENSG)]; sigDE_entrez<-sigDE_entrez[!is.na(sigDE_entrez)]
# write.table(sigDE,row.names=F,col.names=F,quote=F,file=paste("output/Res.KNOWN",adjusted_covariates,dataversion,low_expression_cutoff,"DE_sigGenes_Padj0.1_ensembl","txt",sep="."))
# write.table(sigDE_symbol,row.names=F,col.names=F,quote=F,file=paste("output/Res.KNOWN",adjusted_covariates,dataversion,low_expression_cutoff,"DE_sigGenes_Padj0.1_symbol","txt",sep="."))
# write.table(sigDE_entrez,row.names=F,col.names=F,quote=F,file=paste("output/Res.KNOWN",adjusted_covariates,dataversion,low_expression_cutoff,"DE_sigGenes_Padj0.1_entrez","txt",sep="."))



