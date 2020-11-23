###########################
# define useful functions #
###########################

# filtering lowly expressed genes
getGeneFilteredGeneExprMatrix <- function(gene_expression_counts,MIN_GENE_CPM=minGeneCpm,MIN_SAMPLE_PERCENT_WITH_MIN_GENE_CPM=minSamplePercentWithMinGeneCpm){
	        # Make edgeR object:
	        MATRIX.ALL_GENES = DGEList(counts=gene_expression_counts, genes=rownames(gene_expression_counts))
	
	        # Keep genes with at least MIN_GENE_CPM count-per-million reads (cpm) in at least (MIN_SAMPLE_PERCENT_WITH_MIN_GENE_CPM)% of the samples:
	        # Gabriel: version-2: new low-expression gene cutoff: this is due to PCA test on different cutoff to see the batch effect converge
	        fracSamplesWithMinCPM = rowSums(cpm(MATRIX.ALL_GENES) >MIN_GENE_CPM)
	        isNonLowExpr = fracSamplesWithMinCPM >= MIN_SAMPLE_PERCENT_WITH_MIN_GENE_CPM*ncol(gene_expression_counts)        
	        
	        MATRIX.NON_LOW_GENES = MATRIX.ALL_GENES[isNonLowExpr, ]
	        cat(paste("Will normalize expression counts for ", nrow(MATRIX.NON_LOW_GENES), " genes (those with a minimum of ", MIN_GENE_CPM, " CPM in at least ", sprintf("%.2f", 100 * MIN_SAMPLE_PERCENT_WITH_MIN_GENE_CPM), "% of the ", ncol(MATRIX.NON_LOW_GENES), " samples).\n", sep=""))
	
	        return(MATRIX.NON_LOW_GENES)
}

# normalise data
doNorma<-function(expDat,covDat=NULL,form=NULL,blockVar=NULL,saveNorma=FALSE,outPrefix=""){
	# expDat = count / expression data matrix (as a DGEList object!); rows = genes, cols = samples; should have been filtered to the samples and genes relevant for analyses
	# covDat = covariate data matrix with rows = samples matching the columns from expDat
	# form = R formula indicating which covariates to adjust for (NB nothing gets adjusted for here; but covariates are used in the lowess fit of the standard deviation)
	# blockVar = variable that should be used for blocking; if NULL, no blocking will be taken into account
	# outPrefix = path + filename prefix for output files (only used if saveNorm==TRUE)
	# saveNorma = TRUE/FALSE whether or not the normalised data should be written to output or not (a plot will also be produced)

	if(!is.null(form)){if(is.character(form)){form<-as.formula(form)}}

	if(!is.null(covDat) & !is.null(form)){
		designMat<-model.matrix(form,data=covDat)
		newnames=unlist(lapply(strsplit(colnames(designMat)," "),function(x){if(length(x)>1){paste(x,collapse=".")}else{paste(x,sep="")}}))
		colnames(designMat)=newnames
		
		print(dim(designMat))
	}else{
		designMat<-NULL
	}
	
	print(dim(expDat))

	TRUE.GENE_EXPRESSION_DGELIST_MAT.NORM <-calcNormFactors(expDat,method='TMM')
	print(dim(TRUE.GENE_EXPRESSION_DGELIST_MAT.NORM))
	
	if(!is.null(blockVar)){
		vobj<-voom(TRUE.GENE_EXPRESSION_DGELIST_MAT.NORM, design=designMat, plot=F)
		corfit<-duplicateCorrelation(vobj, design=designMat, block=blockVar)
		if(saveNorma){
			pdf(paste(sep="",outPrefix,"_normalisedNotAdjusted.pdf"))
			vobj = voom(TRUE.GENE_EXPRESSION_DGELIST_MAT.NORM, design=designMat, block=blockVar, correlation=corfit$consensus, plot=T)
			dev.off()
			save(vobj,file=paste(sep="",paste(sep="",outPrefix,"_normalisedNotAdjusted.RData")))
		}else{
			vobj<-voom(TRUE.GENE_EXPRESSION_DGELIST_MAT.NORM, design=designMat, block=blockVar, correlation=corfit$consensus, plot=F)
		}
	}else{
		if(saveNorma){
			pdf(paste(sep="",outPrefix,"_normalisedNotAdjusted.pdf"))
			vobj = voom(TRUE.GENE_EXPRESSION_DGELIST_MAT.NORM, design=designMat, plot=T)
			dev.off()
			save(vobj,file=paste(sep="",paste(sep="",outPrefix,"_normalisedNotAdjusted.RData")))
		}else{
			vobj<-voom(TRUE.GENE_EXPRESSION_DGELIST_MAT.NORM, design=designMat, plot=F)
		}
	}
	return(vobj)
}

# computing residuals
getResiduals<-function(expDat,covDat,form,outPrefix,adjCovs=NULL,voomFirst=FALSE,useVarPar=FALSE,blockVar=NULL,returnFit=FALSE,saveNorma=F,dupCorDoneAlready=FALSE){
	# expDat = voom object or count / expression data matrix (the latter only if voom is to be run first); rows = genes, cols = samples; should have been filtered to the samples and genes relevant for analyses
	# covDat = covariate data matrix with rows = samples matching the columns from expDat
	# form = R formula indicating which covariates to adjust for
	# outPrefix = path + filename prefix for output files
	# adjCovs = character string indicating what covariates were adjusted for (only used to include in output filenames)
	# voomFirst = TRUE/FALSE indicating whether voom normalisation should be run first or not
	# useVarPar = TRUE/FALSE indicating whether residuals should be computed using linear mixed models via the variancePartition package or not
	# blockVar = variable that should be used for blocking; if NULL, no blocking will be taken into account; supported only if useVarPart==FALSE
	# returnFit = TRUE / FALSE whether fit object from lmFit() should be returned as well; if TRUE the returned object is a list with 2 elements: $resid for the residuals and $fit for the fit
	# saveNorma = TRUE/FALSE whether or not the normalised data should be written to output or not (a plot will also be produced); only used if voomFirst==TRUE
	# dupCorDoneAlready = TRUE / FALSE whether or not the data have already been normalised by blocking and taking duplicate correlation into account

	if(is.null(adjCovs)){adjCovs<-""}else{adjCovs<-paste(sep="","_",adjCovs)}
	
	if(voomFirst | !useVarPar){
		# build the design matrix
		designMat<-model.matrix(form,data=covDat)
		newnames=unlist(lapply(strsplit(colnames(designMat)," "),function(x){if(length(x)>1){paste(x,collapse=".")}else{paste(x,sep="")}}))
		colnames(designMat)=newnames
	}
	
	if(voomFirst){
		vobj<-doNorma(expDat=expDat,covDat=covDat,form=form,outPrefix=outPrefix,saveNorma=saveNorma,blockVar=NULL) # note that we set blockVar==NULL here; if blockVar is specified non-null, this will be caught later; Gordon Smyth's advice is to run vobj and duplicateCorrelation twice [https://support.bioconductor.org/p/59700/]
	}else{
		vobj<-expDat
	}
	
	if(!useVarPar){
		# use only lmFit to compute residuals
		if(is.null(blockVar)){
			fit<-lmFit(vobj,designMat)
		}else{
			if(!dupCorDoneAlready){
				dupcor<-duplicateCorrelation(vobj,design=designMat,block=blockVar) # this takes quite some time to run; note that Gordon Smyth recommends running both vobj and duplicateCorrelation twice [https://support.bioconductor.org/p/59700/; hence the seemingly duplicated lines below]
				vobj<-voom(vobj,design=designMat,block=blockVar,correlation=dupcor$consensus.correlation)
			}
			dupcor<-duplicateCorrelation(vobj,design=designMat,block=blockVar) # this takes quite some time to run...
			fit<-lmFit(vobj, designMat, block=blockVar, correlation=dupcor$consensus.correlation)
		}
		resid<-residuals(fit,vobj)
		save(resid,file=paste(sep="",outPrefix,"_residuals-lmFit",adjCovs,".RData"))
		res<-resid

		if(returnFit){
			save(fit,file=paste(sep="",outPrefix,"_residuals-lmFit",adjCovs,"_fit.RData"))
			res<-list(resid=res,fit=fit)
		}

	}else{
		# use only variancePartition to compute residuals; NB blocking not supported
		fit.VarPart<-fitVarPartModel(exprObj=vobj, form, data=covDat, useWeights=TRUE)
		resid.VarPart<-residuals(fit.VarPart)
		colnames(resid.VarPart)<-colnames(vobj$E)
                rownames(resid.VarPart)<-rownames(vobj$E)
		save(resid.VarPart,file=paste(sep="",outPrefix,"_residuals-VP",adjCovs,".RData"))
		res<-resid.VarPart
		fit<-fit.VarPart

		if(returnFit){
			save(fit,file=paste(sep="",outPrefix,"_residuals-VP",adjCovs,"_fit.RData"))
			res<-list(resid=res,fit=fit)
		}

	}
		
	return(res)
}

# doing DE analysis
doDE<-function(...,expDat,covDat,form,outPrefix,blockVar=NULL,returnTopSet=FALSE,ENSG=T){
	# ... = expression or character strings which can be parsed to an expression, specifying the desired contrast
	# expDat = voom object (residuals or normalised data); rows = genes, cols = samples
	# covDat = covariate data matrix with rows = samples matching the columns from expDat
	# form = R formula indicating which covariates to adjust for; needs to be intercept-free formula (i.e. starting with ~0+)
	# outPrefix = path + filename prefix for output files
	# blockVar = variable that should be used for blocking; if NULL, no blocking will be taken into account
	# returnTopSet = TRUE/FALSE whether or not to return the table of DE results produced by topTable()
        # ENSG = TRUE / FALSE whether the gene names are in Ensembl ID and should be converted to gene symbols; if TRUE then a HGNC table (object hgncTable) needs to be loaded BEFORE this function is called
    
	designDE<-model.matrix(form,data=covDat)	
	cm<-makeContrasts(...,levels=designDE)
	if(is.null(blockVar)){
		fit.DE<-lmFit(expDat, designDE)
	}else{
		dupcor<-duplicateCorrelation(expDat,design=designDE,block=blockVar) # this takes quite some time to run...
		fit.DE<-lmFit(expDat, designDE, block=blockVar, correlation=dupcor$consensus.correlation)
	}
	fit2.DE<-contrasts.fit(fit.DE,cm)
	fit2.DE<-eBayes(fit2.DE)
	topSetTmp<-topTable(fit2.DE,number=nrow(fit.DE))
	print(length(which(topSetTmp$adj.P.Val<0.05)))

        if(ENSG){
            genesENSG<-rownames(topSetTmp)
            genesSymbol<-hgncTable$symbol[match(genesENSG,hgncTable$ensembl_gene_id)]
            topSetTmp<-data.frame(topSetTmp,geneSymbol=genesSymbol)
        }
	
	write.table(topSetTmp,file=paste(sep="",outPrefix,".txt"),sep="\t",col.names=T,row.names=T,quote=F)
	write.table(topSetTmp[topSetTmp$adj.P.Val<0.01,],file=paste(sep="",outPrefix,"_Padj0.01.txt"),sep="\t",col.names=T,row.names=T,quote=F)
	write.table(topSetTmp[topSetTmp$adj.P.Val<0.05,],file=paste(sep="",outPrefix,"_Padj0.05.txt"),sep="\t",col.names=T,row.names=T,quote=F)
	write.table(topSetTmp[topSetTmp$adj.P.Val<0.1,],file=paste(sep="",outPrefix,"_Padj0.1.txt"),sep="\t",col.names=T,row.names=T,quote=F)
	write.table(rownames(topSetTmp)[topSetTmp$adj.P.Val<0.01],file=paste(sep="",outPrefix,"_geneIDsPAdj0.01.txt"),sep="\t",col.names=F,row.names=F,quote=F)
	write.table(rownames(topSetTmp)[topSetTmp$adj.P.Val<0.05],file=paste(sep="",outPrefix,"_geneIDsPAdj0.05.txt"),sep="\t",col.names=F,row.names=F,quote=F)
	write.table(rownames(topSetTmp)[topSetTmp$adj.P.Val<0.1],file=paste(sep="",outPrefix,"_geneIDsPAdj0.1.txt"),sep="\t",col.names=F,row.names=F,quote=F)
	write.table(rownames(expDat),file=paste(sep="",outPrefix,"_FilteredGenesListForMsigDBBackground.txt"),sep="\t",col.names=F,row.names=F,quote=F)
	
	if(returnTopSet){
		return(topSetTmp)
	}
}

# small helper function
is.factorCovar<-function(x){
	res<-rep(FALSE,length(x))
	if(is.null(dim(x))){
		if(is.factor(x) | is.character(x)){res<-TRUE}
	}else{
		for(j in 1:ncol(x)){
			if(is.factor(x[,j]) | is.character(x[,j])){res[j]<-TRUE}
		}
	}
	return(res)
}

# produce PCA plots
doPcaPlots<-function(expDat,covDat,outFile,covarList=colnames(covDat),covarTypeList=ifelse(is.factorCovar(covDat),"factor","numeric")){
	# expDat = expression data matrix (residuals or normalised data); rows = genes, cols = samples
	# covDat = covariate data matrix with rows = samples matching the columns from expDat
	# outFile = path + filenames of output pdf file
	
	SampleByVariable=t(expDat)
	clonename<-rownames(SampleByVariable)
	pca <- prcomp(SampleByVariable, scale=T)
	pdf(outFile,width=5,height=5)
	doPcaPlot(pcaObject=pca,covInfo=covDat,voomMat=t(expDat),sampleLabels=NULL,covList=covarList,labelList=covarList,typeList=covarTypeList,libGG=TRUE)
	dev.off()
}

# run variance partition
doVarPart<-function(expDat,covDat,form,outPrefix){
	# expDat = voom object (residuals or normalised data); rows = genes, cols = samples
	# covDat = covariate data matrix with rows = samples matching the columns from expDat
	# form = R formula indicating which covariates to adjust for
	# outPrefix = path + filename prefix for output files
	
	if(is.character(form)){form<-as.formula(form)}
	
	results<-fitVarPartModel(expDat, form, covDat)
	varPart<-extractVarPart(results)
	pdf(paste(sep="",outPrefix,".pdf"))
	plotVarPart(varPart)
	dev.off()
	save(varPart,file=paste(sep="",outPrefix,".RData"))
}

# do hierarchical clustering
doHierarchClust<-function(expDat,outFile,colLabels=FALSE,colInfo=NULL){
	# expDat = voom object (residuals or normalised data); rows = genes, cols = samples
	# outFile = path + filenames of output pdf file
        # colLabels = TRUE/FALSE if labels should be colored according to values in colInfo 
	
	dist.res<-dist(t(expDat))
	hclust.normAdj<-hclust(dist.res)

       if(colLabels){
           require(dendextend)
           hclust.normAdj<-as.dendrogram(hclust.normAdj)
           labCols<-as.numeric(colInfo)[order.dendrogram(hclust.normAdj)]
           labels_colors(hclust.normAdj)<-labCols
       }
    
	pdf(outFile,width=20,height=10)
	plot(hclust.normAdj,cex=0.6)
	dev.off()
}

# paired t test for comparing two gene expression matrices for the same genes and samples (e.g. in 2 different tissues)
# assumes normally distributed gene expression data (just like all tests in functions above), but since the only thing we actually use is the mean of the expression differences, if sample size is large enough, normality will anyway hold at least appoximately
pairedSamplesTestForGeneExpMatrices<-function(G1,G2,tail="twoSided"){
	# G1 = first kxn matrix genes on rows, samples on columns (k = # genes, n = # samples)
	# G2 = first kxn matrix genes on rows, samples on columns (k = # genes, n = # samples)
	# tail = one of "twoSided" (default), "oneSidedLargerThan0", "oneSidedSmallerThan0" depending on whether the test should be two-sided or one-sided and, in the latter case, whether the difference is tested to be larger or smaller than 0
	
	if(nrow(G1)!=nrow(G2)){stop("Both matrices must have the same number of rows/genes/transcripts/probes.")}
	if(ncol(G1)!=ncol(G2)){stop("Both matrices must have the same number of columns/samples.")}
	if(sum(rownames(G1)!=rownames(G2))>0){warning("The rownames of the two expression matrices do not match. Will use the names of the first matrix. You have been warned.")}
	k<-nrow(G1)
	n<-ncol(G1)
	
	D<-G1-G2
	d<-rowMeans(D) # if the data are on the log scale (e.g. log2-cpm as per voom) then this is the average fold change
	sd<-apply(X=D,MARGIN=1,FUN=sd) # could also use function rowSds() from package matrixStats [probably faster]
	
	t<-d/(sd/sqrt(n))
	
	if(tail=="twoSided"){
		p<-2*pt(-abs(t),df=n-1) # more usual formula would be 2*(1-pt(abs(t),df=n-1)); but for very low p-values that formula in R would round them to exactly 0
	}else if(tail=="oneSidedLargerThan0"){
		p<-pt(-t,df=n-1) # as per above, the more usual formula for this would be 1-pt(t,df=n-1)
	}else if(tail=="oneSidedSmallerThan0"){
		p<-pt(t,df=n-1)
	}else{stop("Parameter \"tail\" needs to be one of \"twoSided\", \"oneSidedLargerThan0\", \"oneSidedSmallerThan0\".")}
	
	p.adj.BH<-p.adjust(p,method="BH")

	res<-data.frame(id=rownames(G1),meanDiff=d,sdDiff=sd,n=rep(n,length=k),t.stat=t,p.val=p,p.val.adj.BH=p.adj.BH)
	res<-res[order(decreasing=F,res$p.val),]
	return(res)
}

