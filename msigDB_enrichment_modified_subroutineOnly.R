# original code by Noam Beckmann for lists of DE genes
# modified September 2015 by Marc Henrion so as to allow inputting distinct lists of selected genes (pre-determined by user, not determined by code based on p-value) and background genes
# modified December 2016 to translate C5 GO terms to GO IDs using RamiGO

library(HTSanalyzeR)
library(GSEABase)
library(snow)
library(gage)
library(RamiGO)

gsea_enrichment<-function(geneList,backgroundList,list.gs){
    # performs genes set enrichment analysis
    # Input:
    # geneList = list of gene IDs that should be tested for enrichment (typically Entrez IDs, but depends on the IDs used in list.gs)
    # backgroundList = lists of genes to use as background; needs to include all genes in geneList; usually all Human genes is a good guess on what to use here
    # list.gs = a GMT format list from MSigDB (will be read in with readList from library(GAGE))
    # Output:
    # a table of significant enrichment per module
   
    list.gs<-readList(list.gs)
    
    #set cluster for faster processing
    options(cluster=makeCluster(3, "SOCK"))

	#perform enrichment analysis
    gsca<-multiHyperGeoTest(collectionOfGeneSets=list.gs, universe=as.character(backgroundList), hits=as.character(geneList), minGeneSetSize=15, pAdjustMethod = "BH", verbose = TRUE)
	gsca<-cbind(rownames(gsca),gsca)
	colnames(gsca)[1]<-"Term"

	if(length(which(gsca[,"Adjusted.Pvalue"]<1))>0){
	    res.final<-gsca[which(gsca[,"Adjusted.Pvalue"]<1),]
	
	    if (length(res.final)>8){
	        if (exists("res")){
	            res = rbind(res, (res.final))
	        }else{
	            res = (res.final) # creating res table with 1st column module color and rest res.final
	        }
	    }else if (length(res.final)>1 & length(res.final)<9){
	        if (exists("res")){
	            res = rbind(res, (t(res.final)))
	        }else{
	            res = (t(res.final)) # creating res table with 1st column module color and rest res.final
	        }
			colnames(res)[1]="Term"
			if(nrow(res)==0){cat("res is empty")}
	    }
	    
	    #stop cluster after analysis
	    if(is(getOption("cluster"), "cluster")) {
	        stopCluster(getOption("cluster"))
	        options(cluster=NULL)
	    }
	    
	    res.mod = cbind(res, as.numeric(res[,"Observed Hits"])/as.numeric(res[,"Expected Hits"]))# adding a column of ratio of significant over expected (fold enrichment)
	    colnames(res.mod)[c(1,ncol(res.mod))] = c("GO_Term", "fold_enrichment")
	    res.mod = res.mod[,c("GO_Term", "Universe Size", "Gene Set Size", "Total Hits", "Expected Hits", "Observed Hits", "fold_enrichment","Pvalue","Adjusted.Pvalue")]#,"BH")]
	    rownames(res.mod)=c()

	    if(!is.null(dim(res.mod)[1])){
	      data(c5.go.mapping)
	      GO<-c5.go.mapping$goid[match(res.mod[,1],c5.go.mapping$description)]
	      res.mod<-data.frame(Term=res.mod[,1],GO,res.mod[,-1])
	    }
	}else{
		res.mod<-"no enrichment with P_adj<1"
	}
    return(res.mod)
}
