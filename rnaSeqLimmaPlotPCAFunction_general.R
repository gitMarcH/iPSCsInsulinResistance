library(ggplot2)
library(scales)

doPcaPlot<-function(pcaObject,covInfo,voomMat=NULL,sampleLabels=NULL,covList,labelList,typeList,libGG=TRUE){
	# pcaObject = object returned by prcomp()
	# covInfo = matrix with the covariates as a data frame (rows=samples, cols=covariates)
	# voomMat = expression matrix that was used (as transpose) in the call to prcomp(); only the column / sample labels will be used; if sampleLabels is provided, voomMat can be NULL
	# sampleLabels = column labels of voomMat, resp. vector with sample IDs that match the rows of pcaObject; can be NULL if voomMat is specified or the row names of pcaObj should be used
	# covList = vector with names of covariates 
	# labelList = vector with labels to be used for covariates (obviously can be the same as covList)
	# typeList = vector giving the type of each covariate; ideally restricted to "factor" and "numeric"; otherwise "character" will be treated as "factor" and "integer" and "double as "numeric"
	# libGG = whether library ggplot2 should be used or not; NB this is included for legacy reasons; must be set to TRUE
	
	pve=round(digits=2,100*pcaObject$sdev^2/sum(pcaObject$sdev^2))
	if(ncol(pcaObject$x)>2){
		pcaObj<-data.frame(pc1=pcaObject$x[,1],pc2=pcaObject$x[,2],pc3=pcaObject$x[,3])
	}else{
		pcaObj<-data.frame(pc1=pcaObject$x[,1],pc2=pcaObject$x[,2])
	}
	
	if(is.null(voomMat) & is.null(sampleLabels)){
		sampleLabels<-rownames(pcaObj)
	}else if(is.null(sampleLabels)){
		sampleLabels<-colnames(voomMat)
	}
	
	covList<-unlist(strsplit(split=",",covList))
	labelList<-unlist(strsplit(split=",",labelList))
	typeList<-unlist(strsplit(split=",",typeList))
	
	pcaObj<-data.frame(pcaObj)
	pcaObj$sampleLabels<-rownames(pcaObj)
	
	if(length(covList)!=length(labelList) | length(covList)!=length(typeList)){stop("lists of covariates, labels and types need to be of same length!")}

	for(i in 1:length(covList)){
		covName<-covList[i]
				
		if(libGG){
			print(paste(sep="",i," of ",length(covList),"; name: ",covName,", type: ",typeList[i]))
					
			if(typeList[i]=="factor" | typeList[i]=="character"){		
				pcaObj$col<-factor(gsub(as.character(covInfo[[covName]]),pattern=" ",replacement="_"))
				
				# PC1 vs PC2
				print(ggplot(pcaObj, aes(x=pc1, y=pc2, colour=col, label=sampleLabels))
					+ geom_point()
					+ geom_text(aes(label=sampleLabels),hjust=0,vjust=0,size=2)
					+ labs(title=labelList[i],colour=covName)
					+ xlab(label=paste(sep="","PC1 (PVE ",pve[1],"%)"))
					+ ylab(label=paste(sep="","PC2 (PVE ",pve[2],"%)"))
					+ theme(legend.position=c(0,1),legend.justification=c(0,1),legend.box.just="top",legend.background=element_rect(fill=alpha("white",0.4))) )

				if(ncol(pcaObject$x)>2){						
				# PC1 vs PC3
				print(ggplot(pcaObj, aes(x=pc1, y=pc3, colour=col, label=sampleLabels))
					+ geom_point()
					+ geom_text(aes(label=sampleLabels),hjust=0,vjust=0,size=2)
					+ labs(title=labelList[i],colour=covName)
					+ xlab(label=paste(sep="","PC1 (PVE ",pve[1],"%)"))
					+ ylab(label=paste(sep="","PC3 (PVE ",pve[3],"%)"))
					+ theme(legend.position=c(0,1),legend.justification=c(0,1),legend.box.just="top",legend.background=element_rect(fill=alpha("white",0.4))) )

				# PC2 vs PC3			
				print(ggplot(pcaObj, aes(x=pc2, y=pc3, colour=col, label=sampleLabels))
					+ geom_point()
					+ geom_text(aes(label=sampleLabels),hjust=0,vjust=0,size=2)
					+ labs(title=labelList[i],colour=covName)
					+ xlab(label=paste(sep="","PC2 (PVE ",pve[2],"%)"))
					+ ylab(label=paste(sep="","PC3 (PVE ",pve[3],"%)"))
					+ theme(legend.position=c(0,1),legend.justification=c(0,1),legend.box.just="top",legend.background=element_rect(fill=alpha("white",0.4))) )
				}					
			}else if(typeList[i]=="numeric" | typeList[i]=="integer" | typeList[i]=="double"){
				pcaObj$col<-as.numeric(covInfo[[covName]])
						
				# PC1 vs PC2			
				print(ggplot(pcaObj, aes(x=pc1, y=pc2, colour=col, label=sampleLabels))
					+ geom_point()
					+ geom_text(aes(label=sampleLabels),hjust=0,vjust=0,size=2)
					+ labs(title=labelList[i],colour=covName)
					+ xlab(label=paste(sep="","PC1 (PVE ",pve[1],"%)"))
					+ ylab(label=paste(sep="","PC2 (PVE ",pve[2],"%)"))
					+ theme(legend.position=c(0,1),legend.justification=c(0,1),legend.box.just="top",legend.background=element_rect(fill=alpha("white",0.4))) )

				if(ncol(pcaObject$x)>2){
				# PC1 vs PC3
				print(ggplot(pcaObj, aes(x=pc1, y=pc3, colour=col, label=sampleLabels))
					+ geom_point()
					+ geom_text(aes(label=sampleLabels),hjust=0,vjust=0,size=2)
					+ labs(title=labelList[i],colour=covName)
					+ xlab(label=paste(sep="","PC1 (PVE ",pve[1],"%)"))
					+ ylab(label=paste(sep="","PC3 (PVE ",pve[3],"%)"))
					+ theme(legend.position=c(0,1),legend.justification=c(0,1),legend.box.just="top",legend.background=element_rect(fill=alpha("white",0.4))) )
				
				# PC2 vs PC3
				print(ggplot(pcaObj, aes(x=pc2, y=pc3, colour=col, label=sampleLabels))
					+ geom_point()
					+ geom_text(aes(label=sampleLabels),hjust=0,vjust=0,size=2)
					+ labs(title=labelList[i],colour=covName)
					+ xlab(label=paste(sep="","PC2 (PVE ",pve[2],"%)"))
					+ ylab(label=paste(sep="","PC3 (PVE ",pve[3],"%)"))
					+ theme(legend.position=c(0,1),legend.justification=c(0,1),legend.box.just="top",legend.background=element_rect(fill=alpha("white",0.4))) )
				}
			}
		}
	}
}
