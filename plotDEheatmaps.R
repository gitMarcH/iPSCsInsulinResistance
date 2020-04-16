library(gplots)
library(RColorBrewer)

# input arguments
args<-commandArgs(TRUE)
#residFile<-args[1]
#deResultsFile<-args[2]
#covFile<-args[3]
#outPrefix<-args[4]
#condense<-args[5]

logFile<-paste(sep="",outPrefix,"_R.log")
pdfFile<-paste(sep="",outPrefix,".pdf")

cat(file=logFile,append=F,paste(sep="","This is plotDEheatmaps.R.\n\nInput parameters:\n\tresidFile = < ",residFile," >\n\tdeResultsFile = < ",deResultsFile," >\n\tcovFile = < ",covFile," >\n\toutPrefix = < ",outPrefix," >\n\tcondense = < ",condense," >.\n\n"))

# load the data
resid<-read.table(sep="\t",header=T,stringsAsFactors=F,residFile)
de<-read.table(sep="\t",header=T,stringsAsFactors=F,deResultsFile)[,1]
cov<-read.table(sep="\t",header=T,stringsAsFactors=F,covFile)

if(condense){
    covTmp<-data.frame(Sample.Name=ifelse(unique(cov$PatientID)=="845A","845A",formatC(as.integer(unique(cov$PatientID)),format="d",width=3,flag="0")),PatientID=unique(cov$PatientID),State=NA)
    for(i in 1:nrow(covTmp)){
        covTmp$State[i]<-cov$State[cov$PatientID==covTmp$PatientID[i]][1]
    }
    cov<-covTmp
}

# extract the DE genes from the residual matrix
dat<-as.matrix(resid[match(de,rownames(resid)),])
datIR<-dat[,is.element(el=gsub(pattern="X",replacement="",colnames(dat)),set=cov$Sample.Name[cov$State=="Insulin resistant"])]
datIS<-dat[,is.element(el=gsub(pattern="X",replacement="",colnames(dat)),set=cov$Sample.Name[cov$State=="Insulin sensitive"])]

# plot the heatmap
dist.pear <- function(x) as.dist(1-cor(t(x)))
hclust.ave <- function(x) hclust(x, method="average")
dendIR<-as.dendrogram(hclust(dist.pear(t(datIR))),method="average")
dendIS<-as.dendrogram(hclust(dist.pear(t(datIS))),method="average")
dend<-merge(dendIR,dendIS)

dat<-dat[,match(c(colnames(datIR)[order.dendrogram(dendIR)],colnames(datIS)[order.dendrogram(dendIS)]),colnames(dat))]

colCols<-ifelse(cov$State[match(gsub(pattern="X",replacement="",colnames(dat)),cov$Sample.Name)]=="Insulin sensitive",brewer.pal(8,"Dark2")[1],brewer.pal(8,"Dark2")[3])
#heatCols<-colorRampPalette(brewer.pal(9, "GnBu"))(200)
blues<-rep(NA,10)
oranges<-rep(NA,10)
tmpBlue<-as.vector(rgb2hsv(col2rgb("blue")))
tmpOrange<-as.vector(rgb2hsv(col2rgb("orange")))
for(i in 1:10){
    blues[i]<-hsv(tmpBlue[1],tmpBlue[2],tmpBlue[3]/sqrt(i))
    oranges[i]<-hsv(tmpOrange[1],tmpOrange[2],tmpOrange[3]/sqrt(i))
}
#heatCols<-colorRampPalette(c("darkblue4","darkblue3","darkblue2","darkblue1","blue","blue","blue","blue","white","orange","orange","orange","orange","darkorange","darkorange2","darkorange3","darkorange4"))(200)
heatCols<-colorRampPalette(c(blues[10:1],"white",oranges))(200)

pdf(paste(sep="",outPrefix,".pdf"),width=10,height=10)
heatmap.2(dat,trace="none",density.info="none",col=heatCols,distfun=dist.pear,hclustfun=hclust.ave,ColSideColors=colCols,labRow=NA,labCol=NA,Colv=dend,key.title=NA,key.xlab="residual expression level",key=F)
dev.off()

png(paste(sep="",outPrefix,".png"),width=10,height=10,res=300,units="in")
heatmap.2(dat,trace="none",density.info="none",col=heatCols,distfun=dist.pear,hclustfun=hclust.ave,ColSideColors=colCols,labRow=NA,labCol=NA,Colv=dend,key.title=NA,key.xlab="residual expression level",key=F)
dev.off()

cat(file=logFile,append=T,"This is the end.\n")

