rm(list = ls())
library(coexpp)
library(corrplot)

options(object.size=200000000)
options(stringsAsFactors = FALSE);

threads=12
coexppSetThreads(threads)

# input
args<-commandArgs(TRUE)
workDir<-args[1] # working directory 
expDat<-args[2] # full path (or relative to workDir) of input expression data (usually residuals); needs to be RData object
sample_name<-args[3] # will appear in the output filenames
date<-args[4] # will appear in the output filenames
beta<-args[5] # (optional); exponent for power adjacency function; if not specified (or set to "FALSE", "NA" or "NULL") the code will determine an appropriate value
do.heatmap<-as.logical(args[6]) # (optional); TRUE/FALSE should the heatmap be plotted
outlier.thr<-as.numeric(args[7]) # (optional); 600 by default; the higher the number the fewer outlying samples will be removed

if(beta=="FALSE" | beta=="NULL" | beta=="NA"){beta<-NULL}else{
	beta<-as.numeric(beta)
	if(is.na(beta) | beta<=0){beta<-NULL}
}	
if(is.na(do.heatmap)){do.heatmap<-FALSE}
if(is.na(outlier.thr)){outlier.thr<-600}

cat(paste(sep="","This is coexppMinervaSpecifiableParameters.R\n\nInput arguments:\n\tworkDir = < ",workDir," >\n\texpDat = < ",expDat," >\n\tsample_name = < ",sample_name," >\n\tdate = < ",date," >\n\tbeta = < ",beta," >\n\tdo.heatmap = < ",do.heatmap," >\n\toutlier.thr = < ",outlier.thr," >.\n\n"))

setwd(workDir)


# define new heatmap plotting function
heatmapCustom <- function(m, colors, plot.raise_power=100, ...) {
        # Following http://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/Simulated-08-Visualization.pdf
        # we set the diagonal of m to NA and raise m to the specified power:
        diag(m) = NA
        m = m ^ plot.raise_power

  heatmap(
    m,
    Rowv=NA, Colv=NA, scale="none", revC=TRUE, symm=TRUE, labRow="", labCol="",   
    ColSideColors=as.character(colors),
    RowSideColors=as.character(colors),
    ...
  )
}

plotTOMHeatmap2<-
function (coexppClusters, geneModules, samplingThreshold = 2000, 
    ...) 
{
    if (length(coexppClusters) > samplingThreshold) {
        warning("Too many observations to construct heatmap directly. Sampling to tractable number of observations.")
        samples <- clusters(coexppClusters)$order[sort(sample(length(coexppClusters), 
            samplingThreshold))]
	samples<-samples[geneModules[samples]!="grey"]
	tomDissSamp<-sampleMatrix(coexppClusters, samples, kind = "tom")
	geneModsSamp<-geneModules[samples]
	save(tomDissSamp,file="coexppPlotTOMDissSamp.RData")
	save(geneModsSamp,file="coexppPlotGeneModsSamp.RData")
	save(samples,file="coexppPlotSamples.RData")
        heatmapCustom(tomDissSamp, 
            geneModsSamp, plot.raise_power=8,...)
    }
    else {
        samples <- clusters(coexppClusters)$order
        heatmapCustom(tom(coexppClusters)[samples, samples], geneModules[samples], 
            ...)
    }
}


#####################################################################
####STEP-1: LOAD IN DATA, Filtering probes, Filtering samples########
#####################################################################

cat(".. reading the data...\n")

allObjs<-ls()
load(expDat)
allObjsNew<-ls()
t<-get(setdiff(allObjsNew,c(allObjs,"allObjs")))
rm(allObjs,allObjsNew)
print(dim(t))

# #ONLY USE FOR GENE CO-EXPP, TURN OFF for PROTEIN-CO-EXPP
# #Threshold for gene
# read.count.threshold = 40 #cut genes with positive values lower than 40 samples.
# genes = apply(t,1,FUN=function(x) { length(x[x > 0])})
# genes = data.frame(cbind(1:length(genes), genes))
# genes = cbind(genes, rownames(t))
# genes = genes[genes[,2] >= read.count.threshold,]
# t = t[genes[,1],]
# print(dim(t))

cat(".. done.\n\n")
#############

# #============================================================#
# # This is only for the raw count RNA-seq data
# # Turn off when analyzing Illumina data or voomed / residual RNA-seq counts/expression
# # filter patients that have less than 10000 read count total
# #============================================================#
# colsums=colSums(t)
# to_remove=which(colsums<10000)
# t = t[,-to_remove]#/colsums[-to_remove]
# # r = resid(linM)
# t = r

############################################################
####STEP-2: Now, we will build coexpression networks########
############################################################
cat(".. doing some reformatting and editing...\n")

output.network.name = paste(sample_name,"_coexpr_network_",date,".txt",sep="")

## STEP 1: replace the x data.frame with your dataset name

##transpose the orignal data matrix (row=probe, col=sample) to fit the 
##format defined by coexpp.
 
dataset.mod = as.data.frame(t(t)); # here we are grabbing just the columns with expr data

names(dataset.mod) = rownames(t)

anno = data.frame(cbind(1:nrow(t)))

anno = cbind(anno, rownames(t))

## STEP 2: decide how many genes to keep by setting the num.genes.to.include variable; if you want to keep all genes, set this variable to ncol(dataset.mod)
## now we need to trim the number genes we are going to look at since we can't handle tens of thousands

num.genes.to.include <- ncol(dataset.mod)
gene.var = apply(dataset.mod, 2, FUN=function(x, percent.trim) { var(sort(x)[floor(percent.trim*length(x)):ceiling((1-percent.trim)*length(x))]) }, percent.trim=0.05) #computes variance per genes after removing highest and lowest percentage from distribution
gene.var <- cbind(1:length(gene.var), gene.var)
gene.var <- gene.var[order(gene.var[,2],decreasing=TRUE),] # ordering genes on variance, highest on top and keeping track of indices
gene.var <- as.vector(unlist(gene.var[1:num.genes.to.include,1])) # keeping only ordered index
dataset.mod <- dataset.mod[,gene.var] # ordering the gene columns by variance, from highest on left to lowest on right
anno.mod = anno[gene.var,] # ordering anno the same way, saving it to anno.mod

# now that we have the data loaded, we test to see which genes have enough data to proceed#

gene.set = goodSamplesGenes(dataset.mod, verbose = 3);

# if there are genes that have too many missing values, then we need to eliminate those from the mix#

if (!gene.set$allOK)
{
    # Optionally, print the gene and sample names that were removed:
    if (sum(!gene.set$goodGenes)>0)
        printFlush(paste("Removing genes:", paste(names(dataset.mod)[!gene.set$goodGenes], collapse = ", ")));

    if (sum(!gene.set$goodSamples)>0)
        printFlush(paste("Removing samples:", paste(rownames(dataset.mod)[!gene.set$goodSamples], collapse = ", ")));

    # Remove the offending genes and samples from the data:
    dataset.mod= dataset.mod[gene.set$goodSamples, gene.set$goodGenes]
}

cat(".. done.\n\n")

############################################################
## STEP 3: generate the cluster tree for the samples 
## to detect outliers; set the height.threshold variable 
## to the height based on the clustergram to include
## samples that fall below the line; if you want all samples included, 
## set this variable to a number that is higher than the graph indicates
##
# now that we have a reasonable data matrix, 
# cluster the experiments to determine if there are any obvious
# outliers
############################################################
cat(".. constructing the cluster tree...\n")

sampleTree = flashClust(dist(dataset.mod), method = "complete");

# based on the clustering, automatically remove outliers
# Plot a line to show the cut; Remove sample above this line
height.threshold <- outlier.thr
pdf(file = paste(sample_name,"_clusterOutlierCheck_",date,".pdf",sep=""), width = 12, height = 9);
sizeGrWindow(12,9)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,cex.axis = 1.5, cex.main = 2)

abline(h = height.threshold, col = "red");  # we want to cut anything above the line
dev.off()
# Determine cluster under the line
#
clust = cutreeStatic(sampleTree, cutHeight = height.threshold, minSize = 5)
table(clust)
# clust 1 contains the samples we want to keep
#
keepSamples = (clust==1)
dataset.mod2 = dataset.mod[keepSamples, ]
nGenes = ncol(dataset.mod2)
nSamples = nrow(dataset.mod2)

cat(".. done.\n\n")

##
## now that we have a filtered dataset of interest, let's go to work building the network
##
######################################################################################
## STEP 4: determine the scale free parameter to choose to ensure an optimal power law distribution is followed for the network construction;
## here you want to pick the number of the graph that plots out below that is at point in the curve just before it levels out; the power.constant
## variable needs to be set to this number
##
# set the dataset to build network from
######################################################################################

cat(".. computing correlations etc...\n")

x = dataset.mod2;

################
#Calculates (correlation or distance) network adjacency from given expression data or from a similarity.
#Correlation and distance are transformed as follows: for type = "unsigned", adjacency = |cor|^power; for type = "signed", adjacency = (0.5 * (1+cor) )^power; for type = "signed hybrid", adjacency = cor^power if cor>0 and 0 otherwise; and for type = "distance", adjacency = (1-(dist/max(dist))^2)^power.
#http://www.inside-r.org/packages/cran/WGCNA/docs/adjacency
#*****Correlation Network****
correlation.gene=adjacency(as.matrix(x),type="unsigned",power=1) #here it calculate absolute value of pearson correlation;
#save(correlation.gene,file=paste(sample_name,"_corr_",date,".RData",sep=""))
#write.table(correlation.gene,paste(sample_name,"_corr_",date,".txt",sep=""),row.names=T,col.names=T)
#*****Correlation Network****

#Now get the row & column index of the corr.>cut-off
corr.cut=0.9#0.5#0.7,0.9
corr.cut.index=which(correlation.gene>corr.cut, arr.ind = T)
corr.cut.index.noauto.index=which(corr.cut.index[,1]!=corr.cut.index[,2])
corr.cut.index.noauto= corr.cut.index[corr.cut.index.noauto.index,]
row.gene.name.corr.noauto=rownames(correlation.gene)[corr.cut.index.noauto[,1]]
col.gene.name.corr.noauto=colnames(correlation.gene)[corr.cut.index.noauto[,2]]

corr.value.cut=c()
for (rw in 1:length(row.gene.name.corr.noauto)){
	
	corr.value.cut= rbind(corr.value.cut,correlation.gene[corr.cut.index.noauto[rw,1], corr.cut.index.noauto[rw,2]])
}
corr.cut.network=cbind(row.gene.name.corr.noauto, col.gene.name.corr.noauto, corr.value.cut)

save(corr.cut.network,file=paste(sample_name,"_corr.cut.network_",corr.cut,"_",date,".RData",sep=""))
write.table(corr.cut.network,paste(sample_name,"_corr.cut.network_",corr.cut,"_",date,".txt",sep=""),row.names=F,col.names=F)

#In cytoscape, directly load in corr.cut.network.txt file, then transparentize edges according to correlation(3rd column); color node based on the module membership

cat(".. done.\n\n")

################

cat(".. building the coexpression network...\n\n")

# define here the number of threadts to be used, in this case, by default, 8.
coexppSetThreads(6)
net = coexpressionAnalysis(as.matrix(x),beta=beta)


table(net$geneModules)
# now get modules for network
#
save(net,file=paste(sample_name,"_coexpr_network_",date,".RData",sep=""))

cat(".. done.\n\n")

if(is.null(beta)){
	cat(".. plotting the power curve and extracting the beta that was used ...\n")
	sft<-net$clusters@sftStatistics
	png(paste(sample_name,"_Threshold_",date,".png",sep=""),width=6,height=5,units="in",res=300)
	cex1 = 0.9;
	# Scale-free topology fit index as a function of the soft-thresholding power
	plot(sft[,1], -sign(sft[,3])*sft[,2],
	xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
	main = paste("Scale independence"));
	text(sft[,1], -sign(sft[,3])*sft[,2],
	labels=sft$Power,cex=cex1,col="red");
	# this line corresponds to using an R^2 cut-off of h
	abline(h=cex1,col="red")
	
	dev.off()
	
	# beta used for the analysis
	write.table(net$clusters@beta,paste(sample_name,"_beta_",date,".txt",sep=""),row.names=F,col.names=F)
	cat(".. done.\n\n")
}

cat(".. plotting cluster dendrogram ...\n")
n <- net
table(n$geneModules);
moduleLabels = n$geneModules;
moduleColors = labels2colors(n$geneModules);
pdf(paste(sample_name,"_clusterDendrogram_",date,".pdf",sep=""),width=12,height=5)
plotClustering(n$clusters,n$geneModules,dendroLabels=FALSE)
dev.off()
cat(".. done.\n\n")


cat(".. computing eigengenes for all modules ...\n")
eigMods<-moduleEigengenes(expr=net$clusters@data,colors=as.character(net$geneModules))$eigengenes
rownames(eigMods)<-rownames(net$clusters@data)
eigMods<-orderMEs(eigMods)
write.table(eigMods,sep="\t",row.names=T,col.names=T,quote=F,file=gsub(output.network.name,pattern=".txt",replacement="_eigenModules.tab"))

pdf(gsub(output.network.name,pattern=".txt",replacement="_eigenModules_correlations.pdf"),width=10,height=10)
corTmp<-cor(eigMods)
corrplot(corTmp,method="circle",type="upper")
dev.off()
cat(".. done.\n\n")

###############
if(do.heatmap){
	cat(".. producing the heat map (may take some time)...\n")
	# run that part only if you want to get the heatmap done, it takes a long time to run, so I am not putting it in by default.
	png(paste(sample_name,"_heatmap_",date,".png",sep=""),height=8,width=8,units="in",res=450)
	plotTOMHeatmap2(n$clusters,n$geneModules)
	#plotTOMHeatmap(n$clusters,n$geneModules,col=colorRampPalette(c("black","brown","darkred","red","yellow","lightyellow"))(12))
	#plotTOMHeatmap(n$clusters,n$geneModules,keep.dendro=T)
	dev.off()
	cat(".. done.\n\n")
}	

cat(".. writing output ...\n")
################
# the output is the network table of 3 columns, 1st column gene name, second to be ingnored and 3rd module color
res <- cbind(moduleLabels, data.frame(names(x)))
res = merge(res, anno, by=2)
colnames(res)<-c("gene","module","toBeIgnored")

write.table(res[,c(1,3,2)], file=output.network.name, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

write.table(res[,1],row.names=F,col.names=F,quote=F,file=gsub(output.network.name,pattern=".txt",replacement="_allGenesBackground.txt"))

for(mod in unique(res[,2])){
	resTmp<-res[res[,2]==mod,]
	write.table(resTmp[,1],row.names=F,col.names=F,quote=F,file=gsub(output.network.name,pattern=".txt",replacement=paste(sep="","_",mod,"_geneList.txt")))
}

cat(".. done.\n\n")


cat("Done.\n\n")
