# input arguments
args<-commandArgs(TRUE)
enrichFile<-args[1]
selModsFile<-args[2]
outPrefix<-args[3]
main<-args[4]

outLog<-paste(sep="",outPrefix,"_R.log")

cat(file=outLog,append=F,paste(sep="","This is plotGoEnrichmentForCoexpModules.R\n\nInput parameters:\n\tenrichFile = < ",enrichFile," >\n\tselModsFile = < ",selModsFile," >\n\toutPrefix = < ",outPrefix," >\n\tmain = < ",main," >.\n\n"))

# read data
enrich<-read.table(sep="\t",quote="",header=T,stringsAsFactors=F,file=enrichFile)
enrich$classic[enrich$classic=="< 1e-30"]<-"1e-30"
enrich$classic<-as.numeric(enrich$classic)
enrich$BH<-as.numeric(enrich$BH)
selMods<-read.table(selModsFile,header=F,stringsAsFactors=F)[,1]

# build matrix with top GO term per module and p-value
mods<-unique(enrich$module)
goModMat<-data.frame(module=mods,GO=rep(NA,length(mods)),p_value=rep(NA,length(mods)),p_value_adj=rep(NA,length(mods)))
for(i in 1:length(mods)){
    mod<-mods[i]
    tmp<-enrich[enrich$module==mod,]
    tmp<-tmp[order(tmp$classic,1/tmp$fold_enrichment),]
    topGOterm<-tmp$Term[1]
    topGOp<-tmp$classic[1]
    topGOpAdj<-tmp$BH[1]
    goModMat$GO[i]<-topGOterm
    goModMat$p_value[i]<-topGOp
    goModMat$p_value_adj[i]<-topGOpAdj
}

goModMat<-goModMat[order(goModMat$p_value),]

#print(goModMat)
write.table(goModMat,sep="\t",row.names=F,col.names=T,quote=F,file=paste(sep="",outPrefix,"_table.tab"))

idxSelMod<-which(is.element(el=goModMat$module,set=selMods))

png(paste(sep="",outPrefix,".png"),width=12,height=8,units="in",res=450)
par(mar=c(12,5,3,1))
tmpPlot<-barplot(height=-log10(goModMat$p_value),names.arg=NA,col=as.character(goModMat$module),ylab=expression(-log[10](P)),cex.axis=0.75,main=main)
abline(col="darkgrey",lty=2,h=-log10(0.05))
par(xpd=T)

for(i in 1:nrow(goModMat)){
	if(!is.element(el=i,set=idxSelMod)){
		text(cex=0.75,x=tmpPlot[i]-0.15,y=par("usr")[1]+0.25,as.character(goModMat$GO[i]),srt=45,adj=c(1,1))
	}else{
		text(cex=0.75,x=tmpPlot[i]-0.15,y=par("usr")[1]+0.25,as.character(goModMat$GO[i]),srt=45,adj=c(1,1),font=2)
	}
#for(mod in selMods){
#    idx.row<-which(goModMat$module==mod)
#    text(x=tmpPlot[idx.row],y=-log10(goModMat$p_value[idx.row])+0.65,"*",cex=1.5,adj=c(0.5,0.5))
#}
}
par(xpd=F)
dev.off()

pdf(paste(sep="",outPrefix,".pdf"),width=12,height=8)
par(mar=c(12,5,3,1))
tmpPlot<-barplot(height=-log10(goModMat$p_value),names.arg=NA,col=as.character(goModMat$module),ylab=expression(-log[10](P)),cex.axis=0.75,main=main)
abline(col="darkgrey",lty=2,h=-log10(0.05))
par(xpd=T)

for(i in 1:nrow(goModMat)){
	if(!is.element(el=i,set=idxSelMod)){
		text(cex=0.75,x=tmpPlot[i]-0.15,y=par("usr")[1]+0.25,as.character(goModMat$GO[i]),srt=45,adj=c(1,1))
	}else{
		text(cex=0.75,x=tmpPlot[i]-0.15,y=par("usr")[1]+0.25,as.character(goModMat$GO[i]),srt=45,adj=c(1,1),font=2)
	}
#for(mod in selMods){
#    idx.row<-which(goModMat$module==mod)
#    text(x=tmpPlot[idx.row],y=-log10(goModMat$p_value[idx.row])+0.65,"*",cex=1.5,adj=c(0.5,0.5))
#}
}
par(xpd=F)
dev.off()
cat(file=outLog,append=T,paste(sep="","This is the end.\n"))
