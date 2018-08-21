options(stringsAsFactors=FALSE);

##Load libraries
library(xlsx)
library(WGCNA)


###################################
#HY GO #(loops over all input lists)
###################################

inputpath ="/geschwindlabshares/RNAseq_singlecellfetal/analysis/Enrichment_DS2-11_KeepCC/GO_ELITE/input"    ##Put all input lists in this folder
outputpath = "/geschwindlabshares/RNAseq_singlecellfetal/analysis/Enrichment_DS2-11_KeepCC/GO_ELITE/output"  ##Create this empty folder
backgroundpath = "/geschwindlabshares/RNAseq_singlecellfetal/analysis/Enrichment_DS2-11_KeepCC/GO_ELITE/background"

options(stringsAsFactors=F)
setwd(inputpath)
files = dir()

for(i in 1:length(files)){
  geneInfo=read.table(files[i])
  geneInfo$SystemCode =rep("En",length=nrow(geneInfo))
  colnames(geneInfo) <- c("Source Identifier","SystemCode")
  write.table(geneInfo,paste0("./",files[i]),row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")
}

setwd(backgroundpath)
files = dir()
refbackground=read.table(files)
background <- cbind(refbackground,rep("En",length=length(refbackground)))
colnames(background) <- c("Source Identifier","SystemCode")
write.table(background,paste0("./",files),row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")

## Run GO elite as nohupped shell script:
#codedir <-"/home/jason/bin/GO-Elite_v.1.2.5-Py"
codedir <- "/home/ldelatorre/bin/GO-Elite_v.1.2.5-Py"
nperm=as.integer(10000) # or 50000
system(paste("python ",codedir,"/GO_Elite.py --species Hs --mod Ensembl --permutations ",nperm,"  --method \"z-score\" --zscore 1.96 --pval 0.01 --num 5 --input ",inputpath," --denom "
             ,backgroundpath," --output ",outputpath," &",sep=""))

##Check all genes are in background before running

files = dir(inputpath);
background = read.csv(paste(backgroundpath,"background-proteincoding.txt",sep="/"), sep="\t");
#filecheck = matrix("NA", ncol= length(files));

filecheck = list();

for(i in 1:length(files)){
         thislist = read.csv(paste(inputpath,files[i],sep="/"), sep = "\t");
         matchlist = match(thislist[,1],background$Source.Identifier);
         matchgene = thislist[which(is.na(matchlist)),];
         if (length(matchgene)>0){
         	filecheck[i] = matchgene;
         	#filecheck[1,i] = matchgene;
         	}
         	else {
         		next
         	}
         }
str(filecheck);

## Generate bar graphs:

##Top10 Processes

#Specify filename:
filename="GO_ELITE-wCC-Top10-Filtered.pdf"
lists= substr(files, 1, nchar(files)-4);

pdf(file=paste(outputpath,filename,sep="/"),height=4,width=4);
for (i in 1:length(lists)) {
    ## GO ontology - make a bar plot
    GOoutput = read.table(paste(outputpath,"/GO-Elite_results/CompleteResults/ORA_pruned/",lists[i],'-GO_z-score_elite.txt',sep=""),sep="\t",header=T,quote="");
    ##Take only biological process
    #GOoutput = GOoutput[(GOoutput$Ontology.Type=="biological_process" | GOoutput$Ontology.Type=="molecular_function"),];
    ##Filter large categories (uninterpretable) or few genes (not robust)
    GOoutput = GOoutput[(GOoutput$Number.in.Ontology<=1000 & GOoutput$Number.Changed>=5),];
    if (nrow(GOoutput)>0) {
        par(oma=c(2,10,1,1));
        bp = barplot(GOoutput$Z.Score[10:1],main=paste("GO Ontology:",lists[i]),horiz=TRUE,yaxt='n',col="blue",xlab="Z-score",cex.main=0.5);
        axis(2,at=bp,labels=GOoutput$Ontology.Name[10:1],tick=FALSE,las=2,cex.axis=0.6);
        abline(v=1.96,col="red",lwd=2,lty=1);
    }

    ## KEGG pathway - make a bar plot
    KEGGoutput = read.csv(paste(outputpath,"/GO-Elite_results/CompleteResults/ORA_pruned/",lists[i],'-KEGG_z-score_elite.txt',sep=""),sep="\t",header=T,fill=T);
    ##Filter large categories (uninterpretable) or few genes (not robust)
    KEGGoutput = KEGGoutput[(KEGGoutput$Number.in.Gene.Set<=1000 & KEGGoutput$Number.Changed>=10),];
    if (nrow(KEGGoutput)>0) {
        par(oma=c(2,10,1,1));
        bp = barplot(KEGGoutput$Z.Score[10:1],main=paste("KEGG Ontology:",lists[i]),horiz=TRUE,yaxt='n',col="blue",xlab="Z-score",cex.main=0.5);
        axis(2,at=bp,labels=KEGGoutput$Gene.Set.Name[10:1],tick=FALSE,las=2,cex.axis=0.6);
        abline(v=1.96,col="red",lwd=2,lty=1);
    }

    ## Transcription factors - make a bar plot
    TFoutput = read.table(paste(outputpath,"/GO-Elite_results/CompleteResults/ORA_pruned/",lists[i],'-TFTargets_z-score_elite.txt',sep=""),sep="\t",header=T,fill=T);
    if (nrow(TFoutput)>0) {
        par(oma=c(2,10,1,1));
        bp = barplot(TFoutput$Z.Score[10:1],main=paste("TF Enrichment:",lists[i]),horiz=TRUE,yaxt='n',col="blue",xlab="Z-score",cex.main=0.5);
        axis(2,at=bp,labels=TFoutput$Gene.Set.Name[10:1],tick=FALSE,las=2,cex.axis=0.6);
        abline(v=1.96,col="red",lwd=2,lty=1);
    }
}
dev.off();



##Top50 Processes

#Specify filename:
filename="GO_ELITE-wCC-Top25-Filtered.pdf"
lists= substr(files, 1, nchar(files)-4);

pdf(file=paste(outputpath,filename,sep="/"),height=4,width=4);
for (i in 1:length(lists)) {
    ## GO ontology - make a bar plot
    GOoutput = read.table(paste(outputpath,"/GO-Elite_results/CompleteResults/ORA_pruned/",lists[i],'-GO_z-score_elite.txt',sep=""),sep="\t",header=T,quote="");
    ##Take only biological process
    #GOoutput = GOoutput[(GOoutput$Ontology.Type=="biological_process" | GOoutput$Ontology.Type=="molecular_function"),];
    ##Filter large categories (uninterpretable) or few genes (not robust)
    GOoutput = GOoutput[(GOoutput$Number.in.Ontology<=1000 & GOoutput$Number.Changed>=5),];
    if (nrow(GOoutput)>0) {
        par(oma=c(2,10,1,1));
        bp = barplot(GOoutput$Z.Score[25:1],main=paste("GO Ontology:",lists[i]),horiz=TRUE,yaxt='n',col="blue",xlab="Z-score",cex.main=0.5);
        axis(2,at=bp,labels=GOoutput$Ontology.Name[25:1],tick=FALSE,las=2,cex.axis=0.4, mgp=c(0.5, 0.4, 0));
        abline(v=1.96,col="red",lwd=2,lty=1);
    }

    ## KEGG pathway - make a bar plot
    KEGGoutput = read.csv(paste(outputpath,"/GO-Elite_results/CompleteResults/ORA_pruned/",lists[i],'-KEGG_z-score_elite.txt',sep=""),sep="\t",header=T,fill=T);
    ##Filter large categories (uninterpretable) or few genes (not robust)
    KEGGoutput = KEGGoutput[(KEGGoutput$Number.in.Gene.Set<=1000 & KEGGoutput$Number.Changed>=10),];
    if (nrow(KEGGoutput)>0) {
        par(oma=c(2,10,1,1));
        bp = barplot(KEGGoutput$Z.Score[25:1],main=paste("KEGG Ontology:",lists[i]),horiz=TRUE,yaxt='n',col="blue",xlab="Z-score",cex.main=0.5);
        axis(2,at=bp,labels=KEGGoutput$Gene.Set.Name[25:1],tick=FALSE,las=2,cex.axis=0.6);
        abline(v=1.96,col="red",lwd=2,lty=1);
    }

    ## Transcription factors - make a bar plot
    TFoutput = read.table(paste(outputpath,"/GO-Elite_results/CompleteResults/ORA_pruned/",lists[i],'-TFTargets_z-score_elite.txt',sep=""),sep="\t",header=T,fill=T);
    if (nrow(TFoutput)>0) {
        par(oma=c(2,10,1,1));
        bp = barplot(TFoutput$Z.Score[25:1],main=paste("TF Enrichment:",lists[i]),horiz=TRUE,yaxt='n',col="blue",xlab="Z-score",cex.main=0.5);
        axis(2,at=bp,labels=TFoutput$Gene.Set.Name[25:1],tick=FALSE,las=2,cex.axis=0.6);
        abline(v=1.96,col="red",lwd=2,lty=1);
    }
}
dev.off();
