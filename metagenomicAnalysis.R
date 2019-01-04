
#####################################################################################
## Description: Code to analyse Metagenomic data				      #
## Author: Paurush Praveen							      #
## Data set used: ERS075404-ERS075415 						      #
## Platform: Data biom files from EBI Metagenomics			 	      #
## Biology: 16s sequences from infant infant stool samples with two 	  	      #
## diets 1. Breast Fed and 2. Formula Fed. Formula fed assumed as control	      #
#####################################################################################

DIR="/home/praveen/Paurush/COSBI/D6/Data/MetaGSE31075/AllTaxonomy"
library(phyloseq)
library(biom)
library(ggplot2)
library(scales)
library(grid)
setwd(DIR)
myBiomFiles=list.files(path = ".", pattern = NULL, all.files = FALSE, full.names = FALSE, recursive = FALSE,ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)


plot_richness(mybiom, x = "SEX", color = "SOLUBLE_DIETARY_FIBER_G_AVE") + geom_boxplot()

mybiom=list()
for(i in 1:12){
mybiom[[i]]=import_biom(BIOMfilename=myBiomFiles[i], taxaPrefix=NULL, parallel=FALSE)
}
mybiomBF=merge_phyloseq(mybiom[[1]], mybiom[[2]],mybiom[[3]],mybiom[[4]],mybiom[[5]],mybiom[[6]])
mybiomFF=merge_phyloseq(mybiom[[7]], mybiom[[8]],mybiom[[9]],mybiom[[10]],mybiom[[11]],mybiom[[12]])
mybiomAll=merge_phyloseq(mybiom[[1]], mybiom[[2]],mybiom[[3]],mybiom[[4]],mybiom[[5]],mybiom[[6]], mybiom[[7]], mybiom[[8]],mybiom[[9]],mybiom[[10]],mybiom[[11]],mybiom[[12]])


merge_samples(mybiom[[1]], group, fun=mean)

# Using the annotation files from RTM server
# Preprocess data manually for formatting tab
annot=read.csv("/home/praveen/Paurush/COSBI/D6/Data/gse31075/ERR056992/annotations.txt", sep="\t")
org=annot[,6]
as.character(unique(org))	


###########################################################################
micGenes=read.csv('/home/praveen/Paurush/COSBI/D6/Data/MetaGSE31075/Annotations/Functions/merged_func.txt', na.strings=c("", "NA"), head=TRUE, sep="\t" )
micGenes=micGenes[complete.cases(micGenes),]
> micGenes=read.csv('/home/praveen/Paurush/COSBI/D6/Data/MetaGSE31075/Annotations/Functions/merged_func.txt', na.strings=c("", "NA"), head=TRUE, sep="\t" )
> dim(micGenes)
[1] 93259    13
> micGenes=micGenes[complete.cases(micGenes),]
> dim(micGenes)
[1] 93258    13
> rownames(micGenes)=micGenes[,1]
> micGenes=micGenes[,-1]
> dim(micGenes)
[1] 93258    12
> colnames(micGenes)
 [1] "ERR056989_functions" "ERR056990_functions" "ERR056991_functions"
 [4] "ERR056992_functions" "ERR056993_functions" "ERR056994_functions"
 [7] "ERR056995_functions" "ERR056996_functions" "ERR056997_functions"
[10] "ERR056998_functions" "ERR056999_functions" "ERR057000_functions"
> colnames(micGenes)=c(rep("BF",6), rep("FF",6))

dens=data.matrix(micGenes_hypRm)
library(EMA)
cl.sample = clustering(dens, method="ward", metric="pearson")
cl.gene = clustering(t(dens), method="ward", metric="pearson")
pdf(file=file.path(DIR, "HeatmapMicGenes.pdf"))
clustering.plot(tree=cl.sample, tree.sup=cl.gene, data=dens, names.sup=FALSE, lab=data.frame(KOs), title="log-densities")
dev.off()


p=cumNormStatFast(dens)


library(metagenomeSeq)



######################
cumNormStatFast <-function(obj,pFlag = FALSE,rel=.1,...){
if(class(obj)=="MRexperiment"){
mat = MRcounts(obj,norm=FALSE,log=FALSE)
} else {
stop("Object needs to be a MRexperiment object.")
    }
smat = lapply(1:ncol(mat), function(i) {
sort(mat[which(mat[, i]>0),i], decreasing = TRUE)
})
leng = max(sapply(smat,length))
smat2 = array(NA,dim=c(leng,ncol(mat)))
for(i in 1:ncol(mat)){
smat2[leng:(leng-length(smat[[i]])+1),i] = smat[[i]]
}

rmat2 = sapply(1:ncol(smat2),function(i){
quantile(smat2[,i],p=seq(0,1,length.out=nrow(smat2)),na.rm=TRUE)
})
smat2[is.na(smat2)] = 0
ref1 = rowMeans(smat2)

ncols = ncol(rmat2)
diffr = sapply(1:ncols, function(i) {
ref1 - rmat2[,i]
})
diffr1=rowMedians(abs(diffr))
if(pFlag==TRUE){
plot(abs(diff(diffr1))/diffr1[-1],type="h",...)
abline(h=rel)
axis(1,at=seq(0,length(diffr1),length.out=5),labels = seq(0,1,length.out=5))
}
x= which(abs(diff(diffr1))/diffr1[-1] > rel)[1]/length(diffr1)
if(x<=0.50){
warning("Low quantile estimate. Default value being used.")
x = 0.50
}
obj@expSummary$cumNormStat = x;
return(x)
}
#########################
dataDirectory<-system.file("extdata",package="metagenomeSeq")
taxa =read.delim(file.path(dataDirectory,"CHK_otus.taxonomy.csv"),stringsAsFactors= F)[,2]
otu=read.delim(file.path(dataDirectory,"CHK_otus.taxonomy.csv"),stringsAsFactors= F)[,1]
lung = load_meta(file.path(dataDirectory,"CHK_NAME.otus.count.csv"))
dim(lung$counts)

phenotypeData=as(clin,"AnnotatedDataFrame")
clin=load_phenoData(file.path(dataDirectory,"CHK_clinical.csv"),tran=TRUE)
ord=match(colnames(lung$counts),rownames(clin))
clin=clin[ord, ]
OTUdata=as(lung$taxa,"AnnotatedDataFrame")
varLabels(OTUdata)="taxa"

obj=newMRexperiment(lung$counts,phenoData=phenotypeData,featureData=OTUdata)




myclin=data.frame(ID=colnames(micGenes_hypRm),SampleType=c(rep("BreastFed",6), rep("FormulaFed",6)), SiteSampled=rep("Intestine", 12), Status=c(rep("BreastFed",6), rep("FormulaFed",6)))
mytaxa=rownames(micGenes_hypRm)
myotu=1:nrow(micGenes_hypRm)
dat=list(counts=micGenes_hypRm, taxa=mytaxa)



myphenotypeData=as(myclin,"AnnotatedDataFrame")
myOTUdata=as(data.frame(factor(myotu)),"AnnotatedDataFrame")
varLabels(myOTUdata)="mytaxa"
obj=newMRexperiment(data.matrix(dat$counts),phenoData=NULL,featureData=NULL)

> setwd("/home/praveen/Paurush/COSBI/D6/Data/MetaGSE31075/metagenomicAnalysis")
> save(tabphyla, file="phylumInData.RData")
> pdf("phylums.pdf")
> ggplot(tabphyla, aes(x=Feeding, y=Abundance))+geom_boxplot(aes(fill=Feeding))+facet_wrap(~Phyla)+theme_bw()
> dev.off()







