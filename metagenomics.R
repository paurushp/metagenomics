#####################################################################################
## Description: Code to analyse Metagenomic data				    #
## Author: Paurush Praveen
## Contact: praveen@cosbi.eu							    #
## Data set used: ERS075404-ERS075415 						    #
## Platform: Data biom files from EBI Metagenomics			 	    #
## Biology: 16s sequences from infant infant stool samples with two 	  	    #
## diets 1. Breast Fed and 2. Formula Fed. Formula fed assumed as control	    #
#####################################################################################


### Diversity Analysis


DIR="/home/praveen/Paurush/COSBI/D6/Data/MetaGSE31075/metagenomicAnalysis/output/"

allabundance=read.csv(paste(DIR,"Onlyspecies.csv", sep=""), head=TRUE, sep="\t")
allabundanceHighTaxa=read.csv("/home/praveen/Software/metaphlan/output/merged_abundance_table.txt", head=TRUE, sep="\t")
head(allabundance)
allabundance=t(allabundance)
columns=allabundance$ID
allabundance_BF=allabundance[,c(2:7)]
allabundance_FF=allabundance[,c(8:13)]
allabundance_BF=t(allabundance_BF)
allabundance_FF=t(allabundance_FF)
colnames(allabundance_BF)=columns
colnames(allabundance_FF)=columns
allabundance_BF=data.frame(allabundance_BF)
allabundance_FF=data.frame(allabundance_FF)
diversity(allabundance_FF, "inv")


###################################################################
# Outcome
# > diversity(allabundance_FF, "shannon")
#      FF_1      FF_2      FF_3      FF_4      FF_5      FF_6 
# 1.0812252 1.1268256 0.1837957 1.2461021 0.8314139 1.6453186 
# > diversity(allabundance_BF, "shannon")
#      BF_1      BF_2      BF_3      BF_4      BF_5      BF_6 
# 1.0653076 0.7659192 1.7043706 0.8981299 0.5813619 0.7965836 


# > diversity(allabundanceHighTaxa_FF, "shannon")
#     FF_1     FF_2     FF_3     FF_4     FF_5     FF_6 
# 2.775818 2.417726 2.078146 2.950417 2.478781 3.106067 
# > diversity(allabundanceHighTaxa_BF, "shannon")
#     BF_1     BF_2     BF_3     BF_4     BF_5     BF_6 
# 2.281645 2.518249 2.925836 2.357567 2.173701 2.566292 
###################################################################


### Plot Diversity

shan=read.csv("/home/praveen/Paurush/COSBI/D6/Data/MetaGSE31075/metagenomicAnalysis/diversityAnalysis.txt", head=TRUE, sep="\t")
qplot(x= Sample, y = Shannon, data = shan, color = Feeding, geom  = c('point', 'smooth'), facets= ~Taxonomy,shape = factor(Feeding), ylab  = "Shannon Index", xlab  = "Samples")+ theme_bw()

BCI=allabundance_FF
BCI=BCI*100
BCI=round(BCI,0)
     H2 <- diversity(BCI)
     simp <- diversity(BCI, "simpson")
     invsimp <- diversity(BCI, "inv")
     ## Unbiased Simpson of Hurlbert 1971 (eq. 5):
     unbias.simp <- rarefy(BCI, 2) - 1
     ## Fisher alpha
     alpha <- fisher.alpha(BCI)
     ## Plot all
     pairs(cbind(H, simp, invsimp, unbias.simp, alpha), pch="+", col="red")
     ## Species richness (S) and Pielou's evenness (J):
     S <- specnumber(BCI) ## rowSums(BCI > 0) does the same...
     J <- H/log(S)
     ## beta diversity defined as gamma/alpha - 1:
     data(dune)
     data(dune.env)
     alpha <- with(dune.env, tapply(specnumber(dune), Management, mean))
     gamma <- with(dune.env, specnumber(dune, Management))
     gamma/alpha - 1
     ## Rarefaction
     (raremax <- min(rowSums(BCI)))
     Srare <- rarefy(BCI, raremax)
     plot(S, Srare, col = "red",xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")
     abline(0, 1)
     rarecurve(BCI, step = 20, sample = raremax, col = "red", cex = 0.6)


