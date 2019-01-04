setwd("/home/praveen/Paurush/COSBI/GutProject")
segataAll=read.csv(file="SegataNet.csv", head=T, sep=",")
segataNets=segataAll[,c(1:6)]
colnames(segataNets)      
bodySites=unique(segataNets[,2])

segata_salivaNet=subset(segataAll, Bodysite.1=="Saliva", select=c(Clade.1, Clade.2))
segata_stoolNet=subset(segataAll, Bodysite.1=="Stool", select=c(Clade.1, Clade.2))



mynets.sif=list(segata_stoolNet, segata_salivaNet)
mynets.graph=list()
for(i in 1:length(mynets.sif)){
G = new("graphNEL", nodes=unique(c(as.character(mynets.sif[[i]]$Clade.1), as.character(mynets.sif[[i]]$Clade.2))), edgemode="undirected")
G = addEdge(as.character(mynets.sif[[i]]$Clade.1), as.character(mynets.sif[[i]]$Clade.2), G)
G = removeSelfLoops(G)
mynets.graph[[i]]=G
}
names(mynets.graph)=c("STOOL", "SALIVA")



library(qgraph)
for(i in 1:length(mynets.graph)){
mygraph=mynets.graph[[i]]
myNodes=nodes(mygraph)
network=as(mygraph, "matrix")
net=graph.adjacency(network,mode="undirected", diag=FALSE)
myq=qgraph(network, attributes=T)
cent=centrality(myq)
myComp=data.frame(Nodes=dimnames(network)[[1]], Out_Degree=cent$OutDegree, In_Degree=cent$InDegree, Closeness=cent$Closeness, Betrweenness=cent$Betweenness)
filename=paste("centarlities_", names(mynets.graph)[i], ".csv",sep="")
write.csv(myComp,file=filename, sep="\t")
}


for(i in 1:length(mynets.graph)){
filename=paste("Graph_", names(mynets.graph)[i], ".pdf",sep="")
pdf(filename, height=10, width=10)
plotMyGraph(mynets.graph[[i]])
dev.off()
}


