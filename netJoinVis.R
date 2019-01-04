### Get all the networks together
##############################################################################################
# Object Details
G_Mu #The graph of relation between microbes and genes (both conditions)
NEL_BF #The gene coexpression network from data (Breast Fed conditions) 
NEL_FF #The gene coexpression network from data (Formula Fed conditions)
brayGraph_FF #Species network from abundance data condition Formula Fed-BRAY CURTIS
brayGraph_BF #Species network from abundance data condition Breast Fed-BRAY CURTIS
nel.cornet_FF #Species network from abundance data condition Formula Fed-CORRELATION NETWORK
nel.cornet_BF #Species network from abundance data condition Breast Fed-CORRELATION NETWORK
##############################################################################################
combine into one object
MyallGraphs=list(Species_Gene=G_Mu, Gene_GeneBF=NEL_BF, Gene_GeneFF=NEL_FF, Species_BrayFF=brayGraph_FF, Species_BrayBF=brayGraph_BF, Species_CorFF=nel.cornet_FF, Species_CorBF=nel.cornet_BF)

# > MyallGraphs
# $Species_Gene
# A graphNEL graph with directed edges
# Number of Nodes = 843 
# Number of Edges = 418 
# 
# $Gene_GeneBF
# A graphNEL graph with undirected edges
# Number of Nodes = 2172 
# Number of Edges = 23326 
# 
# $Gene_GeneFF
# A graphNEL graph with undirected edges
# Number of Nodes = 2172 
# Number of Edges = 16553 
# 
# $Species_BrayFF
# A graphNEL graph with undirected edges
# Number of Nodes = 35 
# Number of Edges = 98 
# 
# $Species_BrayBF
# A graphNEL graph with undirected edges
# Number of Nodes = 35 
# Number of Edges = 36 
# 
# $Species_CorFF
# A graphNEL graph with directed edges
# Number of Nodes = 35 
# Number of Edges = 61 
# 
# $Species_CorBF
# A graphNEL graph with directed edges
# Number of Nodes = 35 
# Number of Edges = 20 

MyallGraphs=list(Species_Gene=G_Mu, Gene_GeneBF=NEL_BF, Gene_GeneFF=NEL_FF, Species_BrayFF=brayGraph_FF, Species_BrayBF=brayGraph_BF, Species_CorFF=nel.cornet_FF, Species_CorBF=nel.cornet_BF)
save(MyallGraphs, file="allGraphs.rda")

save(MyallGraphs, file="allGraphs.rda")




brayBasedjoinBF=join(MyallGraphs$Species_Gene,  MyallGraphs$Gene_GeneBF)
brayBasedjoinBF=join(brayBasedjoinBF, MyallGraphs$Species_BrayBF)

brayBasedjoinFF=join(MyallGraphs$Species_Gene,  MyallGraphs$Gene_GeneFF)
brayBasedjoinFF=join(brayBasedjoinFF, MyallGraphs$Species_BrayFF)

#######################################################################
#######################################################################
graphNEL2SIF=function(G){
  sif=data.frame()
  myG=as(G, "matrix")
  myG[myG>0]<-1
  for(i in 1:(nrow(myG)-1)){
    n=i+1
    for(j in n:ncol(myG)){
      if(myG[i,j]==1)
	sif=rbind(sif, data.frame(Node1=rownames(myG)[i], Node2=colnames(myG)[j]))
    }
  }
return(sif)
}

brayBasedjoinBF.SIF=graphNEL2SIF(brayBasedjoinBF)
brayBasedjoinFF.SIF=graphNEL2SIF(brayBasedjoinFF)

for(i in 1:length(MyallGraphs)){
  file=paste(names(MyallGraphs)[i], ".csv", sep="")
    write.csv(graphNEL2SIF(MyallGraphs[[i]]), file=file)
    }

write.csv(brayBasedjoinBF.SIF, file="brayBasedjoinBF_SIF.csv")
write.csv(brayBasedjoinFF.SIF, file="brayBasedjoinFF_SIF.csv")

### Compute edge densities
D=c()
for(i in 1:length(MyallGraphs)){
  nd=length(nodes(MyallGraphs[[i]]))
  m=as(MyallGraphs[[i]], "matrix")
  m[m>0]<-1
  ed=sum(m)
  ed=ed/2
  d1=2*ed
  d2=nd*(nd-1)
  D=c(D,(d1/d2))
  rm(nd,ed, d1, d2)
}
  names(D)=names(MyallGraphs)
###########################################################################
# CorBasedjoinBF=join(MyallGraphs$Species_Gene,  MyallGraphs$Gene_GeneBF)
# CorBasedjoinBF=join(CorBasedjoinBF, MyallGraphs$Species_CorBF)
# 
# CorBasedjoinFF=join(MyallGraphs$Species_Gene,  MyallGraphs$Gene_GeneFF)
# CorBasedjoinFF=join(CorBasedjoinFF, MyallGraphs$Species_CorFF)

############################################################################
# install devtools
install.packages("devtools")
 
# load devtools
library(devtools)
 
# install arcdiagram
install_github('arcdiagram', username='gastonstat')
 
# load arcdiagram
library(arcdiagram)



# location of 'gml' file
mis_file = "/Users/gaston/lesmiserables.txt"
 
# read 'gml' file
mis_graph = read.graph(mis_file, format="gml")


 get edgelist
edgelist = get.edgelist(mis_graph)
 
# get vertex labels
vlabels = get.vertex.attribute(mis_graph, "label")
 
# get vertex groups
vgroups = get.vertex.attribute(mis_graph, "group")
 
# get vertex fill color
vfill = get.vertex.attribute(mis_graph, "fill")
 
# get vertex border color
vborders = get.vertex.attribute(mis_graph, "border")
 
# get vertex degree
degrees = degree(mis_graph)
 
# get edges value
values = get.edge.attribute(mis_graph, "value")

# load reshape
library(reshape)
 
# data frame with vgroups, degree, vlabels and ind
x = data.frame(vgroups, degrees, vlabels, ind=1:vcount(mis_graph))
 
# arranging by vgroups and degrees
y = arrange(x, desc(vgroups), desc(degrees))
 
# get ordering 'ind'
new_ord = y$ind


arcplot(edgelist, ordering=new_ord, labels=vlabels, cex.labels=0.8,
        show.nodes=TRUE, col.nodes=vborders, bg.nodes=vfill,
        cex.nodes = log(degrees)+0.5, pch.nodes=21,
        lwd.nodes = 2, line=-0.5,
        col.arcs = hsv(0, 0, 0.2, 0.25), lwd.arcs = 1.5 * values)




ggplot(melt(as(MyallGraphs$Gene_GeneFF, "matrix")), aes(X1, X2, fill = value)) + geom_tile() + scale_fill_gradient(low = "blue",  high = "yellow")




####################################################
library(network)
library(ggplot2)
library(sna)
library(ergm)
 
 
plotg <- function(net, value=NULL) {
    m <- as.matrix.network.adjacency(net) # get sociomatrix
    # get coordinates from Fruchterman and Reingold's force-directed placement algorithm.
    plotcord <- data.frame(gplot.layout.fruchtermanreingold(m, NULL))
    # or get it them from Kamada-Kawai's algorithm:
    # plotcord <- data.frame(gplot.layout.kamadakawai(m, NULL))
    colnames(plotcord) = c("X1","X2")
    edglist <- as.matrix.network.edgelist(net)
    edges <- data.frame(plotcord[edglist[,1],], plotcord[edglist[,2],])
    plotcord$elements <- as.factor(get.vertex.attribute(net, "elements"))
    colnames(edges) <-  c("X1","Y1","X2","Y2")
    edges$midX  <- (edges$X1 + edges$X2) / 2
    edges$midY  <- (edges$Y1 + edges$Y2) / 2
    pnet <- ggplot()  +
            geom_segment(aes(x=X1, y=Y1, xend = X2, yend = Y2),
                data=edges, size = 0.5, colour="grey") +
            geom_point(aes(X1, X2,colour=elements), data=plotcord) +
            scale_colour_brewer(palette="Set1") +
            scale_x_continuous(breaks = NA) + scale_y_continuous(breaks = NA) +
            # discard default grid + titles in ggplot2
            opts(panel.background = theme_blank()) + opts(legend.position="none")+
            opts(axis.title.x = theme_blank(), axis.title.y = theme_blank()) +
            opts( legend.background = theme_rect(colour = NA)) +
            opts(panel.background = theme_rect(fill = "white", colour = NA)) +
            opts(panel.grid.minor = theme_blank(), panel.grid.major = theme_blank())
    return(print(pnet))
}
 
 
g <- network(150, directed=FALSE, density=0.03)
classes <- rbinom(150,1,0.5) + rbinom(150,1,0.5) + rbinom(150,1,0.5)
set.vertex.attribute(g, "elements", classes)
 
plotg(g)



