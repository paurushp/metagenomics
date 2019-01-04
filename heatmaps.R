myImagePlot <- function(x, ...){
 min <- min(x)
     max <- max(x)
     yLabels <- rownames(x)
     xLabels <- colnames(x)
     title <-c()
  # check for additional function arguments
  if( length(list(...)) ){
    Lst <- list(...)
    if( !is.null(Lst$zlim) ){
       min <- Lst$zlim[1]
       max <- Lst$zlim[2]
    }
    if( !is.null(Lst$yLabels) ){
       yLabels <- c(Lst$yLabels)
    }
    if( !is.null(Lst$xLabels) ){
       xLabels <- c(Lst$xLabels)
    }
    if( !is.null(Lst$title) ){
       title <- Lst$title
    }
  }
# check for null values
if( is.null(xLabels) ){
   xLabels <- c(1:ncol(x))
}
if( is.null(yLabels) ){
   yLabels <- c(1:nrow(x))
}
layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(4,1), heights=c(1,1))
 # Red and green range from 0 to 1 while Blue ranges from 1 to 0
ColorRamp <- rgb( seq(0.95,0.99,length=128),  # Red
                   seq(0.95,0.05,length=128),  # Green
                   seq(0.95,0.05,length=128))  # Blue
ColorLevels <- seq(min, max, length=length(ColorRamp))
 # Reverse Y axis
 reverse <- nrow(x) : 1
 yLabels <- yLabels[reverse]
 x <- x[reverse,]
 # Data Map
 par(mar = c(3,5,2.5,2))
 image(1:length(xLabels), 1:length(yLabels), t(x), col=ColorRamp, xlab="", ylab="", axes=FALSE, zlim=c(min,max))
 if( !is.null(title) ){
    title(main=title)
 }
axis(BELOW<-1, at=1:length(xLabels), labels=xLabels, cex.axis=0.7, las=3)
axis(LEFT <-2, at=1:length(yLabels), labels=yLabels, las= HORIZONTAL<-1, cex.axis=0.7)
 # Color Scale
 par(mar = c(3,2.5,2.5,2))
 image(1, ColorLevels, matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1), col=ColorRamp, xlab="",ylab="", xaxt="n")
 layout(1)
}



myImagePlot(sp)

library(EMA)
cl.sample1=clustering(sp1, method="ward", metric="pearson")
cl.gene1=clustering(t(sp1), method="ward", metric="pearson")
clustering.plot(tree=cl.sample1, tree.sup=cl.gene1, data=sp1, names.supp=TRUE, lab=data.frame(colnames(sp1)), scale="none", title="Abundance heatmap")



















