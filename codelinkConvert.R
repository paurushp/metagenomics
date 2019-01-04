Codelink2CodelinkSet <- function (object, annotation = NULL, phenodata = NULL, featuredata = NULL, intensity = "Smean") 
{
    if (class(object) != "Codelink") 
        stop("Codelink-object needed.")

    switch(intensity, 
      "Smean" = int <- object$Smean,
    "Ri" = int <- object$Ri,
      "Ni" = int <- object$Ni,
    )    
    bkg <- object$Bmedian

    if (is.null(phenodata)) {
      phenodata <- data.frame(sample = object$sample)
      phenodata.varMet <- data.frame(labelDescription = "sample names", row.names = "sample")
      phenodata <- new("AnnotatedDataFrame", data = phenodata, varMetadata = phenodata.varMet)
    }
    
    if (is.null(featuredata)) {
      featuredata <- data.frame(probeName = object$name, probeType = object$type, logicalRow = object$logical[, "row"], logicalCol = object$logical[, "col"], meanSNR = rowMeans(object$snr, na.rm = TRUE))
      featuredata.feMet <- data.frame(labelDescription = c("probe names", "probe types", "probe row position", "probe column position", "mean snr"), row.names = c("probeName", "probeType", "logicalRow",      "logicalCol", "meanSNR"))
    featuredata <- new("AnnotatedDataFrame", data = featuredata, varMetadata = featuredata.feMet)

    }
    
    if (is.null(annotation))
      chip <- annotation(object)
    
    if(is.null(object$weight))
      object$weight=createWeights(object)
    
  tmp <- new("CodelinkSet", exprs = int, background = bkg, 
        flag = object$flag, weight=object$weight, snr = object$snr, annotation = chip)
    phenoData(tmp) <- phenodata
    featureData(tmp) <- featuredata
  experimentData(tmp)@preprocessing <- object$method
  experimentData(tmp)@other <- list("product" = object$product)
    tmp
}
