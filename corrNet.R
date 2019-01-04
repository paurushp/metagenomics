########################################################
### WGCNA ANALYSIS
########################################################

cor.prob <- function(X, dfr = nrow(X) - 2) {
  R <- cor(X, method="spearman")
  above <- row(R) < col(R)
  r2 <- R[above]^2
  Fstat <- r2 * dfr / (1 - r2)
  R[above] <- 1 - pf(Fstat, 1, dfr)
  cor.mat <- t(R)
  cor.mat[upper.tri(cor.mat)] <- NA
  cor.mat
}

data_corNet=exprs(mySet)
data_corNet.BF=data_corNet[,c(1:6)]
data_corNet.FF=data_corNet[,c(6:12)]
data_corNet.BF=t(data_corNet.BF)
data_corNet.FF=t(data_corNet.FF)
colnames(data_corNet.FF)=Anot$PROBE_NAME
colnames(data_corNet.BF)=Anot$PROBE_NAME

set.seed(123)
data <- matrix(rnorm(100), 20, 5)
cornet_BF=cor.prob(data_corNet.BF)
cornet_FF=cor.prob(data_corNet.FF)

