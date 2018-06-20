#ref : http://www.sthda.com/english/wiki/correlation-matrix-a-quick-start-guide-to-analyze-format-and-visualize-a-correlation-matrix-using-r-software

y=read.table("intersect_tagClusters_Gm12878_K562_H1hesc_mapped_with_gene_TC_same_dominant_cluster", header = T, sep = "\t")

y = y[order(y$Gm12878),]
row.names(y) = y$Tag_id
y=y[,2:4]
#y_matrix=data.matrix(y)
y.scale=scale(y)

library(Hmisc)
library(corrplot)
m=rcorr(y.scale, type = "spearman") 
            ###OR####
m=rcorr(y.scale, type = "pearson")

corrplot.mixed(m$r) #generate plot

library(PerformanceAnalytics)
chart.Correlation(y, histogram=TRUE, pch=19)


###ScatterPlot analysis###################
res <- cor(y)  #optional
round(res, 2)  #optinal

# ++++++++++++++++++++++++++++
# flattenCorrMatrix
# ++++++++++++++++++++++++++++
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}
library(Hmisc)
res<-rcorr(as.matrix(y))
flattenCorrMatrix(res$r, res$P)

