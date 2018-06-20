##http://www.sthda.com/english/wiki/print.php?id=241

y=read.table("TPM_gene_enhancer_interactions_Allchr_diff_dominant_non_redundent", header = T, sep = "\t")

library(factoextra)
library(cluster)
library(fpc)
library(NbClust)
#y=y[,2:5]
####optinal for removing redundent data with respect of interaction
mn <- pmin(y$Enhancer_TPM, y$Promoter_TPM)
mx <- pmax(y$Enhancer_TPM, y$Promoter_TPM)
int <- as.numeric(interaction(mn, mx))
y_nr = y[match(unique(int), int),] ###End of process to remove redundent data
m = y
y= y_nr
y=unique(y)
y = y[order(y$Gm12878),]
row.names(y) = y$Tag_id
y.scaled=scale(y[,2:4])
# Compute the number of clusters
nb <- NbClust(y.scaled, distance = "euclidean", min.nc = 2, max.nc = 10, method = "complete", index ="all")
# Visualize the result

#library(factoextra)
fviz_nbclust(nb) + theme_minimal()
# K-means clustering
km.res <- eclust(y.scaled, "kmeans", k = 3, nstart = 25, graph = FALSE) #change k according to graph
# k-means group number of each observation
km.res$cluster
# Visualize k-means clusters
fviz_cluster(km.res, geom = "point", ellipse.type = "ellipse")
# PAM clustering
pam.res <- eclust(y.scaled, "pam", k = 3, graph = FALSE)
pam.res$cluster
# Visualize pam clusters
fviz_cluster(pam.res, geom = "point", ellipse.type = "ellipse")
# Enhanced hierarchical clustering
res.hc <- eclust(y.scaled, "hclust", k = 3, method = "complete", graph = FALSE)
head(res.hc$cluster, 15)
# Dendrogram
fviz_dend(res.hc, rect = TRUE, show_labels = FALSE)
# Silhouette coefficient of observations

#library("cluster")
sil <- silhouette(km.res$cluster, dist(y.scaled))
head(sil[, 1:3], 10)
# Silhouette plot
plot(sil, main ="Silhouette plot - K-means")
fviz_silhouette(sil)
# Summary of silhouette analysis
si.sum <- summary(sil)
# Average silhouette width of each cluster
si.sum$clus.avg.widths
# The total average (mean of all individual silhouette widths)
si.sum$avg.width
# The size of each clusters
si.sum$clus.sizes
# Default plot
fviz_silhouette(km.res)
# Change the theme and color
fviz_silhouette(km.res, print.summary = FALSE) + scale_fill_brewer(palette = "Dark2") + scale_color_brewer(palette = "Dark2") + theme_minimal()+theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
# Silhouette information
silinfo <- km.res$silinfo
names(silinfo)
# Silhouette widths of each observation
head(silinfo$widths[, 1:3], 10)
# Average silhouette width of each cluster
silinfo$clus.avg.widths
# The total average (mean of all individual silhouette widths)
silinfo$avg.width
# The size of each clusters
km.res$size
fviz_silhouette(pam.res)
fviz_silhouette(res.hc)
# Silhouette width of observation
sil <- res.hc$silinfo$widths[, 1:3]
# Objects with negative silhouette
neg_sil_index <- which(sil[, 'sil_width'] < 0)
sil[neg_sil_index, , drop = FALSE]
write.table(sil[neg_sil_index, , drop = FALSE] ,file = "negative_cluster_HC", row.names = TRUE, col.names = TRUE)

# Silhouette width of observation for K-means clustering
sil <- km.res$silinfo$widths[, 1:3]
# Objects with negative silhouette
neg_sil_index <- which(sil[, 'sil_width'] < 0)
sil[neg_sil_index, , drop = FALSE]
write.table(sil[neg_sil_index, , drop = FALSE] ,file = "negative_cluster_km", row.names = TRUE, col.names = TRUE)

# Silhouette width of observation for PAM
sil <- pam.res$silinfo$widths[, 1:3]
# Objects with negative silhouette
neg_sil_index <- which(sil[, 'sil_width'] < 0)
sil[neg_sil_index, , drop = FALSE]
write.table(sil[neg_sil_index, , drop = FALSE] ,file = "negative_cluster_pam", row.names = TRUE, col.names = TRUE)

#library(fpc)
# Compute pairwise-distance matrices
dd <- dist(y.scaled, method ="euclidean")
# Statistics for k-means clustering
km_stats <- cluster.stats(dd,  km.res$cluster)
# (k-means) within clusters sum of squares
km_stats$within.cluster.ss
# (k-means) cluster average silhouette widths
km_stats$clus.avg.silwidths
# Display all statistics
km_stats
# Statistics for pam clustering
pam_stats <- cluster.stats(dd,  pam.res$cluster)
# (pam) within clusters sum of squares
pam_stats$within.cluster.ss
# (pam) cluster average silhouette widths
pam_stats$clus.avg.silwidths
# Statistics for hierarchical clustering
hc_stats <- cluster.stats(dd,  res.hc$cluster)
# (HCLUST) within clusters sum of squares
hc_stats$within.cluster.ss
# (HCLUST) cluster average silhouette widths
hc_stats$clus.avg.silwidths
library("fpc")
# Compute cluster stats
species <- as.numeric(y$Cell_line)
clust_stats <- cluster.stats(d = dist(y.scaled), species, km.res$cluster)
# Corrected Rand index
clust_stats$corrected.rand
# VI
clust_stats$vi
# Agreement between Cellines and K-mean clusters
table(y$Cell_line, km.res$cluster)
write.table(table(y$Cell_line, km.res$cluster),file = "Agreement_betw_Cell_lines_km_cluster", row.names = TRUE, col.names = TRUE)
# Agreement between Cellines and pam clusters
table(y$Cell_line, pam.res$cluster)
write.table(table(y$Cell_line, pam.res$cluster),file = "Agreement_betw_Cell_lines_pam_cluster", row.names = TRUE, col.names = TRUE)
# Agreement between species and HC clusters
table(y$Cell_line, res.hc$cluster)
write.table(table(y$Cell_line, res.hc$cluster),file = "Agreement_betw_Cell_lines_hc_cluster", row.names = TRUE, col.names = TRUE)

write.table(km.res$cluster,file = "intragene_km.res_cluster", row.names = TRUE, col.names = TRUE)
write.table(res.hc$cluster,file = "intragene_res.hc_cluster", row.names = TRUE, col.names = TRUE)
write.table(pam.res$cluster,file = "intragene_pam.res_cluster", row.names = TRUE, col.names = TRUE)
