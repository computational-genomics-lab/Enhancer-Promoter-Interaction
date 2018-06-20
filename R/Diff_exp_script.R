library("edgeR")
y=read.csv("FPKM_All_Celline", sep = "\t", header = T)
z=y[1:4] #subject to change
z=z[,-3] #subject to change
row.names(z) = y$Gene_id
group <- c(rep("C",2),rep("N",1))
d <- DGEList(counts = z, group=group)
d <- calcNormFactors(d)
d <- estimateCommonDisp(d, verbose=T)
d <- estimateTagwiseDisp(d)
de.tgw <- exactTest(d)
summary(decideTestsDGE(de.tgw, p.value=0.05))
topTags(de.tgw)
View(de.tgw$table)
write.table(de.tgw$table, file = "diff_exp", row.names = T, col.names = T)
