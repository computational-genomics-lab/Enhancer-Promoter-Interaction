# Enhancer-Promoter-Interaction
Step 1. Download CAGE DATA from ENCODE
Link : https://genome.ucsc.edu/ENCODE/dataMatrix/encodeDataMatrixHuman.html
For Gm12878
wgEncodeRikenCageGm12878CellPapAlnRep1.bam
wgEncodeRikenCageGm12878CellPapAlnRep2.bam
For K562
wgEncodeRikenCageK562CellPapAlnRep1.bam
wgEncodeRikenCageK562CellPapAlnRep2.bam
For H1hESC
wgEncodeRikenCageH1hescCellPapAlnRep1.bam
wgEncodeRikenCageH1hescCellPapAlnRep2.bam

Run CAGE_analysis.R on data set of an individual cell to get TPM

**Step 2**. Download enhancer-promoter interaction from following links
http://www.cbrc.kaust.edu.sa/dendb/src/enhancer_interactions.tsv.zip

Step 3. Download Gencode file for human from the following links:
https://www.gencodegenes.org/gencodeformat.html
Step 4. Map the promoter region with Gencode file to get the gene information of the promoter
Step 5. Filter the enhancer promoter interaction having gene information
Step 6. Map the CTSS get in the in Step 1 in both enhancer and promoter region. Extract those interactions having TSS mapped with both enhancer and promoter regions
Step 7. Extract those CTSS mapped with enhancers and having a same genomic location. Then separate out CTSS mapped with enhancer and having a different genomic location. Likewise for genic regions.
Step 8. Next separate all the EPI mapped with CTSSs for all cell lines. (Processed file name “TPM_gene_enhancer_interactions_Allchr_diff_dominant_non_redundent” given) 
Step 9. Run Clustering_analysis_&_validation.R script on TPM_gene_enhancer_interactions_Allchr_diff_dominant_non_redundent
Step 10. After successful run of above script one will get corrected cluster file (Corrected_cluster_1, Corrected_cluster_2,  Corrected_cluster_3 given)
Step 11. Download RNA-Seq from ENCODE (Gm12878_1_1, Gm12878_2_1, K562_1_1, K562_2_1, H1hESC_1_1) and extract the FPKM value from all the file to make a cumulative file (FPKM_All_Cellines). Run Diff_exp_script.R on the file (FPKM_All_Cellines) to get Fold_change. Save the result accordingly. Like (K562_Gm12878_edgerAnalysis, K562-H1hesc_edgerAnalysis, Gm12878-H1hesc_edgerAnalysis) given in additional files.
Step 12. For SNP analysis, Download RegulomeDB from the following link:
http://www.regulomedb.org/downloads
Extract the the SNPs information for the corresponding cell lines
In the next step extract the EPI of those 71 genes from the file used in Step 8 (“TPM_gene_enhancer_interactions_Allchr_diff_dominant_non_redundent”)
Extend 2000 bp from the TSS location mapped at enhancer and promoter.
Map individual SNPs from all the cell lines
The Final result has been tabulated in “Contingency_table”.
