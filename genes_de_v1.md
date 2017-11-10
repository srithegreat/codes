---
title: "Differentially expression analysis of genes in mouse genome"
author: "Philge Philip"
date: "11/7/2017"
---

## Collecting genes reads count

```r
library(edgeR)
library(ggplot2)
library(cqn)
library(RColorBrewer)
library(pheatmap)
library(biomaRt)
library(org.Mm.eg.db)
library(clusterProfiler)
library(pathview)

# Collect expression data
A_1 <- read.table("./1.genes.results", sep='\t', header=T, stringsAsFactors=F, check.names=FALSE)
A_2 <- read.table("./2.genes.results", sep='\t', header=T, stringsAsFactors=F, check.names=FALSE)
B_1 <- read.table("./13.genes.results", sep='\t', header=T, stringsAsFactors=F, check.names=FALSE)
B_2 <- read.table("./14.genes.results", sep='\t', header=T, stringsAsFactors=F, check.names=FALSE)
# Extract expected reads count
exp <- merge(A_1[,c("gene_id", "expected_count")], A_2[,c("gene_id", "expected_count")], by="gene_id")
colnames(exp) <- c("gene_id", "A_1","A_2")
exp <- merge(exp, B_1[,c("gene_id", "expected_count")], by="gene_id")
exp <- merge(exp, B_2[,c("gene_id", "expected_count")], by="gene_id")
colnames(exp) <- c("gene_id", "A_1","A_2", "B_1", "B_2")
rownames(exp) <- exp[,1]
exp <- exp[,-1]
exp.full <- exp
# genes filter
keep <- rowSums(cpm(exp)>1) >= 2
exp <- exp[keep,]
```

#### 12463 genes filtered which has a cpm of greater than 1 for at least two samples.

## Genes Density

```r
nsamples <- ncol(exp)
col <- brewer.pal(nsamples, "Paired")
lcpm <- cpm(exp.full, log=TRUE)
lcpm1 <- cpm(exp, log=TRUE)
png(filename="./figures/data_filtering.png")
{par(mfrow=c(1,2))
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.6), las=2, 
     main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
 den <- density(lcpm[,i])
 lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", colnames(exp.full), text.col=col, bty="n")
plot(density(lcpm1[,1]), col=col[1], lwd=2, ylim=c(0,0.2), las=2, 
     main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
   den <- density(lcpm1[,i])
   lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", colnames(exp), text.col=col, bty="n")}
dev.off()
```
![alt text](./figures/data_filtering.png "Logo Title Text 1")
### Figure: The density of log-CPM values for raw pre-filtered data (A) and post-filtered data (B)

## PCA to explore outliers

```r
PC <- prcomp(voom(t(exp))$E, scale.=T, center = T)
samples <- rownames(PC$x)
png(filename="./figures/PCA_outlier_analysis.png")
ggplot(data.frame(PC$x),aes(PC1,PC2,color=samples)) + 
  geom_point(size = 5) +
  ggtitle("PCA after voom normalization") + 
  theme_bw() + 
  theme(plot.title = element_text(face="bold", color="black", size=14, hjust=0.5))
dev.off()
```
![alt text](./figures/PCA_outlier_analysis_v1.png "Logo Title Text 2")
#### Samples in condition B are dissimilar though belonging to same group; decided to proceed with them coz of no other replicates

## Calculating genes GC content and length

### Extracting gene coordinates

```perl
open(FILE, "<", "gencode.vM14.chr_patch_hapl_scaff.annotation.gtf");
open(OUT, ">", "gencode.vM14.chr_patch_hapl_scaff.annotation.bed");
while($line = <FILE>){
	if($line =~ /\tgene\t/){
		@cols = split('\t', $line);
		if($cols[8] =~ /gene_id\s\"(.*)\";\sgene_type.*/){
			print OUT "$cols[0]\t$cols[3]\t$cols[4]\t$1\n";
		}
	}
}
close(FILE);
close(OUT);
exit;
```
### Using bedtools and emboss tool to calculate percentage GC content and length

```bash
bedtools getfasta -fi GRCm38.p5.genome.fa -bed gencode.vM14.chr_patch_hapl_scaff.annotation.bed\
-name -fo GRCm38_genes.fa
infoseq -auto -only -name -length -pgc GRCm38_genes.fa >GRCm38_length_gc.tsv
perl -pe 's/ +/\t/g' GRCm38_length_gc.tsv >GRCm38_length_gc_v1.tsv
```

## Conditional Quantile Normalisation (CQN)

### Initial normalisation usign cqn which remove the influence of GC content on counts, smooth the effect of gene length.

```r
# Collect gene length and GC content
GENE.LEN.GC <- read.table("./GRCm38_length_gc_v1.tsv", sep='\t', header=T, row.names = 1,\
stringsAsFactors=F, check.names=FALSE)
# filtering genes in expression data
GENE.GC <- data.frame(pgc_content = GENE.LEN.GC[rownames(exp), '%GC'], row.names = rownames(exp))
GENE.LEN <- data.frame(gene.length = GENE.LEN.GC[rownames(exp), 'Length'], row.names = rownames(exp))
CQN.GENE_EXPRESSION <- cqn(exp, 
                          x = GENE.GC$pgc_content,
                          lengths = GENE.LEN$gene.length,
                          lengthMethod = "smooth", 
                          verbose = FALSE)
CQN.GENE_EXPRESSION$E <- CQN.GENE_EXPRESSION$y + CQN.GENE_EXPRESSION$offset
png(filename="./figures/data_normalization.png")
{par(mfrow=c(1,2))
lcpm <- cpm(exp, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="A. Unnormalised",ylab="Log-cpm")
boxplot(CQN.GENE_EXPRESSION$E, las=2, col=col, main="")
title(main="B. CQN Normalised",ylab="Log-cpm")}
dev.off()
```
![alt text](./figures/data_normalization.png "Logo Title Text 3")
### Figure: Boxplots of log-CPM values showing expression distributions for unnormalised data (A) and normalised data (B)

## Differential expression (DE) analysis

### Dispersion control with voom

```r
# Construct design
group <- as.factor(c("A", "A", "B", "B"))
design <- model.matrix(~0+group)
colnames(design) <- gsub("group", "", colnames(design))
png(filename="./figures/mean-variance.png")
{par(mfrow=c(1,2))
# Estimate voom weights for dispersion control
VOOM.WEIGHTS <- voom(exp, design, plot=T)
# Fit linear model using weights and design
VOOM.WEIGHTS$E <- CQN.GENE_EXPRESSION$E
FIT <- lmFit(VOOM.WEIGHTS)
# Fit contrast
contrast <- makeContrasts(AvsB=A-B, levels = colnames(design))
FIT.CONTR <- contrasts.fit(FIT, contrasts=contrast)
FIT.CONTR <- eBayes(FIT.CONTR)
plotSA(FIT.CONTR, main="Final model: Mean−variance")}
dev.off()
# Get differential expression
DE <- topTable(FIT.CONTR, number = Inf, adjust.method = "BH", sort.by = "p")
```
![alt text](./figures/mean-variance.png "Logo Title Text 4")
#### The mean-variance relationship of log-CPM values for this dataset is shown in the left-hand panel. Typically, the voom-plot shows a decreasing trend between the means and variances resulting from a combination of technical variation in the sequencing experiment and biological variation amongst the replicate samples from different cell populations. Experiments with high biological variation usually result in flatter trends, where variance values plateau at high expression values. Experiments with low biological variation tend to result in sharp decreasing trends. The model’s residual variances are plotted against average expression values in the right-hand panel.

```r
# Get top 500 genes
DE_500 <- head(DE, 500)
DE.topgenes <- rownames(DE_500)
# Remove gene version number
DE.topgenes_1 <- sapply(strsplit(DE.topgenes, "[.]"), `[`, 1)
# Collect gene symbols from Ensembl
ensembl = useEnsembl(biomart="ensembl", version=89, dataset="mmusculus_gene_ensembl")
g_symb <- getBM(filters= "ensembl_gene_id", values=DE.topgenes_1, attributes=  c("ensembl_gene_id",\
"external_gene_name"),
                 mart= ensembl)
rownames(g_symb) <- g_symb$ensembl_gene_id
# Get genes in order as DE_500
g_symb <- g_symb[DE.topgenes_1,]
# Add gene name to DE_500
DE_500 <- cbind(ensembl_gene_id = rownames(DE_500), external_gene_name = g_symb$external_gene_name, DE_500)
# print top 10 genes
head(DE_500, 10)
write.table(DE_500, file = "./top_500_deg.txt", quote = FALSE, sep = "\t", row.names = FALSE)
```
ensembl_gene_id|external_gene_name|logFC|AveExpr|t|P.Value|adj.P.Val|B
---|---|---|---|---|---|---|---
ENSMUSG00000045573.9|Penk|2.542816496|10.30388144|195.5607573|7.36E-16|4.63E-12|27.19223002
ENSMUSG00000062070.12|Pgk1|-1.65249878|11.23646437|-195.3497142|7.42E-16|4.63E-12|27.3819314
ENSMUSG00000027720.7|Il2|4.064011576|8.47169247|180.7514066|1.37E-15|5.70E-12|25.81473792
ENSMUSG00000037071.2|Scd1|-3.041164966|8.43935736|-142.5512761|8.97E-15|2.79E-11|24.53029122
ENSMUSG00000049775.16|Tmsb4x|-1.003137874|12.564248|-130.220929|1.83E-14|4.06E-11|24.06088064
ENSMUSG00000025161.16|Slc16a3|-3.003606063|7.751178809|-127.4049851|2.18E-14|4.06E-11|23.48073218
ENSMUSG00000064363.1|mt-Nd4|1.25126026|11.60642329|124.8993995|2.55E-14|4.06E-11|23.80767785
ENSMUSG00000028965.13|Tnfrsf9|-2.325830783|8.77747308|-124.5765966|2.60E-14|4.06E-11|23.78605475
ENSMUSG00000026249.10|Serpine2|4.423638752|6.970224139|122.0486693|3.06E-14|4.24E-11|22.63097952
ENSMUSG00000050912.15|Tmem123|-1.903255652|8.16040951|-111.8056117|6.12E-14|7.49E-11|22.93191957
### Table: Top 10 differentially expressed genes (deg); top 500 deg can be found in ./codes/data/top_500_deg.txt

### Heat map showing the normalized expression of top 500 DE genes

```r
Samples = c("A", "A", "B", "B")
pheatmap(VOOM.WEIGHTS$E[DE.topgenes,],  color = colorRampPalette((brewer.pal(n = 9, name = "OrRd")))(9),\
border_color = NA, cellwidth = 20, cellheight = 3, show_colnames = T, annotation_col = data.frame(Samples,\
row.names = c("A_1", "A_2", "B_1", "B_2")), annotation_names_col = FALSE, labels_row = \
g_symb$external_gene_name, labels_col = group, filename = "./figures/heatmap.png", fontsize_row = 3,\
clustering_distance_cols = "euclidean", clustering_distance_rows = "euclidean")
```
![alt text](./figures/heatmap.png "Logo Title Text 5")

## Gene ontology of the differentially expressed genes

```r
# Collect fold change of DE genes
fold <- setNames(DE_500[,"logFC"], DE.topgenes_1)
fold <- sort(fold, decreasing = T)
# Prepare Universal set
DE.genes <- rownames(DE)
# Remove gene version number
DE.genes_1 <- sapply(strsplit(DE.genes, "[.]"), `[`, 1)
# GO over-representation test
ego <- enrichGO(gene          = DE.topgenes_1,
                universe      = DE.genes_1,
                keytype       = 'ENSEMBL',
                OrgDb         = org.Mm.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                minGSSize     = 10,
                maxGSSize     = 500,
                readable      = TRUE)
# simplify enriched GO result
ego2 <- simplify(ego, cutoff=0.7, by="p.adjust", select_fun=min)
head(ego2, 5)
write.csv(ego2, file = "./gene_ontology.txt", quote = FALSE, row.names = FALSE)
```
ID|Description|GeneRatio|BgRatio|pvalue|p.adjust|qvalue|geneID|Count
---|---|---|---|---|---|---|---|---
GO:0046649|lymphocyte activation|64/488|461/11099|9.69798632017757e-17|3.79094285255741e-13|2.55618502586575e-13|Il2/Cd2/Cd6/Ctla4/Hspd1/Irf4/Tnfrsf13b/Il2ra/Cdkn1a/Traf6/Cd5/Cd47/Prnp/Themis/Sdc4/Ap3b1/Xbp1/Cd44/Lgals1/Rab27a/Bcl2/Tnfsf11/Nfatc1/Stat3/Ptprc/Rbpj/Ebi3/Bax/Ccr6/Prex1/Ptk2b/Egr1/Cd81/Icos/Cd28/Treml2/Tnfsf4/Lcp1/Thy1/Fam49b/Itm2a/Foxp1/Spn/Il18r1/Zeb1/Sema4a/Nfkbid/Pik3r1/Rora/Foxp3/Itgb2/Il2rg/Stat5b/Irs2/Cd83/Il2rb/Slfn1/Pou2f2/Lrrc32/H2-T23/Hdac7/Cd4/Btla/Vsir|64
GO:0042110|T cell activation|52/488|333/11099|6.85166095016204e-16|1.33915713270917e-12|9.02976790484513e-13|Il2/Cd2/Cd6/Ctla4/Hspd1/Irf4/Il2ra/Traf6/Cd5/Cd47/Prnp/Themis/Sdc4/Ap3b1/Xbp1/Cd44/Lgals1/Rab27a/Bcl2/Tnfsf11/Stat3/Ptprc/Ebi3/Bax/Ccr6/Prex1/Egr1/Icos/Cd28/Treml2/Tnfsf4/Lcp1/Thy1/Fam49b/Foxp1/Spn/Il18r1/Zeb1/Sema4a/Nfkbid/Rora/Foxp3/Itgb2/Il2rg/Stat5b/Cd83/Slfn1/Lrrc32/H2-T23/Cd4/Btla/Vsir|52
GO:0001816|cytokine production|59/488|421/11099|1.12278993515471e-15|1.46299528550659e-12|9.86479297413123e-13|Il2/Tnfrsf9/Rnf19b/Cd2/Furin/Cd6/Ltb/Hspd1/Lta/Irf4/Ccl3/Traf6/Hif1a/Ndrg2/Trim30a/Prnp/Ncl/Xbp1/Il9/Ccl4/Il12rb2/Rbpj/Eif2ak2/Arrb1/Egr1/Hilpda/Irf8/Cd28/Tnfsf4/Il6ra/Pdcd4/Irf7/Zc3hav1/Fam49b/Isg15/Foxp1/Spn/Fn1/Tspo/Il17f/Il18r1/Tnf/Cadm1/Pik3r1/Rora/Anxa4/Foxp3/Slc37a4/Stat5b/S1pr3/Mapk3/Runx3/Cd83/Srgn/Lrrc32/H2-T23/Hdac7/Cd4/Vsir|59
GO:0001817|regulation of cytokine production|54/488|377/11099|8.00268031391472e-15|6.53226795508263e-12|4.40462602079598e-12|Il2/Tnfrsf9/Cd2/Furin/Cd6/Ltb/Hspd1/Lta/Irf4/Ccl3/Traf6/Hif1a/Ndrg2/Trim30a/Prnp/Ncl/Xbp1/Il9/Ccl4/Il12rb2/Eif2ak2/Arrb1/Egr1/Hilpda/Irf8/Cd28/Tnfsf4/Il6ra/Pdcd4/Irf7/Zc3hav1/Fam49b/Isg15/Spn/Fn1/Tspo/Il17f/Il18r1/Tnf/Cadm1/Pik3r1/Rora/Anxa4/Foxp3/Slc37a4/Stat5b/S1pr3/Mapk3/Cd83/Srgn/Lrrc32/H2-T23/Hdac7/Vsir|54
GO:0098602|single organism cell adhesion|56/488|402/11099|8.3554207662863e-15|6.53226795508263e-12|4.40462602079598e-12|Il2/Serpine2/Cd2/Cd6/Ctla4/Hspd1/Il2ra/Traf6/Cd5/Cd47/Prnp/Sdc4/Ap3b1/Xbp1/Cd44/Lgals1/Perp/Pkp4/Bcl2/Tnfsf11/Selplg/Ptprc/Map2k1/Cyr61/Ets1/Prex1/Nlgn2/Icos/Cd28/Tnfsf4/Myl12b/Vegfa/Plek/Thy1/Flna/Fam49b/Spn/Fn1/Tnf/Calr/Zeb1/Nfkbid/Pik3r1/Foxp3/Itgb2/Il2rg/Stat5b/Cd83/St6gal1/Slfn1/Lrrc32/Pvr/H2-T23/Cd4/Btla/Vsir|56
### Table: Top 5 over represented gene ontology (biological process); all over represented gene ontology can be found in ./codes/data/gene_ontology.txt

```r
cnetplot(ego2, categorySize="pvalue", foldChange=fold, fixed = F)
```
![alt text](./figures/gene_ontology.png "Logo Title Text 6")
### Figure: Association of differentially expressed genes with gene ontology

## Pathway analysis on the differentially expressed genes

```r
# get uniprot ids
up_1 = bitr(DE.topgenes_1, fromType="ENSEMBL", toType=c("UNIPROT"), OrgDb="org.Mm.eg.db")
up = bitr(DE.genes_1, fromType="ENSEMBL", toType=c("UNIPROT"), OrgDb="org.Mm.eg.db")
# KEGG over-representation test
kk <- enrichKEGG(gene          =  up_1$UNIPROT,
                 universe      = up$UNIPROT,
                 organism      = 'mmu',
                 keyType       = "uniprot",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05,
                 minGSSize     = 10,
                 maxGSSize     = 500)
head(kk, 5)
write.csv(head(kk, 5), file = "./pathway_analysis.txt", quote = FALSE, row.names = FALSE)
```
ID|Description|GeneRatio|BgRatio|pvalue|p.adjust|qvalue|geneID|Count
---|---|---|---|---|---|---|---|---
mmu04066|HIF-1 signaling pathway|37/452|122/6476|4.09063034112917e-15|1.04311073698794e-12|6.84642341304776e-13|P09411/D2KHZ9/P16858/Q5FWB7/P05064/P39689/Q564P6/P17809/Q61221/Q3UCW2/P12382/Q8CD98/Q3TTB4/Q6GQU1/Q91YE3/P10417/P42227/Q3ULI4/A0A0R4J0H9/Q91UZ4/Q8BFP9/P31938/Q3TMJ8/P46414/Q60876/O08528/Q3URV7/P22272/Q00731/P17182/Q5FW97/Q3TP23/Q8C7P2/P26450/Q564E2/P06151/Q63844|37
mmu04060|Cytokine-cytokine receptor interaction|49/452|224/6476|1.9678269026006e-13|2.50897930081576e-11|1.64676040796576e-11|P04351/P20334/Q3U3R1/A0A0U5JAA2/P41155/P09225/Q542S2/A5D8Y6/Q9ET35/P01590/Q544I2/P10855/Q5QNW0/P09056/P15247/P14097/Q5QNV9/D3Z6H5/P97378/P56484/Q3ZB17/O35235/P01586/Q5SX77/P51680/O54689/Q542B6/Q3UDV8/Q3UDL4/Q61727/B6DXE3/P43488/Q3URV7/P22272/Q00731/O55237/Q05A52/Q7TNI7/Q61098/P06804/Q3U593/O35714/Q540M6/Q3UPA9/P34902/Q9QZM4/Q04998/Q3UY39/P16297|49
mmu05230|Central carbon metabolism in cancer|28/452|86/6476|1.43094835867036e-12|1.21630610486981e-10|7.98318557995044e-11|P57787/Q3TMA0/Q3UDP9/P17809/Q61221/Q3UCW2/P12382/Q8CD98/Q3TTB4/Q6GQU1/O08586/Q3UUT8/Q9Z127/Q8C605/Q9WUA3/Q3U7Z6/Q9DBJ1/Q8BFP9/P31938/Q3TMJ8/O08528/Q3TP23/Q8C7P2/P26450/Q564E2/P06151/Q63844/Q571F8|28
mmu00010|Glycolysis / Gluconeogenesis|25/452|70/6476|2.21778087592155e-12|1.41383530839999e-10|9.27966208609282e-11|P09411/D2KHZ9/P16858/Q5FWB7/P05064/P06745/P17751/P12382/Q8CD98/Q3TTB4/Q6GQU1/Q3U6X6/Q9D0F9/Q8C605/Q9WUA3/P05063/Q3U7Z6/Q9DBJ1/A0A0R4J0G0/Q8BH04/O08528/P17182/Q5FW97/Q564E2/P06151|25
mmu04514|Cell adhesion molecules (CAMs)|30/452|119/6476|3.04684214895518e-10|1.55388949596714e-08|1.01989031933447e-08|P08920/Q549Q4/Q61003/P09793/Q6GTR6/P14429/O35988/Q3U5S6/Q02242/Q544F3/Q3TA56/Q8C6Q7/S4R1S4/Q69ZK9/Q6PHN2/Q3V3X2/Q5SUZ7/Q9WVS0/P31041/Q8CDB3/P15702/Q544C5/Q8R5M8/P11835/Q542I8/Q8K094/P06339/Q3V014/P06332/Q3TSV7|30
### Table: Top 5 enriched pathways

```r
mmu04066 <- pathview(gene.data  = fold,
                     pathway.id = "mmu04066",
                     species    = "mmu",
                     limit      = list(gene=max(abs(fold)), cpd=1),
                     gene.idtype= "ENSEMBL")
```
![alt text](./figures/mmu04066.pathview.png "Logo Title Text 7")
### Figure: Differentially expressed genes in enriched HIF-1 signalling pathway 
