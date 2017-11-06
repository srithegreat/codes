```r
suppressMessages({
  library('cowplot')
  library('sleuth')
})
# Collect Metadata
METADATA<-read.table("./Metadata.txt",sep='\t', header=T, stringsAsFactors=F, row.names = 1, check.names=FALSE)
# Collect samples failed QC
SAMPLES.EXCLUDE = c((data.table::fread("MayoRNAseq_RNAseq_TCX_QCdetails.txt", data.table=F, header=T))$`Sample Name`)
# Collect samples with Salmon tool output
sample_id <- dir(file.path("/newvolume/RNAseq_Reprocessesing/processed_data/mayo_tcx/salmon_quant/"))
# IMP: NEED TO CHECK HOW MANY SAMPLES ARE IN sample_id and edit 276
sample_id <- sample_id[1:276]
# Filter QC passed samples
sample_id <- setdiff(sample_id, SAMPLES.EXCLUDE)
# Process Metadata
sample_to_condition <- METADATA %>% dplyr::rename(sample = ID)
rownames(sample_to_condition) <- sample_to_condition$sample
# Filter Metadata with samples that have Salmon output
sample_to_condition <- sample_to_condition[sample_id,]
# Collect Sleuth tool input from Salmon output
kal_dirs <- file.path("/newvolume/RNAseq_Reprocessesing/processed_data/mayo_tcx/salmon_quant", sample_id)
s2c <- dplyr::mutate(sample_to_condition, path = kal_dirs)

# Create Sleuth object
so <- sleuth_prep(s2c, extra_bootstrap_summary = TRUE, num_cores = 2)
new_position_theme <- theme(legend.position = c(0.80, 0.90))
png(filename="./figures/all_samples_density.png")
plot_group_density(so, use_filtered = TRUE, units = "est_counts",
  trans = "log", grouping = "Diagnosis", offset = 1) +
  new_position_theme
dev.off()
png(filename="./figures/all_samples_pca.png")
plot_pca(so, color_by = 'Diagnosis', text_labels = TRUE) +
  new_position_theme
dev.off()
```


```r
plot_group_density(so, use_filtered = TRUE, units = "est_counts",
  trans = "log", grouping = "Diagnosis", offset = 1) +
  new_position_theme
```


```r
plot_pca(so, color_by = 'Diagnosis', text_labels = TRUE) +
  new_position_theme
```


```r
# Filter PCA outliers
s2c_a_pc <- dplyr::filter(s2c, !sample %in% c('7095_TCX', '1940_TCX', '1959_TCX', '1950_TCX', '1932_TCX', '11471_TCX', '11327_TCX', '1924_TCX', '11423_TCX'))
# Create Sleuth object
so_a_pc <- sleuth_prep(s2c_a_pc, extra_bootstrap_summary = TRUE, num_cores = 2)
png(filename="./figures/filtered_samples_pca.png")
plot_pca(so_a_pc, color_by = 'Diagnosis', text_labels = TRUE) +
  new_position_theme
dev.off()
```


```r
plot_pca(so_a_pc, color_by = 'Diagnosis', text_labels = TRUE) +
  new_position_theme
```


```r
# Filter First Case and Control
s2c_a_pc_AD_C <- dplyr::filter(s2c_a_pc, Diagnosis %in% c('AD', 'CONTROL'))
# Create Sleuth object
so_AD_C <- sleuth_prep(s2c_a_pc_AD_C, extra_bootstrap_summary = TRUE, num_cores = 2)
# Covariate adjustment
so_AD_C <- sleuth_fit(so_AD_C, ~Source + RIN + RIN2 + Diagnosis + Gender + AgeAtDeath + FLOWCELL + PMI, 'full')
so_AD_C <- sleuth_fit(so_AD_C, ~Source + RIN + RIN2 + Gender + AgeAtDeath + FLOWCELL + PMI, 'reduced')
# Differential expression analysis
so_AD_C <- sleuth_lrt(so_AD_C, 'reduced', 'full')
full_results <- sleuth_results(so_AD_C, 'reduced:full', 'lrt', show_all = FALSE)
sleuth_significant <- dplyr::filter(full_results, qval <= 0.05)
# Prepare result with ID as transcript_id|gene_name
sl_sig <- data.frame("target_id"=matrix(unlist(lapply(strsplit(sleuth_significant$target_id, "|", fixed = TRUE), function(x) paste(x[1], x[6], sep = '|')))), stringsAsFactors=FALSE)
sl_sig_rem <- (sleuth_significant[,2:13])
sl_sig <- cbind(sl_sig, sl_sig_rem)
write.table(sl_sig, file = "./sleuth/AD_C_signi_v1.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(head(sl_sig, 20), file = "./sleuth/AD_C_signi20.txt", quote = FALSE, sep = "\t", row.names = FALSE)

# Filter Second Case and Control
s2c_a_pc_PA_C <- dplyr::filter(s2c_a_pc, Diagnosis %in% c('PA', 'CONTROL'))
# Create Sleuth object
so_PA_C <- sleuth_prep(s2c_a_pc_PA_C, extra_bootstrap_summary = TRUE, num_cores = 2)
# Covariate adjustment
so_PA_C <- sleuth_fit(so_PA_C, ~Source + RIN + RIN2 + Diagnosis + Gender + AgeAtDeath + FLOWCELL + PMI, 'full')
so_PA_C <- sleuth_fit(so_PA_C, ~Source + RIN + RIN2 + Gender + AgeAtDeath + FLOWCELL + PMI, 'reduced')
# Differential expression analysis
so_PA_C <- sleuth_lrt(so_PA_C, 'reduced', 'full')
full_resultsPA <- sleuth_results(so_PA_C, 'reduced:full', 'lrt', show_all = FALSE)
sleuth_significantPA <- dplyr::filter(full_resultsPA, qval <= 0.05)
# Prepare result with ID as transcript_id|gene_name
sl_sigPA <- data.frame("target_id"=matrix(unlist(lapply(strsplit(sleuth_significantPA$target_id, "|", fixed = TRUE), function(x) x[1]))), stringsAsFactors=FALSE)
sl_sig_remPA <- (sleuth_significantPA[,2:13])
sl_sigPA <- cbind(sl_sigPA, sl_sig_remPA)
write.table(sl_sigPA, file = "./sleuth/PA_C_signi_v1.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(head(sl_sigPA, 20), file = "./sleuth/PA_C_signi20.txt", quote = FALSE, sep = "\t", row.names = FALSE)

# Filter Third Case and Control
s2c_a_pc_PSP_C <- dplyr::filter(s2c_a_pc, Diagnosis %in% c('PSP', 'CONTROL'))
# Create Sleuth object
so_PSP_C <- sleuth_prep(s2c_a_pc_PSP_C, extra_bootstrap_summary = TRUE, num_cores = 2)
# Covariate adjustment
so_PSP_C <- sleuth_fit(so_PSP_C, ~Source + RIN + RIN2 + Diagnosis + Gender + AgeAtDeath + FLOWCELL + PMI, 'full')
so_PSP_C <- sleuth_fit(so_PSP_C, ~Source + RIN + RIN2 + Gender + AgeAtDeath + FLOWCELL + PMI, 'reduced')
# Differential expression analysis
so_PSP_C <- sleuth_lrt(so_PSP_C, 'reduced', 'full')
full_resultsPSP <- sleuth_results(so_PSP_C, 'reduced:full', 'lrt', show_all = FALSE)
sleuth_significantPSP <- dplyr::filter(full_resultsPSP, qval <= 0.05)
# Prepare result with ID as transcript_id|gene_name
sl_sigPSP <- data.frame("target_id"=matrix(unlist(lapply(strsplit(sleuth_significantPSP$target_id, "|", fixed = TRUE), function(x) x[1]))), stringsAsFactors=FALSE)
sl_sig_remPSP <- (sleuth_significantPSP[,2:13])
sl_sigPSP <- cbind(sl_sigPSP, sl_sig_remPSP)
write.table(sl_sigPSP, file = "./sleuth/PSP_C_signi_v1.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(head(sl_sigPSP, 20), file = "./sleuth/PSP_C_signi20.txt", quote = FALSE, sep = "\t", row.names = FALSE)

# Get count of isoforms differentially expressed in each comparison
tot1 <- length(sl_sig$target_id)
tot2 <- length(sl_sigPA$target_id)
tot3 <- length(sl_sigPSP$target_id)
# Get count of isoforms intersecting in each comparisons
inter12 <- length(intersect(sl_sig$target_id, sl_sigPA$target_id))
inter13 <- length(intersect(sl_sig$target_id, sl_sigPSP$target_id))
inter23 <- length(intersect(sl_sigPA$target_id, sl_sigPSP$target_id))
inter123 <- length(intersect(intersect(sl_sig$target_id, sl_sigPA$target_id), sl_sigPSP$target_id))

# Plot Venn Diagram
library(VennDiagram)
png(filename="./figures/comp_venn_v1.png")
draw.triple.venn(area1 = tot1, area2 = tot2, area3 = tot3, n12 = inter12, n23 = inter23, n13 = inter13, 
    n123 = inter123, category = c("AD vs CON", "PA vs CON", "PSP vs CON"), lty = "blank", 
    fill = c("skyblue", "pink1", "mediumorchid"))
dev.off()
```

```r
draw.triple.venn(area1 = tot1, area2 = tot2, area3 = tot3, n12 = inter12, n23 = inter23, n13 = inter13, 
    n123 = inter123, category = c("AD vs CON", "PA vs CON", "PSP vs CON"), lty = "blank", 
    fill = c("skyblue", "pink1", "mediumorchid"))
```
