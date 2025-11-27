# ==============================================================================
# Load Required Libraries
# ==============================================================================
library(Seurat)
#library(SeuratDisk) # Use when RAM is insufficient, it uses disk space to compensate

# ==============================================================================
# PROCESS CONTROL CONDITION SAMPLES
# ==============================================================================
sample1 <- Read10X(data.dir = "samples/Normal/GSM7290760")
sample1 <- CreateSeuratObject(counts = sample1, project = "N1")
sample1

sample2 <- Read10X(data.dir = "samples/Normal/GSM7290764")
sample2 <- CreateSeuratObject(counts = sample2, project = "N2")
sample2

sample3 <- Read10X(data.dir = "samples/Normal/GSM7290765")
sample3 <- CreateSeuratObject(counts = sample3, project = "N3")
sample3

sample4 <- Read10X(data.dir = "samples/Normal/GSM7290766")
sample4 <- CreateSeuratObject(counts = sample4, project = "N4")
sample4

sample5 <- Read10X(data.dir = "samples/Normal/GSM7290770")
sample5 <- CreateSeuratObject(counts = sample5, project = "N5")
sample5

# ==============================================================================
# Merge Samples and Add Metadata
# ==============================================================================
normal_merged <- merge(sample1, y = c(sample2, sample3, sample4, sample5),
                         add.cell.ids = c("N1", "N2", "N3", "N4", "N5"))

# Check merged samples
head(colnames(normal_merged))
tail(colnames(normal_merged))
unique(sapply(X = strsplit(colnames(normal_merged), split = "_"), FUN = "[", 1))
table(normal_merged$orig.ident)

# Add condition metadata for normal samples
normal_merged <- AddMetaData(object = normal_merged,
                             metadata = "Normal", col.name = "condition")
table(normal_merged@meta.data$condition)

# ==============================================================================
# PROCESS DISEASED CONDITION SAMPLES
# ==============================================================================
sample6 <- Read10X(data.dir = "samples/colon_cancer_liver_metastasis/GSM7290761")
sample6 <- CreateSeuratObject(counts = sample6, project = "LM1")
sample6

sample7 <- Read10X(data.dir = "samples/colon_cancer_liver_metastasis/GSM7290767")
sample7 <- CreateSeuratObject(counts = sample7, project = "LM2")
sample7

sample8 <- Read10X(data.dir = "samples/colon_cancer_liver_metastasis/GSM7290775")
sample8 <- CreateSeuratObject(counts = sample8, project = "LM3")
sample8

sample9 <- Read10X(data.dir = "samples/colon_cancer_liver_metastasis/GSM7290778")
sample9 <- CreateSeuratObject(counts = sample9, project = "LM4")
sample9

sample10 <- Read10X(data.dir = "samples/colon_cancer_liver_metastasis/GSM7290779")
sample10 <- CreateSeuratObject(counts = sample10, project = "LM5")
sample10

# ==============================================================================
# Merge Samples and Add Metadata
# ==============================================================================
LM_merged <- merge(sample6, y = c(sample7, sample8, sample9, sample10),
                   add.cell.ids = c("LM1", "LM2", "LM3", "LM4", "LM5"))

# Add condition metadata for liver metastasis samples
LM_merged <- AddMetaData(object = LM_merged, 
                         metadata = "LM", col.name = "condition")

# ==============================================================================
# Merge Control and Diseased Samples
# ==============================================================================
merged_samples <- merge(normal_merged, LM_merged)

# Clean up intermediate objects
rm(sample1, sample2, sample3, sample4, sample5, sample6,
   sample7, sample8, sample9, sample10, normal_merged, LM_merged)

# ==============================================================================
# Verify Final Merged Object
# ==============================================================================
head(colnames(merged_samples))
tail(colnames(merged_samples))
table(merged_samples$condition)
table(merged_samples$orig.ident)
dim(merged_samples)     # 38224 features (genes) and 36289 cells in the merged samples

# Save Final Merged Object
saveRDS(merged_samples, file = "merged_samples.rds")

# ==============================================================================
# END OF SCRIPT
# ==============================================================================