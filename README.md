![cover](data/00-coverpage.png)


## 📑 Table of contents

- [📝 Overview](#-overview)
- [Installation](#installation)
- [Input Files Format](#input-files-format)
- [01. Merging Samples](#01-merging-samples)
- [02. Single Cell RNA Sequencing Analysis](#02-single-cell-rna-sequencing-analysis)
- [03. Multiple Datasets Integration](#03-multiple-datasets-integration)
- [04. 3D UMAP and t-SNE Plots](#04-3d-umap-and-t-sne-plots)
- [05. Volcano Plot](#05-volcano-plot)
- [06. Heatmap](#06-heatmap)
- [07. Boxplot](#07-boxplot)
- [📚 Citation](#-citation)
- [🤝 Acknowledgements](#-acknowledgements)


## 📝 Overview
<p align="justify"> <b>scRNA-seq-navigator</b> contains a comprehensive single-cell RNA sequencing (scRNA-seq) analysis pipeline designed for processed data following alignment. The pipeline leverages the robust <a href="https://doi.org/10.1016/j.cell.2019.05.031">Seurat (Stuart et al., 2019)</a> toolkit for essential analytical steps, including data integration, quality control, and batch-effect correction, and incorporates the <a href="https://doi.org/10.1038/s41467-022-28803-w">ScType database (Ianevski et al., 2022)</a> for refined cell-type annotation. It further facilitates differential gene expression analysis and multi-dimensional visualization to extract meaningful biological insights from complex single-cell datasets. </p>

- The analyses provide end-to-end solutions for:
    - **Data Integration & Processing:** Merging multiple datasets and experimental conditions
    - **Batch Effect Correction:** Addressing technical variations across different platforms and studies
    - **Cell Type Identification:** Automated annotation using reference gene signatures
    - **Advanced Visualization:** 2D/3D embeddings and interactive plots
    - **Differential Expression:** Statistical analysis of gene expression changes
    - **Publication-ready Outputs:** High-quality figures and comprehensive results
- <p align="justify"> Each script represents a modular component of the analysis workflow, allowing researchers to process single-cell data from raw counts to biological interpretation </p>

| Script | Description |
|--------|-------------|
| `01-samples_merger.R` | Merges multiple scRNA-seq samples from different conditions into a unified dataset for analysis |
| `02-scRNA-seqr.R` | Performs comprehensive scRNA-seq analysis including QC, clustering, cell type annotation, and differential expression analysis |
| `03-multiple-datasets-integrator.R` | Integrates multiple scRNA-seq datasets using Harmony to remove batch effects |
| `04-umap-tsne-3d-plotr.R` | Creates interactive 3D visualizations of UMAP and t-SNE embeddings for enhanced data exploration |
| `05-volcano_plotr.R` | Generates publication-ready volcano plots to visualize differential gene expression results |
| `06-heatmapr.R` | Produces clustered heatmaps showing gene expression patterns across cell types and conditions |
| `07-boxplotr.R` | Creates comparative box plots to analyze gene expression differences across experimental conditions and/or cell types |

## Installation

- Navigate to the directory containing `install_packages.R` file
- Open terminal and run:

```r
Rscript install_packages.R
```

- Alternatively, if the above does not work, open `install_packages.R` file in RStudio
- Click the "Source" button in the top-right of the script editor OR in R console run:

```
# Navigate to the script directory first, then:
source("install_packages.R")
```

## Input Files Format

- The analysis workflow is designed for processed scRNA-seq data following alignment.
- Input should be in 10X Genomics format, consisting of the three standard files for each sample directory:
    - `barcodes.tsv.gz`
    - `features.tsv.gz`
    - `matrix.mtx.gz`

## 01. Merging Samples

<p align="justify"> This script performs the initial data integration step for single-cell RNA sequencing (scRNA-seq) analysis by merging multiple samples from two experimental conditions (e.g Normal and Parkinson's Disease) into a unified Seurat object. </p>

### Purpose

- Load individual scRNA-seq samples from both control and diseased conditions
- Merge samples within each condition group
- Combine both condition groups into a final merged dataset
- Add appropriate metadata for sample tracking and condition identification
- Export the consolidated data for downstream analysis

### Workflow

### A. Reading and Creating Seurat Object

```r
sample1 <- Read10X(data.dir = "samples/Normal/GSM7290760")
sample1 <- CreateSeuratObject(counts = sample1, project = "N1")
```

- `Read10X()` Reads the 10X Genomics format data from the specified directory containing three essential files: 
    - `barcodes.tsv.gz` (cell identifiers)
    - `features.tsv.gz` (gene information)
    - `matrix.mtx.gz`   (count matrix)
- `CreateSeuratObject()` Converts the raw count data into a Seurat object, which is the standard data structure

### B. Merging Multiple Samples

- Merging all control samples together:

```r
normal_merged <- merge(sample1, y = c(sample2, sample3, sample4, sample5),
                       add.cell.ids = c("N1", "N2", "N3", "N4", "N5"))
```

- Merging all diseased samples together:

```r
LM_merged <- merge(sample6, y = c(sample7, sample8, sample9, sample10),
                   add.cell.ids = c("LM1", "LM2", "LM3", "LM4", "LM5"))
```

### C. Adding Condition Metadata

- Control samples:

```r
normal_merged <- AddMetaData(object = normal_merged,
                             metadata = "Normal", col.name = "condition")
```

- Diseased samples:

```r
LM_merged <- AddMetaData(object = LM_merged, 
                         metadata = "LM", col.name = "condition")
```

### D. Final Merging

- Merging control and diseased samples together:

```r
merged_samples <- merge(normal_merged, LM_merged)
```

## 02. Single Cell RNA Sequencing Analysis

<p align="justify"> This script performs comprehensive single-cell RNA sequencing analysis from quality control through cell type annotation and differential expression analysis. </p>

### Workflow

### A. Data Loading & Quality Control

- Loads merged dataset from previous script (samples merger)

```r
rawdata <- merged_samples
dim(rawdata)  # 38224 features and 36289 cells
rm(merged_samples)
```

- Calculates mitochondrial gene percentage

```r
rawdata[["percent.mt"]] <- PercentageFeatureSet(rawdata, pattern = "^(?i)MT-")
```

- Visualizes QC metrics using violin plots

```r
Plot1 <- VlnPlot(rawdata, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                 group.by = "condition", ncol = 3)
Plot1
```

- Applies quality filters based on:
    - Gene counts
    - RNA counts
    - Mitochondrial content

```r
filtered_data <- subset(rawdata, subset = nFeature_RNA > 200
                        & nFeature_RNA < 6000 & nCount_RNA > 2000 & percent.mt < 10)
dim(filtered_data) # 38224 features and 16195 cells in filtered data
```

<div style="display:flex; justify-content:space-between;">
  <img src="data/01-rawdata.png" width="49%">
  <img src="data/02-filtered_data.png" width="49%">
</div>

**Note:**
- **Before saving plots, adjust the plot viewing window size in RStudio to match the intended resolution and aspect ratio.**
- **After resizing the plot window, save the plot using `ggsave()` with the following format:**

```r
ggsave(file = "raw_data.png", Plot1, dpi = 600, bg = "white")
```

### B. Data Normalization & Feature Selection

- Normalizes data using log normalization

```r
normalized_data <- NormalizeData(filtered_data, normalization.method = "LogNormalize", scale.factor = 10000)
```

- Identifies highly variable genes

```r
normalized_data <- FindVariableFeatures(normalized_data, selection.method = "vst", nfeatures = 2000)
```

- Visualizes top variable features

```r
top10 <- head(VariableFeatures(normalized_data), 10)
top10
```

<img src="data/03-hvfs.png" width="400px" />

### C. Dimensionality Reduction & Clustering

- Performs scaling of highly variable features and performs PCA

```r
scaled_data <- ScaleData(normalized_data)
scaled_data <- RunPCA(scaled_data, features = VariableFeatures(object = scaled_data))
pca_plot <- DimPlot(scaled_data, reduction = "pca")
pca_plot
```

**Note:**
- **If you want to perform scaling of all genes, then run the following command instead of the previous one:**

```r
all.genes <- rownames(normalized_data)
scaled_data <- ScaleData(normalized_data, features = all.genes)
scaled_data <- RunPCA(scaled_data, features = all.genes)
pca_plot <- DimPlot(scaled_data, reduction = "pca")
pca_plot
```

- Determines optimal dimensions using elbow plot

```r
elbow <- ElbowPlot(scaled_data, ndim = 50)
elbow
```

- Clusters cells using graph-based clustering

```r
data <- FindNeighbors(scaled_data, dims = 1:20)
data <- FindClusters(data, resolution = 0.6, algorithm = 4)
levels(data)
```

- Performs non-linear dimensionality reduction (UMAP & t-SNE)

```r
data <- RunUMAP(data, dims = 1:20)
data <- RunTSNE(data, dims = 1:20)
```

- Visualizes clusters with and without conditions

```r
umap <- DimPlot(data, reduction = "umap", label = TRUE)
umap

tSNE <- DimPlot(data, reduction = "tsne", label = TRUE)
tSNE

umap_condition <- DimPlot(data, reduction = "umap", label = TRUE,
                          split.by = "condition")
umap_condition

tSNE_condition <- DimPlot(data, reduction = "tsne", label = TRUE,
                          split.by = "condition")
tSNE_condition

```

---

### PCA and Elbow Plot

<div style="display:flex; justify-content:space-between;">
  <img src="data/04-pca.png" width="49%">
  <img src="data/05-elbow.png" width="49%">
</div>

---

### UMAP and t-SNE - Cell Clusters

<div style="display:flex; justify-content:space-between;">
  <img src="data/06-umap.png" width="49%">
  <img src="data/08-tsne.png" width="49%">
</div>

---

### UMAP - Clusters Grouped By Condition

![umap-condition](data/07-umap_condition.png)

---

### t-SNE - Clusters Grouped By Condition

![tsne-condition](data/09-tsne_condition.png)

### D. Cell Type Annotation

- Uses ScType database for automated cell type identification
- Annotates clusters based on tissue markers
- Assigns cell types and handles low-confidence annotations
- Updates metadata with cell type classifications

```r
# load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

# DB file
db_ <- "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
tissue <- "Immune system" # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus 

# prepare gene sets
gs_list <- gene_sets_prepare(db_, tissue)

# ==============================================================================
# Run scType Annotation
# ==============================================================================
# check Seurat object version (scRNA-seq matrix extracted differently in Seurat v4/v5)
seurat_package_v5 <- isFALSE('counts' %in% names(attributes(data[["RNA"]])));
print(sprintf("Seurat object %s is used", ifelse(seurat_package_v5, "v5", "v4")))

# extract scaled scRNA-seq matrix
scRNAseqData_scaled <- if (seurat_package_v5) as.matrix(data[["RNA"]]$scale.data) else as.matrix(data[["RNA"]]@scale.data)

# run ScType
es.max <- sctype_score(scRNAseqData = scRNAseqData_scaled, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)

# merge by cluster
cL_resutls <- do.call("rbind", lapply(unique(data@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(data@meta.data[data@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(data@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores <- cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] <- "Unknown"
print(sctype_scores[,1:3])
write.csv(sctype_scores, file = "Celltypes.csv")

# Convert sctype_scores to a data frame if it isn't already
sctype_scores <- as.data.frame(sctype_scores)

# Extract the cluster and cell type columns
cluster_ids <- sctype_scores$cluster
cell_types <- sctype_scores$type

# Create a named vector for the new cluster identities
new_cluster_ids <- cell_types
names(new_cluster_ids) <- cluster_ids

# Rename the cluster identities in the object
data <- RenameIdents(data, new_cluster_ids)
levels(data)

# Overlay the identified cell types on UMAP plot
data@meta.data$sctype_classification = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  data@meta.data$sctype_classification[data@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}
```

### E. Visualization & Analysis

- Generates annotated UMAP/t-SNE plots
- Creates condition-split visualizations
- Counts cell distribution across types and conditions
- Saves annotated dataset

```r
# ==============================================================================
# Count Cells by Cell Type and Condition
# ==============================================================================
# Number of cells in each condition
# Create the table
num_of_cells <- table(data@meta.data$sctype_classification, data@meta.data$condition)
# Convert to data frame
num_of_cells_df <- as.data.frame.matrix(num_of_cells)
num_of_cells_df
# Save as CSV
write.csv(num_of_cells_df, file = "number_of_cells_in_conditions.csv", row.names = TRUE)

# ==============================================================================
# VISUALIZE ANNOTATED DATA
# ==============================================================================

# umap and tsne plots
umap_annt <- DimPlot(data, reduction = "umap", label = FALSE,
                     repel = TRUE, group.by = 'sctype_classification')        
umap_annt

tsne_annt <- DimPlot(data, reduction = "tsne", label = FALSE,
                     repel = TRUE, group.by = 'sctype_classification')        
tsne_annt

umap_annt_condition <- DimPlot(data, reduction = "umap", label = FALSE,
                               repel = TRUE, group.by = 'sctype_classification',
                               split.by = "condition")        
umap_annt_condition

tsne_annt_condition <- DimPlot(data, reduction = "tsne", label = FALSE,
                               repel = TRUE, group.by = 'sctype_classification',
                               split.by = "condition")        
tsne_annt_condition
```

---

### UMAP and t-SNE Annotated Cell Clusters

<div style="display:flex; justify-content:space-between;">
  <img src="data/10-umap_annotated.png" width="49%">
  <img src="data/12-tsne_annotated.png" width="49%">
</div>

---

### UMAP - Annotated Cell Clusters Grouped By Condition

![umap-annotated-condition](data/11-umap_annotated_condition.png)

---

### t-SNE - Annotated Cell Clusters Grouped By Condition

![tsne-annotated-condition](data/13-tsne_annotated_condition.png)

---

### F. Differential Expression Analysis

- Subsets specific cell types if needed

```r
# cells to keep
cells <- c("Memory CD4+ T cells",
           "Naive CD4+ T cells",
           "Cancer cells",
           "CD8+ NKT-like cells",
           "Natural killer  cells")
# Subset based on metadata column
data_subset <- subset(data, subset = sctype_classification %in% cells)
# Verify the subset
levels(data_subset)
table(data_subset@meta.data$sctype_classification)
table(data_subset@meta.data$condition)
saveRDS(data_subset, file = "data_subset.rds")
```

- Performs differential expression between:
- Same cell type across conditions (Normal vs PD)

```r
      data$celltype.condition <- paste(data$sctype_classification, data$condition, sep = "_")
Idents(data) <- "celltype.condition"
levels(data)

dff_markers <- FindMarkers(data, ident.1 = "Memory CD4+ T cells_LM",
                           ident.2 = "Memory CD4+ T cells_Normal", min.pct = 0.1,
                           test.use = "wilcox", verbose = TRUE, only.pos = FALSE)
# Apply FDR
dff_markers$fdr = p.adjust(dff_markers$p_val, method='fdr')
write.csv(dff_markers, file = "Memory CD4+ T cells_DEGs.csv") # Required for volcano plot

dff_up_markers <- subset(dff_markers, avg_log2FC > 0.25)
dff_dn_markers <- subset(dff_markers, avg_log2FC < -0.25)
# Save degs results
write.csv(dff_up_markers, file = "Memory CD4+ T cells_up_markers.csv")
write.csv(dff_dn_markers, file = "Memory CD4+ T cells_dn_markers.csv")
```
      
- Different cell types regardless of condition

```r
Idents(data) <- "sctype_classification"
levels(data)
degs <- FindMarkers(data, ident.1 = "Cancer cells",
                    ident.2 = NULL, test.use = "wilcox",
                    only.pos = FALSE)
# Apply FDR
degs$fdr = p.adjust(degs$p_val, method='fdr')
write.csv(degs, file = "Cancer cells_DEGs.csv") # Required for volcano plot

degs_up <- subset(degs, avg_log2FC > 0.25)
degs_down <- subset(degs, avg_log2FC < -0.25)
# Save degs results
write.csv(degs_up, file = "Cancer_cells_up_markers.csv")
write.csv(degs_down, file = "Cancer_cells_dn_markers.csv")
```

- Different cell types in different conditions

```r
data$celltype.condition <- paste(data$sctype_classification, data$condition, sep = "_")
Idents(data) <- "celltype.condition"
levels(data)
# DEGs of Cancer cells in diseased condition vs All other cell types in Normal condition
degs <- FindMarkers(data, ident.1 = "Cancer cells_LM",
                    ident.2 = c("Macrophages_Normal",
                                "Natural killer  cells_Normal",
                                "Basophils_Normal",
                                "Non-classical monocytes_Normal",
                                "CD8+ NKT-like cells_Normal",
                                "Memory CD4+ T cells_Normal",
                                "Myeloid Dendritic cells_Normal",
                                "Pre-B cells_Normal",
                                "Endothelial_Normal",
                                "Classical Monocytes_Normal",
                                "Plasmacytoid Dendritic cells_Normal",
                                "γδ-T cells_Normal",
                                "Memory B cells_Normal",
                                "Naive CD4+ T cells_Normal"
                    ), test.use = "wilcox",
                    only.pos = FALSE)
# Apply FDR
degs$fdr = p.adjust(degs$p_val, method='fdr')
write.csv(degs, file = "Cancer cells_DEGs.csv") # Required for volcano plot

degs_up <- subset(degs, avg_log2FC > 0.25)
degs_down <- subset(degs, avg_log2FC < -0.25)
# Save degs results
write.csv(degs_up, file = "Cancer_cells_up_markers.csv")
write.csv(degs_down, file = "Cancer_cells_dn_markers.csv")
```

- Identifies up/down-regulated genes
- Applies statistical filtering and FDR correction
- Exports marker gene lists

### Outputs

- The script produces:
    - Quality control plots
    - Clustering visualizations
    - Annotated cell type maps
    - Differential expression results
    - Processed datasets for downstream analysis (e.g Cell-cell communication, Trajectory analysis etc)

## 03. Multiple Datasets Integration

<p align="justify"> This script performs comprehensive integration of multiple single-cell RNA sequencing datasets to address technical variations while preserving biological signals. </p>

### Workflow

### A. Data Preparation

- Loads and combines samples from multiple experimental datasets
- Adds metadata labels for dataset origin and experimental conditions
- Merges all samples into a unified analysis object

### B. Quality Control & Preprocessing

- Clusters data without batch correction
- Generates low-dimensional embeddings (UMAP/t-SNE) showing raw cluster patterns
- Visualizes inherent batch effects and sample separation

### C. Batch Effect Correction

- Applies advanced integration methods to remove technical variability
- Creates harmonized dimensional space while preserving biological variation
- Maintains dataset-specific biological signals while minimizing batch effects

### E. Integrated Data Analysis

- Re-clusters cells using integrated dimensions
- Generates new embeddings on corrected data
- Creates comparative visualizations showing integration effectiveness

---

### UMAP - Before and After Batch-Effect Correction

![umap-harmony](data/14-umap-batch.png)

---

### t-SNE - Before and After Batch-Effect Correction

![tsne-harmony](data/15-tsne-batch.png)

### F. Cell Type Annotation

- Automates cell type identification using reference gene signatures
- Assigns cell types to clusters based on marker expression patterns
- Verifies annotation quality through multiple visualization approaches

### G. Data Exploration & Analysis

- Quantifies cell distribution across types and conditions
- Enables flexible subsetting for focused analysis of specific populations
- Exports fully annotated datasets for further investigation

### H. Differential Expression Analysis

- Performs comprehensive gene expression comparisons:
    - Within cell types across experimental conditions
    - Between different cell types
    - Complex multi-factor comparisons
- Identifies statistically significant differentially expressed genes

### Outputs
- This script produes:
    - Publication-ready visualizations
    - Fully annotated datasets
    - Comprehensive differential expression results

## 04. 3D UMAP and t-SNE Plots

This script creates interactive 3D visualizations of single-cell data using UMAP and t-SNE dimensionality reduction methods.

### Workflow

### A. 3D UMAP and t-SNE Generation

- Computes UMAP and t-SNE with three components for 3D coordinates
- Extracts embedding data and cluster information
- Creates interactive 3D scatter plots colored by cell clusters
- Exports visualizations in multiple formats (PNG, HTML)

### B. Gene Expression Visualization

- Maps gene expression patterns onto 3D UMAP and t-SNE space
- Uses `HLA-DQA1` gene as example with adjustable expression scaling
- Creates color-coded plots showing expression levels
- Generates interactive plots with cell-level information

### C. View 3D UMAP and t-SNE Plots

| Description | Link |
|---|---|
| **3D umap** | [3d-umap](https://usman4373.github.io/scRNA-seq-navigator/data/16-3d_umap.html) |
| **Gene Expression in 3D umap** | [3d-gene-umap](https://usman4373.github.io/scRNA-seq-navigator/data/17-gene-3d_umap.html) |
| **3D tsne** | [3d-tsne](https://usman4373.github.io/scRNA-seq-navigator/data/18-3d_tsne.html) |
| **Gene Expression in 3D tsne** | [3d-gene-tsne](https://usman4373.github.io/scRNA-seq-navigator/data/19-gene-3d_tsne.html) |

### Outputs

- The script produces:
    - Interactive 3D UMAP and t-SNE plots
    - Gene expression overlays in 3D space
    - Multiple export formats (static PNG, interactive HTML)
    - Enhanced visualization capabilities for data exploration

## 05. Volcano Plot

<p align="justify"> This script creates volcano plots of specific cell types to visualize differentially expressed genes between experimental conditions. </p>

### Workflow

### A. Data Preparation

- Loads differential expression results from CSV files
- Formats data for volcano plot visualization
- Sets up gene identifiers and statistical columns

### B. Plot Generation

- Uses EnhancedVolcano package for advanced plotting
- Maps log2 fold changes against adjusted p-values
- Applies statistical cutoffs for significance
- Customizes colors for different significance categories

### C. Visualization & Export

- Creates publication-quality volcano plots
- Highlights significantly up/down-regulated genes
- Saves high-resolution images for analysis

![volcano](data/20-volcano.png)

---

## 06. Heatmap

This script generates gene-specific heatmaps to visualize gene expression patterns across cell types and conditions.

### Workflow

### A. Data Extraction

- Selects genes of interest for visualization
- Extracts expression data from Seurat object
- Combines expression data with metadata

### B. Expression Analysis

- Calculates average expression per cell type and condition
- Creates expression matrices for heatmap plotting
- Saves expression data for documentation

### C. Plot Generation

- Prepares annotation colors for conditions and cell types
- Uses pheatmap for clustered visualization
- Applies color schemes and formatting
- Exports high-quality heatmap images

<div style="display:flex; justify-content:space-between;">
  <img src="data/21-heatmap.png" width="49%">
  <img src="data/22-heatmap.png" width="49%">
</div>

---

## 07. Boxplot

This script creates box plots to compare gene expression across conditions and specific cell types.

### Workflow

### Condition-Based Comparison

- Extracts expression data for specific genes
- Compares expression between Control and Diseased conditions
- Creates box plots with jittered points
- Applies custom color schemes

### Cell Type-Specific Analysis

- Filters data for specific cell types
- Compares expression within cell types across conditions
- Generates focused visualizations

### Visualization Customization

- Uses ggplot2 for advanced plotting
- Applies consistent theming and styling
- Creates publication-ready figures

<div style="display:flex; justify-content:space-between;">
  <img src="data/23-boxplot.png" width="49%">
  <img src="data/24-boxplot.png" width="49%">
</div>

---

## 📚 Citation

- If you use this analysis pipeline in your work, please cite the repository:

```
scRNA-seq-navigator. GitHub: https://github.com/usman4373/scRNA-seq-navigator
```

## 🤝 Acknowledgements

- <p align="justify"> This pipeline implements established single-cell analysis methods, including <a href="https://doi.org/10.1016/j.cell.2019.05.031">Seurat (Stuart et al., 2019)</a>, <a href="https://doi.org/10.1038/s41592-019-0619-0">Harmony (Korsunsky et al., 2019)</a>, and <a href="https://doi.org/10.1038/s41467-022-28803-w">ScType (Ianevski et al., 2022)</a>. Please also cite the original publications of these tools when using their respective functionalities in your research. </p>
- Please also cite other tools used in the analysis pipeline:
    - dplyr - For efficient data manipulation and transformation
    - patchwork - For seamless arrangement of multiple ggplot2 plots
    - ggplot2 - For creating elegant and customizable data visualizations
    - openxlsx - For reading and writing Excel files in R
    - HGNChelper - For handling and validating HGNC gene symbols
    - plotly - For creating interactive web-based visualizations
    - ComplexHeatmap - For generating sophisticated and annotated heatmaps
    - circlize - For circular visualization and heatmap enhancements
    - paletteer - For accessing comprehensive color palette collections
    - EnhancedVolcano - For publication-ready volcano plot generation
    - viridis - For colorblind-friendly color scales in visualizations
    - htmlwidgets - For embedding interactive JavaScript visualizations in R
    - webshot2 - For converting interactive plots to static images
    - readr - For fast and tidy data import functionality
    - Other open-source R packages and dependencies

---
