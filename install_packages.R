#!/usr/bin/env Rscript
# ==============================================================================
# INSTALL REQUIRED R PACKAGES FOR SCRNA-SEQ ANALYSIS PIPELINE
# ==============================================================================

# ==============================================================================
# CRAN PACKAGES
# ==============================================================================

cat("Installing CRAN packages...\n")

cran_packages <- c(
  "Seurat",
  "dplyr", 
  "patchwork",
  "ggplot2",
  "plotly",
  "readr",
  "openxlsx",
  "paletteer",
  "htmlwidgets",
  "webshot2",
  "viridis",
  "harmony",
  "circlize",
  "ComplexHeatmap"
)

for (pkg in cran_packages) {
  cat("Installing", pkg, "...\n")
  tryCatch({
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      install.packages(pkg, repos = "https://cloud.r-project.org/")
      cat("✓", pkg, "installed successfully\n")
    } else {
      cat("✓", pkg, "already installed\n")
    }
  }, error = function(e) {
    cat("✗ ERROR installing", pkg, ":", e$message, "\n")
  })
}

# ==============================================================================
# BIOCONDUCTOR PACKAGES
# ==============================================================================

cat("\nInstalling Bioconductor packages...\n")

bioc_packages <- c("HGNChelper", "EnhancedVolcano")

# Install BiocManager if not already installed
if (!require("BiocManager", quietly = TRUE)) {
  cat("Installing BiocManager...\n")
  tryCatch({
    install.packages("BiocManager", repos = "https://cloud.r-project.org/")
    cat("✓ BiocManager installed successfully\n")
  }, error = function(e) {
    cat("✗ ERROR installing BiocManager:", e$message, "\n")
  })
}

# Install Bioconductor packages
for (pkg in bioc_packages) {
  cat("Installing", pkg, "...\n")
  tryCatch({
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      BiocManager::install(pkg, ask = FALSE)
      cat("✓", pkg, "installed successfully\n")
    } else {
      cat("✓", pkg, "already installed\n")
    }
  }, error = function(e) {
    cat("✗ ERROR installing", pkg, ":", e$message, "\n")
  })
}

# ==============================================================================
# GITHUB PACKAGES (for scType)
# ==============================================================================

cat("\nLoading scType functions (will be sourced during analysis)...\n")
cat("Note: scType functions are sourced directly from GitHub during analysis\n")
cat("No installation required for scType\n")

# ==============================================================================
# OPTIONAL/SPECIALIZED PACKAGES
# ==============================================================================

cat("\nOptional packages (commented out in scripts)...\n")

optional_packages <- c("SeuratDisk", "BPCells", "orca", "pheatmap")

for (pkg in optional_packages) {
  cat("Optional package available:", pkg, "\n")
}

# ==============================================================================
# VERIFICATION
# ==============================================================================

cat("\nVerifying installations...\n")

all_packages <- c(cran_packages, bioc_packages)

missing_packages <- c()
for (pkg in all_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    missing_packages <- c(missing_packages, pkg)
  }
}

if (length(missing_packages) == 0) {
  cat("✓ All required packages installed successfully!\n")
} else {
  cat("✗ The following packages failed to install:\n")
  for (pkg in missing_packages) {
    cat("  -", pkg, "\n")
  }
  quit(status = 1)
}

cat("\nInstallation completed successfully!\n")
cat("You can now run the scRNA-seq analysis scripts.\n")