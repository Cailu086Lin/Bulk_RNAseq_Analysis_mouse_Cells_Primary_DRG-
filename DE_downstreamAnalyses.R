###############################################################################
# Mouse neuronal bulk RNA-seq analysis
#
# A1–A3 = ESP treatment
# B1–B3 = SGE treatment
# C1–C3 = Control
#
# Input:
#   Salmon GRCm39 decoy-aware quantifications
#
# Analyses:
#   tximport
#   DESeq2
#   PCA and clustering QC
#   Differential expression
#   Heatmaps
#   Volcano plots
#   Gene Ontology enrichment
###############################################################################

# =============================================================================
# 1. Packages
# =============================================================================

library(tximport)
library(DESeq2)
library(rtracklayer)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(clusterProfiler)

library(dplyr)
library(tibble)
library(tidyr)
library(readr)
library(ggplot2)
library(ggrepel)
library(pheatmap)

# Optional, for publication heatmaps
library(ComplexHeatmap)
library(circlize)
library(grid)
library(EnhancedVolcano)

# =============================================================================
# 2. Project directories
# =============================================================================

project_dir <- "/mnt/easystore/DBH_26/Stephenson"

quant_dir <- file.path(
  project_dir,
  "salmon_GRCm39_decoy_quant_idfixed"
)

ref_dir <- file.path(
  project_dir,
  "reference/GRCm39"
)

# Change this filename if your GTF has a different name
gtf_file <- file.path(
  ref_dir,
  "gencode.vM39.annotation.gtf"
)

results_dir <- file.path(
  project_dir,
  "downstream_GRCm39_decoy"
)

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

dir.create(
  file.path(results_dir, "DEG"),
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  file.path(results_dir, "Figures"),
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  file.path(results_dir, "GO"),
  recursive = TRUE,
  showWarnings = FALSE
)


# =============================================================================
# 3. Sample information
# =============================================================================

sample_table <- data.frame(
  sample = c(
    "A1", "A2", "A3",
    "B1", "B2", "B3",
    "C1", "C2", "C3"
  ),
  
  condition = c(
    rep("ESP", 3),
    rep("SGE", 3),
    rep("Control", 3)
  ),
  
  replicate = rep(1:3, 3),
  
  stringsAsFactors = FALSE
)

sample_table$condition <- factor(
  sample_table$condition,
  levels = c("Control", "ESP", "SGE")
)

rownames(sample_table) <- sample_table$sample

sample_table
table(sample_table$condition)

write.csv(
  sample_table,
  file.path(results_dir, "sample_table.csv"),
  row.names = FALSE
)


# =============================================================================
# 4. Salmon quantification files
# =============================================================================

files <- file.path(
  quant_dir,
  paste0(sample_table$sample, "_quant"),
  "quant.sf"
)

names(files) <- sample_table$sample

file_check <- data.frame(
  sample = names(files),
  file = unname(files),
  exists = file.exists(files)
)

print(file_check)

if (!all(file_check$exists)) {
  
  missing_files <- file_check$file[!file_check$exists]
  
  stop(
    paste(
      "Missing Salmon quant.sf files:\n",
      paste(missing_files, collapse = "\n")
    )
  )
}


# =============================================================================
# 5. Create transcript-to-gene mapping from GENCODE M39 GTF
# =============================================================================

if (!file.exists(gtf_file)) {
  stop("GTF file was not found: ", gtf_file)
}

gtf <- rtracklayer::import(gtf_file)

tx2gene <- as.data.frame(gtf) %>%
  filter(
    type == "transcript",
    !is.na(transcript_id),
    !is.na(gene_id)
  ) %>%
  transmute(
    transcript_id = sub("\\.[0-9]+$", "", transcript_id),
    gene_id = sub("\\.[0-9]+$", "", gene_id),
    gene_name = as.character(gene_name)
  ) %>%
  distinct(transcript_id, gene_id, .keep_all = TRUE)

head(tx2gene)
dim(tx2gene)
sapply(tx2gene, class)

write.csv(
  tx2gene,
  file.path(results_dir, "GRCm39_M39_tx2gene.csv"),
  row.names = FALSE
)


# =============================================================================
# 6. Import Salmon data
# =============================================================================

txi <- tximport(
  files = files,
  type = "salmon",
  tx2gene = tx2gene[, c("transcript_id", "gene_id")],
  countsFromAbundance = "lengthScaledTPM",
  ignoreTxVersion = TRUE,
  dropInfReps = TRUE
)

dim(txi$counts)

head(txi$counts[, 1:3])

stopifnot(
  identical(
    colnames(txi$counts),
    rownames(sample_table)
  )
)


# =============================================================================
# 7. Create DESeq2 object
# =============================================================================

# No batch variable was provided, so the appropriate design is condition only.
dds <- DESeqDataSetFromTximport(
  txi = txi,
  colData = sample_table,
  design = ~ condition
)

# Keep genes with at least 10 estimated counts in at least 3 samples
keep <- rowSums(counts(dds) >= 10) >= 3

table(keep)

dds <- dds[keep, ]

dds <- DESeq(dds)

saveRDS(
  dds,
  file.path(results_dir, "dds_GRCm39_decoy.rds")
)


# =============================================================================
# 8. Normalized expression tables
# =============================================================================

normalized_counts <- counts(
  dds,
  normalized = TRUE
)

write.csv(
  data.frame(
    gene_id = rownames(normalized_counts),
    normalized_counts,
    check.names = FALSE
  ),
  file.path(results_dir, "normalized_counts.csv"),
  row.names = FALSE
)

# Variance-stabilized expression for PCA and heatmaps
vsd <- vst(
  dds,
  blind = FALSE
)

vsd_mat <- assay(vsd)

write.csv(
  data.frame(
    gene_id = rownames(vsd_mat),
    vsd_mat,
    check.names = FALSE
  ),
  file.path(results_dir, "VST_expression.csv"),
  row.names = FALSE
)
vsd <- varianceStabilizingTransformation(dds, blind = TRUE)
plotPCA(vsd, intgroup = "condition")

# =============================================================================
# 9. Gene annotation
# =============================================================================

gene_annotation <- tx2gene %>%
  select(gene_id, gene_name) %>%
  filter(!is.na(gene_name), gene_name != "") %>%
  distinct(gene_id, .keep_all = TRUE)

# Add Entrez IDs for GO analysis
entrez_map <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys = unique(gene_annotation$gene_id),
  keytype = "ENSEMBL",
  columns = c("SYMBOL", "ENTREZID")
) %>%
  as.data.frame() %>%
  rename(
    gene_id = ENSEMBL,
    annotation_symbol = SYMBOL,
    entrez_id = ENTREZID
  ) %>%
  distinct(gene_id, .keep_all = TRUE)

gene_annotation <- gene_annotation %>%
  left_join(entrez_map, by = "gene_id") %>%
  mutate(
    gene_symbol = ifelse(
      is.na(gene_name) | gene_name == "",
      annotation_symbol,
      gene_name
    )
  )

write.csv(
  gene_annotation,
  file.path(results_dir, "mouse_gene_annotation.csv"),
  row.names = FALSE
)


# =============================================================================
# 10. Library-size plot
# =============================================================================

library_sizes <- data.frame(
  sample = colnames(txi$counts),
  estimated_counts = colSums(txi$counts)
) %>%
  left_join(sample_table, by = "sample")

p_library <- ggplot(
  library_sizes,
  aes(
    x = sample,
    y = estimated_counts / 1e6,
    fill = condition
  )
) +
  geom_col(width = 0.75) +
  labs(
    title = "Estimated library size",
    x = NULL,
    y = "Estimated counts (millions)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "top"
  )

ggsave(
  file.path(results_dir, "Figures", "library_size.pdf"),
  p_library,
  width = 7,
  height = 5
)

p_library


# =============================================================================
# 11. PCA
# =============================================================================

pca_data <- plotPCA(
  vsd,
  intgroup = "condition",
  returnData = TRUE
)

percent_var <- round(
  100 * attr(pca_data, "percentVar"),
  1
)

p_pca <- ggplot(
  pca_data,
  aes(
    x = PC1,
    y = PC2,
    color = condition,
    label = name
  )
) +
  geom_point(size = 4) +
  geom_text_repel(
    size = 4,
    max.overlaps = Inf
  ) +
  labs(
    title = "PCA: mouse neuronal cultures",
    x = paste0("PC1: ", percent_var[1], "% variance"),
    y = paste0("PC2: ", percent_var[2], "% variance"),
    color = "Condition"
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "top"
  )

ggsave(
  file.path(results_dir, "Figures", "PCA_all_samples.pdf"),
  p_pca,
  width = 7,
  height = 6
)

p_pca


# =============================================================================
# 12. Sample-to-sample distance heatmap
# =============================================================================

sample_dist <- dist(t(vsd_mat))

sample_dist_matrix <- as.matrix(sample_dist)

annotation_col <- data.frame(
  Condition = sample_table[colnames(sample_dist_matrix), "condition"]
)

rownames(annotation_col) <- colnames(sample_dist_matrix)

pdf(
  file.path(
    results_dir,
    "Figures",
    "sample_distance_heatmap.pdf"
  ),
  width = 7,
  height = 6
)

pheatmap(
  sample_dist_matrix,
  annotation_col = annotation_col,
  annotation_row = annotation_col,
  main = "Sample-to-sample distances",
  border_color = NA
)

dev.off()


# =============================================================================
# 13. Differential-expression contrasts
# =============================================================================

comparison_list <- list(
  ESP_vs_Control = c("ESP", "Control"),
  SGE_vs_Control = c("SGE", "Control"),
  ESP_vs_SGE = c("ESP", "SGE")
)

extract_deseq_result <- function(
    dds_object,
    numerator,
    denominator,
    comparison_name,
    annotation_table,
    padj_cutoff = 0.05,
    lfc_cutoff = 1
) {
  
  res <- results(
    dds_object,
    contrast = c(
      "condition",
      numerator,
      denominator
    ),
    alpha = padj_cutoff
  )
  
  res_df <- as.data.frame(res) %>%
    rownames_to_column("gene_id") %>%
    left_join(
      annotation_table,
      by = "gene_id"
    ) %>%
    mutate(
      Comparison = comparison_name,
      numerator = numerator,
      denominator = denominator,
      
      direction = case_when(
        !is.na(padj) &
          padj < padj_cutoff &
          log2FoldChange >= lfc_cutoff ~
          paste0("Up_in_", numerator),
        
        !is.na(padj) &
          padj < padj_cutoff &
          log2FoldChange <= -lfc_cutoff ~
          paste0("Down_in_", numerator),
        
        TRUE ~ "NS"
      )
    ) %>%
    arrange(padj)
  
  res_df
}


deg_list <- list()

for (comparison_name in names(comparison_list)) {
  
  numerator <- comparison_list[[comparison_name]][1]
  denominator <- comparison_list[[comparison_name]][2]
  
  deg_list[[comparison_name]] <- extract_deseq_result(
    dds_object = dds,
    numerator = numerator,
    denominator = denominator,
    comparison_name = comparison_name,
    annotation_table = gene_annotation
  )
  
  write.csv(
    deg_list[[comparison_name]],
    file.path(
      results_dir,
      "DEG",
      paste0(comparison_name, "_all_genesNoB1.csv")
    ),
    row.names = FALSE
  )
}

deg_all <- bind_rows(deg_list)

write.csv(
  deg_all,
  file.path(results_dir, "DEG", "all_comparisons.csv"),
  row.names = FALSE
)


# =============================================================================
# 14. DEG summary
# =============================================================================

deg_summary <- deg_all %>%
  group_by(Comparison, direction) %>%
  summarise(
    n_genes = n(),
    .groups = "drop"
  ) %>%
  arrange(Comparison, direction)

print(deg_summary)

write.csv(
  deg_summary,
  file.path(results_dir, "DEG", "DEG_summary.csv"),
  row.names = FALSE
)

# Summary using padj < 0.05 and |log2FC| >= 1
significant_deg_summary <- deg_all %>%
  filter(
    !is.na(padj),
    padj < 0.05,
    abs(log2FoldChange) >= 1
  ) %>%
  group_by(Comparison) %>%
  summarise(
    up = sum(log2FoldChange >= 1),
    down = sum(log2FoldChange <= -1),
    total = n(),
    .groups = "drop"
  )

print(significant_deg_summary)

write.csv(
  significant_deg_summary,
  file.path(results_dir, "DEG", "significant_DEG_summary.csv"),
  row.names = FALSE
)


# =============================================================================
# 15. Volcano plots
# =============================================================================
make_volcano <- function(
    result_df,
    title,
    padj_cutoff = 0.05,
    lfc_cutoff = 1,
    max_labels = 20
){
  
  plot_df <- result_df %>%
    mutate(
      plot_padj = ifelse(is.na(padj), 1, padj),
      plot_padj = pmax(plot_padj, .Machine$double.xmin),
      
      label = case_when(
        !is.na(gene_symbol) & gene_symbol != "" ~ gene_symbol,
        !is.na(gene_id) & gene_id != "" ~ gene_id,
        TRUE ~ rownames(result_df)
      )
    ) %>%
    filter(
      !is.na(log2FoldChange),
      is.finite(log2FoldChange)
    )
  
  ## Top labels from both directions
  n_up   <- ceiling(max_labels/2)
  n_down <- floor(max_labels/2)
  
  selected_labels <- bind_rows(
    
    plot_df %>%
      filter(
        plot_padj < padj_cutoff,
        log2FoldChange >= lfc_cutoff
      ) %>%
      arrange(plot_padj) %>%
      slice_head(n = n_up),
    
    plot_df %>%
      filter(
        plot_padj < padj_cutoff,
        log2FoldChange <= -lfc_cutoff
      ) %>%
      arrange(plot_padj) %>%
      slice_head(n = n_down)
    
  ) %>%
    pull(label) %>%
    unique()
  
  EnhancedVolcano(
    plot_df,
    
    lab = plot_df$label,
    
    x = "log2FoldChange",
    y = "plot_padj",
    
    title = title,
    subtitle = NULL,
    caption = NULL,
    
    pCutoff = padj_cutoff,
    FCcutoff = lfc_cutoff,
    
    selectLab = selected_labels,
    
    drawConnectors = TRUE,
    widthConnectors = 0.4,
    labSize = 3.5,
    pointSize = 2,
    
    legendPosition = "none",
    
    gridlines.major = FALSE,
    gridlines.minor = FALSE,
    
    col = c(
      "grey80",   # NS
      "grey80",# FC only
      "grey80",# p only
      "red3"      # significant
    ),
    
    border = "partial"
  )
}
dir.create(
  file.path(results_dir, "Figures"),
  recursive = TRUE,
  showWarnings = FALSE
)

for (comparison_name in names(deg_list)) {
  
  p_volcano <- make_volcano(
    result_df = deg_list[[comparison_name]],
    title = gsub("_", " ", comparison_name),
    padj_cutoff = 0.05,
    lfc_cutoff = 1,
    max_labels = 15
  )
  
  ggsave(
    filename = file.path(
      results_dir,
      "Figures",
      paste0("Volcano_", comparison_name, ".pdf")
    ),
    plot = p_volcano,
    width = 7,
    height = 6,
    device = cairo_pdf
  )
}



# =============================================================================
# 16. Heatmap of the most variable genes
# =============================================================================
# ============================================================
# Top 50 most variable genes, excluding genes beginning with Gm
# ============================================================

# Remove Ensembl version suffixes, if present
matrix_gene_ids <- sub(
  "\\.[0-9]+$",
  "",
  rownames(vsd_mat)
)

annotation_gene_ids <- sub(
  "\\.[0-9]+$",
  "",
  gene_annotation$gene_id
)

# Match Ensembl IDs to gene symbols
all_gene_symbols <- gene_annotation$gene_symbol[
  match(
    matrix_gene_ids,
    annotation_gene_ids
  )
]

# Replace missing or empty symbols with Ensembl gene IDs
missing_symbol <- is.na(all_gene_symbols) |
  trimws(all_gene_symbols) == ""

all_gene_symbols[missing_symbol] <-
  matrix_gene_ids[missing_symbol]

# Create a visualization copy with gene-symbol row names
# make.unique() prevents duplicated row names
vsd_mat_symbol <- vsd_mat
rownames(vsd_mat_symbol) <- make.unique(all_gene_symbols)

# Calculate variance for every gene
gene_variance <- apply(
  vsd_mat_symbol,
  1,
  var,
  na.rm = TRUE
)

# Retain genes with finite variance and symbols not beginning with Gm
keep_genes <- is.finite(gene_variance) &
  !grepl(
    "^Gm|^ENSMUSG",
    rownames(vsd_mat_symbol),
    ignore.case = TRUE
  )

gene_variance_filtered <- gene_variance[keep_genes]

# Select the top 50 genes after filtering
n_top <- min(
  50,
  length(gene_variance_filtered)
)

top_variable_genes <- names(
  sort(
    gene_variance_filtered,
    decreasing = TRUE
  )
)[seq_len(n_top)]

# Subset the gene-symbol matrix
top_var_mat <- vsd_mat_symbol[
  top_variable_genes,
  ,
  drop = FALSE
]

# Final check: should return character(0)
grep(
  "^Gm",
  rownames(top_var_mat),
  value = TRUE,
  ignore.case = TRUE
)

# Row-wise Z-score
top_var_z <- t(
  scale(
    t(top_var_mat)
  )
)

# Remove genes with undefined Z-scores
top_var_z <- top_var_z[
  complete.cases(top_var_z),
  ,
  drop = FALSE
]

# Match sample annotation by sample name
annotation_col <- data.frame(
  Condition = sample_table[
    match(
      colnames(top_var_z),
      sample_table$sample
    ),
    "condition"
  ],
  row.names = colnames(top_var_z)
)

# Confirm annotation matched correctly
stopifnot(
  !anyNA(annotation_col$Condition)
)

# Create output directory if needed
dir.create(
  file.path(results_dir, "Figures"),
  recursive = TRUE,
  showWarnings = FALSE
)

# Save heatmap
pdf(
  file.path(
    results_dir,
    "Figures",
    "heatmap_top50_variable_genes_no_Gm.pdf"
  ),
  width = 6,
  height = 9.5
)

pheatmap(
  top_var_z,
  #annotation_col = annotation_col,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 10,
  main = "Top 50 variable genes",
  border_color = NA
)

dev.off()

# =============================================================================
# 17. Heatmap of significant genes for each comparison
# =============================================================================

make_deg_heatmap <- function(
    comparison_name,
    deg_table,
    expression_matrix,
    sample_metadata,
    output_file,
    max_genes = 50,
    padj_cutoff = 0.05,
    lfc_cutoff = 1
) {
  
  # Check required columns
  required_cols <- c(
    "gene_id",
    "gene_symbol",
    "padj",
    "log2FoldChange"
  )
  
  missing_cols <- setdiff(required_cols, colnames(deg_table))
  
  if (length(missing_cols) > 0) {
    stop(
      comparison_name,
      ": missing columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  # Remove Ensembl version suffixes for matching
  deg_table <- deg_table %>%
    mutate(
      gene_id_clean = sub("\\.[0-9]+$", "", as.character(gene_id))
    )
  
  matrix_gene_ids <- sub(
    "\\.[0-9]+$",
    "",
    rownames(expression_matrix)
  )
  
  # Select significant genes
  selected <- deg_table %>%
    filter(
      !is.na(padj),
      is.finite(padj),
      padj < padj_cutoff,
      !is.na(log2FoldChange),
      is.finite(log2FoldChange),
      abs(log2FoldChange) >= lfc_cutoff
    ) %>%
    arrange(padj, desc(abs(log2FoldChange))) %>%
    distinct(gene_id_clean, .keep_all = TRUE) %>%
    slice_head(n = max_genes)
  
  message(
    comparison_name,
    ": ",
    nrow(selected),
    " significant genes before expression-matrix matching."
  )
  
  if (nrow(selected) < 2) {
    message(
      comparison_name,
      ": fewer than two significant genes; heatmap skipped."
    )
    return(NULL)
  }
  
  # Match selected DEG IDs to rows of the expression matrix
  matrix_index <- match(
    selected$gene_id_clean,
    matrix_gene_ids
  )
  
  matched <- !is.na(matrix_index)
  
  message(
    comparison_name,
    ": ",
    sum(matched),
    " genes matched to expression matrix."
  )
  
  selected <- selected[matched, , drop = FALSE]
  matrix_index <- matrix_index[matched]
  
  if (length(matrix_index) < 2) {
    message(
      comparison_name,
      ": fewer than two genes matched the expression matrix; ",
      "heatmap skipped."
    )
    return(NULL)
  }
  
  # Extract expression values in the same order as selected genes
  mat <- expression_matrix[
    matrix_index,
    ,
    drop = FALSE
  ]
  
  # Use gene symbols for display
  symbols <- as.character(selected$gene_symbol)
  
  missing_symbol <- is.na(symbols) |
    trimws(symbols) == ""
  
  symbols[missing_symbol] <-
    selected$gene_id_clean[missing_symbol]
  
  rownames(mat) <- make.unique(symbols)
  
  # Remove genes with zero or undefined variance
  row_variance <- apply(
    mat,
    1,
    var,
    na.rm = TRUE
  )
  
  keep_variable <- is.finite(row_variance) &
    row_variance > 0
  
  mat <- mat[
    keep_variable,
    ,
    drop = FALSE
  ]
  
  if (nrow(mat) < 2) {
    message(
      comparison_name,
      ": fewer than two genes have nonzero variance; heatmap skipped."
    )
    return(NULL)
  }
  
  # Row-wise Z-score
  mat_z <- t(
    scale(
      t(mat)
    )
  )
  
  mat_z <- mat_z[
    apply(mat_z, 1, function(x) all(is.finite(x))),
    ,
    drop = FALSE
  ]
  
  if (nrow(mat_z) < 2) {
    message(
      comparison_name,
      ": fewer than two genes remain after scaling; heatmap skipped."
    )
    return(NULL)
  }
  
  # Match sample metadata
  if ("sample" %in% colnames(sample_metadata)) {
    
    metadata_index <- match(
      colnames(mat_z),
      sample_metadata$sample
    )
    
    anno <- data.frame(
      Condition = sample_metadata$condition[metadata_index],
      row.names = colnames(mat_z)
    )
    
  } else {
    
    metadata_index <- match(
      colnames(mat_z),
      rownames(sample_metadata)
    )
    
    anno <- data.frame(
      Condition = sample_metadata$condition[metadata_index],
      row.names = colnames(mat_z)
    )
  }
  
  if (anyNA(anno$Condition)) {
    stop(
      comparison_name,
      ": sample metadata did not match these expression columns: ",
      paste(
        colnames(mat_z)[is.na(anno$Condition)],
        collapse = ", "
      )
    )
  }
  
  anno$Condition <- factor(anno$Condition)
  
  dir.create(
    dirname(output_file),
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  pdf(
    output_file,
    width = 8,
    height = max(6, nrow(mat_z) * 0.18)
  )
  # Report matrix status before plotting
  message(
    comparison_name,
    ": matrix dimension = ",
    nrow(mat_z), " x ", ncol(mat_z),
    "; range = ",
    paste(round(range(mat_z, finite = TRUE), 3), collapse = " to ")
  )
  
  if (
    nrow(mat_z) < 2 ||
    ncol(mat_z) < 2 ||
    !any(is.finite(mat_z))
  ) {
    message(comparison_name, ": invalid heatmap matrix; skipped.")
    return(NULL)
  }
  
  dir.create(
    dirname(output_file),
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  # Build heatmap without opening an internal device
  ph <- pheatmap::pheatmap(
    mat_z,
    annotation_col = anno,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize_row = 7,
    fontsize_col = 10,
    main = gsub("_", " ", comparison_name),
    border_color = NA,
    silent = TRUE
  )
  
  # Explicitly draw the heatmap into the PDF
  grDevices::cairo_pdf(
    filename = output_file,
    width = 8,
    height = max(6, nrow(mat_z) * 0.18)
  )
  
  grid::grid.newpage()
  grid::grid.draw(ph$gtable)
  
  grDevices::dev.off()
  
  message(
    comparison_name,
    ": heatmap saved with ",
    nrow(mat_z),
    " genes to ",
    output_file
  )
  
  invisible(mat_z)
}


dir.create(
  file.path(results_dir, "Figures"),
  recursive = TRUE,
  showWarnings = FALSE
)

graphics.off()

for (comparison_name in names(deg_list)) {
  
  make_deg_heatmap(
    comparison_name = comparison_name,
    deg_table = deg_list[[comparison_name]],
    expression_matrix = vsd_mat,
    sample_metadata = sample_table,
    output_file = file.path(
      results_dir,
      "Figures",
      paste0("Heatmap_", comparison_name, ".pdf")
    ),
    max_genes = 50,
    padj_cutoff = 0.05,
    lfc_cutoff = 1
  )
}

# =============================================================================
# 18. GO enrichment
# =============================================================================

# Background universe should contain genes tested by DESeq2
background_gene_ids <- rownames(dds)

background_entrez <- gene_annotation %>%
  filter(
    gene_id %in% background_gene_ids,
    !is.na(entrez_id)
  ) %>%
  pull(entrez_id) %>%
  unique()


run_go <- function(
    deg_table,
    comparison_name,
    direction_name,
    lfc_direction = c("up", "down"),
    universe_entrez,
    output_directory
) {
  
  lfc_direction <- match.arg(lfc_direction)
  
  if (lfc_direction == "up") {
    
    selected <- deg_table %>%
      filter(
        !is.na(padj),
        padj < 0.05,
        log2FoldChange >= 1,
        !is.na(entrez_id)
      )
    
  } else {
    
    selected <- deg_table %>%
      filter(
        !is.na(padj),
        padj < 0.05,
        log2FoldChange <= -1,
        !is.na(entrez_id)
      )
  }
  
  genes_entrez <- unique(selected$entrez_id)
  
  if (length(genes_entrez) < 5) {
    
    message(
      comparison_name,
      " ",
      direction_name,
      ": fewer than five mapped genes; GO skipped."
    )
    
    return(NULL)
  }
  
  ego <- enrichGO(
    gene = genes_entrez,
    universe = universe_entrez,
    OrgDb = org.Mm.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.20,
    readable = TRUE
  )
  
  go_df <- as.data.frame(ego)
  
  write.csv(
    go_df,
    file.path(
      output_directory,
      paste0(
        comparison_name,
        "_",
        direction_name,
        "_GO_BP.csv"
      )
    ),
    row.names = FALSE
  )
  
  if (nrow(go_df) > 0) {
    
    p_dot <- dotplot(
      ego,
      showCategory = 15
    ) +
      ggtitle(
        paste(
          gsub("_", " ", comparison_name),
          direction_name,
          "GO Biological Process"
        )
      )
    
    ggsave(
      file.path(
        output_directory,
        paste0(
          comparison_name,
          "_",
          direction_name,
          "_GO_BP_dotplot.pdf"
        )
      ),
      p_dot,
      width = 9,
      height = 6
    )
  }
  
  ego
}


go_results <- list()

for (comparison_name in names(deg_list)[-2]) {
  
  go_results[[paste0(comparison_name, "_Up")]] <- run_go(
    deg_table = deg_list[[comparison_name]],
    comparison_name = comparison_name,
    direction_name = "Up",
    lfc_direction = "up",
    universe_entrez = background_entrez,
    output_directory = file.path(results_dir, "GO")
  )
  
  go_results[[paste0(comparison_name, "_Down")]] <- run_go(
    deg_table = deg_list[[comparison_name]],
    comparison_name = comparison_name,
    direction_name = "Down",
    lfc_direction = "down",
    universe_entrez = background_entrez,
    output_directory = file.path(results_dir, "GO")
  )
}


# =============================================================================
# 19. Optional: heatmap of top treatment-responsive genes across all samples
# =============================================================================

all_sig_genes <- deg_all %>%
  filter(
    !is.na(padj),
    padj < 0.05,
    abs(log2FoldChange) >= 1
  ) %>%
  arrange(padj) %>%
  distinct(gene_id, .keep_all = TRUE) %>%
  slice_head(n = 75)

genes_use <- intersect(
  all_sig_genes$gene_id,
  rownames(vsd_mat)
)

if (length(genes_use) >= 2) {
  
  response_mat <- vsd_mat[
    genes_use,
    ,
    drop = FALSE
  ]
  
  response_symbols <- all_sig_genes$gene_symbol[
    match(
      rownames(response_mat),
      all_sig_genes$gene_id
    )
  ]
  
  response_symbols[
    is.na(response_symbols) |
      response_symbols == ""
  ] <- rownames(response_mat)[
    is.na(response_symbols) |
      response_symbols == ""
  ]
  
  rownames(response_mat) <- make.unique(response_symbols)
  
  response_z <- t(scale(t(response_mat)))
  
  response_z <- response_z[
    complete.cases(response_z),
    ,
    drop = FALSE
  ]
  
  response_anno <- data.frame(
    Condition = sample_table[
      colnames(response_z),
      "condition"
    ]
  )
  
  rownames(response_anno) <- colnames(response_z)
  
  pdf(
    file.path(
      results_dir,
      "Figures",
      "heatmap_top_treatment_responsive_genes.pdf"
    ),
    width = 8,
    height = 10
  )
  
  pheatmap::pheatmap(
    response_z,
    annotation_col = response_anno,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    fontsize_row = 6,
    main = "Top treatment-responsive genes",
    border_color = NA
  )
  
  dev.off()
}


# =============================================================================
# 20. Save session information
# =============================================================================

capture.output(
  sessionInfo(),
  file = file.path(results_dir, "sessionInfo.txt")
)

cat(
  "\nAnalysis completed.\nResults directory:\n",
  results_dir,
  "\n"
)
