## ============================================================
## GO Biological Process enrichment with redundancy reduction
## ============================================================

library(dplyr)
library(clusterProfiler)
library(enrichplot)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(ggplot2)

## ------------------------------------------------------------
## 1. Create output directory
## ------------------------------------------------------------

go_output_dir <- file.path(results_dir, "GO")

dir.create(
  go_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


## ------------------------------------------------------------
## 2. Prepare the background gene universe
##
## Assumption:
## rownames(vsd_mat) contain mouse Ensembl gene IDs.
## ------------------------------------------------------------

background_ensembl <- sub(
  "\\.[0-9]+$",
  "",
  rownames(vsd_mat)
)

background_mapping <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys = unique(background_ensembl),
  keytype = "ENSEMBL",
  columns = c("ENSEMBL", "ENTREZID", "SYMBOL")
) %>%
  filter(
    !is.na(ENTREZID),
    ENTREZID != ""
  ) %>%
  distinct(ENTREZID, .keep_all = TRUE)

background_entrez <- unique(
  as.character(background_mapping$ENTREZID)
)

message(
  "Background universe contains ",
  length(background_entrez),
  " unique Entrez IDs."
)


## ------------------------------------------------------------
## 3. GO enrichment function
## ------------------------------------------------------------

run_go <- function(
    deg_table,
    comparison_name,
    direction_name,
    lfc_direction = c("up", "down"),
    universe_entrez,
    output_directory,
    padj_cutoff = 0.05,
    lfc_cutoff = 1,
    minimum_mapped_genes = 5,
    ontology = "BP",
    simplify_cutoff = 0.7,
    show_categories = 15
) {
  
  lfc_direction <- match.arg(lfc_direction)
  
  required_columns <- c(
    "gene_id",
    "gene_symbol",
    "padj",
    "log2FoldChange"
  )
  
  missing_columns <- setdiff(
    required_columns,
    colnames(deg_table)
  )
  
  if (length(missing_columns) > 0) {
    stop(
      comparison_name,
      " ",
      direction_name,
      ": DEG table is missing columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }
  
  ## Clean Ensembl IDs
  deg_table_clean <- deg_table %>%
    mutate(
      gene_id_clean = sub(
        "\\.[0-9]+$",
        "",
        as.character(gene_id)
      )
    )
  
  ## Select genes by direction
  if (lfc_direction == "up") {
    
    selected_genes <- deg_table_clean %>%
      filter(
        !is.na(padj),
        is.finite(padj),
        padj < padj_cutoff,
        !is.na(log2FoldChange),
        is.finite(log2FoldChange),
        log2FoldChange >= lfc_cutoff
      )
    
  } else {
    
    selected_genes <- deg_table_clean %>%
      filter(
        !is.na(padj),
        is.finite(padj),
        padj < padj_cutoff,
        !is.na(log2FoldChange),
        is.finite(log2FoldChange),
        log2FoldChange <= -lfc_cutoff
      )
  }
  
  selected_genes <- selected_genes %>%
    arrange(padj, desc(abs(log2FoldChange))) %>%
    distinct(gene_id_clean, .keep_all = TRUE)
  
  message(
    comparison_name,
    " ",
    direction_name,
    ": ",
    nrow(selected_genes),
    " significant genes before ID mapping."
  )
  
  if (nrow(selected_genes) < minimum_mapped_genes) {
    message(
      comparison_name,
      " ",
      direction_name,
      ": fewer than ",
      minimum_mapped_genes,
      " significant genes; GO skipped."
    )
    
    return(NULL)
  }
  
  ## Map Ensembl IDs to Entrez IDs
  gene_mapping <- AnnotationDbi::select(
    org.Mm.eg.db,
    keys = unique(selected_genes$gene_id_clean),
    keytype = "ENSEMBL",
    columns = c(
      "ENSEMBL",
      "ENTREZID",
      "SYMBOL"
    )
  ) %>%
    filter(
      !is.na(ENTREZID),
      ENTREZID != ""
    ) %>%
    distinct(ENTREZID, .keep_all = TRUE)
  
  selected_entrez <- unique(
    as.character(gene_mapping$ENTREZID)
  )
  
  ## Keep only genes represented in the specified universe
  selected_entrez <- intersect(
    selected_entrez,
    universe_entrez
  )
  
  message(
    comparison_name,
    " ",
    direction_name,
    ": ",
    length(selected_entrez),
    " genes mapped to Entrez IDs and background universe."
  )
  
  if (length(selected_entrez) < minimum_mapped_genes) {
    message(
      comparison_name,
      " ",
      direction_name,
      ": fewer than ",
      minimum_mapped_genes,
      " mapped genes; GO skipped."
    )
    
    return(NULL)
  }
  
  ## Run GO enrichment
  ego <- clusterProfiler::enrichGO(
    gene = selected_entrez,
    universe = universe_entrez,
    OrgDb = org.Mm.eg.db,
    keyType = "ENTREZID",
    ont = ontology,
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.20,
    minGSSize = 10,
    maxGSSize = 500,
    readable = TRUE
  )
  
  ## Check whether enrichment returned terms
  if (
    is.null(ego) ||
    nrow(as.data.frame(ego)) == 0
  ) {
    message(
      comparison_name,
      " ",
      direction_name,
      ": no significant GO terms."
    )
    
    return(NULL)
  }
  
  ## Save complete, unsimplified enrichment table
  original_table <- as.data.frame(ego)
  
  original_file <- file.path(
    output_directory,
    paste0(
      comparison_name,
      "_",
      direction_name,
      "_GO_",
      ontology,
      "_all_terms.csv"
    )
  )
  
  write.csv(
    original_table,
    original_file,
    row.names = FALSE
  )
  
  ## Remove semantically redundant terms
  ego_simplified <- tryCatch(
    {
      clusterProfiler::simplify(
        ego,
        cutoff = simplify_cutoff,
        by = "p.adjust",
        select_fun = min,
        measure = "Wang"
      )
    },
    error = function(e) {
      warning(
        comparison_name,
        " ",
        direction_name,
        ": simplify() failed: ",
        conditionMessage(e),
        ". Original enrichment results will be used."
      )
      
      ego
    }
  )
  
  simplified_table <- as.data.frame(
    ego_simplified
  )
  
  if (nrow(simplified_table) == 0) {
    message(
      comparison_name,
      " ",
      direction_name,
      ": no GO terms remained after simplify()."
    )
    
    return(
      list(
        original = ego,
        simplified = ego_simplified,
        mapping = gene_mapping
      )
    )
  }
  
  ## Save simplified table
  simplified_file <- file.path(
    output_directory,
    paste0(
      comparison_name,
      "_",
      direction_name,
      "_GO_",
      ontology,
      "_simplified.csv"
    )
  )
  
  write.csv(
    simplified_table,
    simplified_file,
    row.names = FALSE
  )
  
  ## Create simplified GO dot plot
  plot_title <- paste(
    gsub("_", " ", comparison_name),
    direction_name,
    "GO Biological Process"
  )
  
  p_dot <- enrichplot::dotplot(
    ego_simplified,
    showCategory = min(
      show_categories,
      nrow(simplified_table)
    ),
    orderBy = "GeneRatio",
    font.size = 11
  ) +
    ggtitle(plot_title) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold"
      ),
      axis.text.y = element_text(
        size = 10
      ),
      panel.grid.minor = element_blank()
    )
  
  dotplot_file <- file.path(
    output_directory,
    paste0(
      comparison_name,
      "_",
      direction_name,
      "_GO_",
      ontology,
      "_simplified_dotplot.pdf"
    )
  )
  
  ggsave(
    filename = dotplot_file,
    plot = p_dot,
    width = 10,
    height = 7,
    device = cairo_pdf
  )
  
  ## Optional PNG for easy viewing
  png_file <- file.path(
    output_directory,
    paste0(
      comparison_name,
      "_",
      direction_name,
      "_GO_",
      ontology,
      "_simplified_dotplot.png"
    )
  )
  
  ggsave(
    filename = png_file,
    plot = p_dot,
    width = 10,
    height = 7,
    dpi = 300
  )
  
  message(
    comparison_name,
    " ",
    direction_name,
    ": ",
    nrow(original_table),
    " original terms; ",
    nrow(simplified_table),
    " terms after simplify()."
  )
  
  return(
    list(
      original = ego,
      simplified = ego_simplified,
      original_table = original_table,
      simplified_table = simplified_table,
      mapping = gene_mapping,
      plot = p_dot
    )
  )
}


## ------------------------------------------------------------
## 4. Run GO analysis
##
## SGE_vs_Control is excluded.
## ------------------------------------------------------------

comparisons_to_run <- setdiff(
  names(deg_list),
  "SGE_vs_Control"
)

go_results <- list()

for (comparison_name in comparisons_to_run) {
  
  message(
    "\n==============================\n",
    "Running: ",
    comparison_name,
    "\n=============================="
  )
  
  go_results[[paste0(
    comparison_name,
    "_Up"
  )]] <- run_go(
    deg_table = deg_list[[comparison_name]],
    comparison_name = comparison_name,
    direction_name = "Up",
    lfc_direction = "up",
    universe_entrez = background_entrez,
    output_directory = go_output_dir,
    padj_cutoff = 0.05,
    lfc_cutoff = 1,
    minimum_mapped_genes = 5,
    simplify_cutoff = 0.7,
    show_categories = 15
  )
  
  go_results[[paste0(
    comparison_name,
    "_Down"
  )]] <- run_go(
    deg_table = deg_list[[comparison_name]],
    comparison_name = comparison_name,
    direction_name = "Down",
    lfc_direction = "down",
    universe_entrez = background_entrez,
    output_directory = go_output_dir,
    padj_cutoff = 0.05,
    lfc_cutoff = 1,
    minimum_mapped_genes = 5,
    simplify_cutoff = 0.7,
    show_categories = 15
  )
}


## ------------------------------------------------------------
## 5. Summarize which analyses produced results
## ------------------------------------------------------------

go_summary <- data.frame(
  Analysis = names(go_results),
  Completed = !vapply(
    go_results,
    is.null,
    logical(1)
  )
)

print(go_summary)

write.csv(
  go_summary,
  file.path(
    go_output_dir,
    "GO_analysis_summary.csv"
  ),
  row.names = FALSE
)

##
for (nm in names(go_results)) {
  
  x <- go_results[[nm]]
  
  if (!is.null(x)) {
    cat(
      nm,
      ":",
      nrow(x$original_table),
      "original terms ->",
      nrow(x$simplified_table),
      "simplified terms\n"
    )
  }
}
