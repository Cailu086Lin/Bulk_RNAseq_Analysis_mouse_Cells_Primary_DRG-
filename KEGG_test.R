library(dplyr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(enrichplot)
library(ggplot2)

kegg_output_dir <- file.path(results_dir, "KEGG")

dir.create(
  kegg_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

run_kegg <- function(
    deg_table,
    comparison_name,
    direction_name,
    lfc_direction = c("up", "down"),
    universe_entrez,
    output_directory,
    padj_cutoff = 0.05,
    lfc_cutoff = 1,
    minimum_genes = 5,
    show_categories = 15
) {
  
  lfc_direction <- match.arg(lfc_direction)
  
  message(
    "\nRunning KEGG: ",
    comparison_name,
    " ",
    direction_name
  )
  
  required_columns <- c(
    "gene_id",
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
      ": missing columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }
  
  deg_clean <- deg_table %>%
    mutate(
      gene_id_clean = sub(
        "\\.[0-9]+$",
        "",
        as.character(gene_id)
      )
    ) %>%
    filter(
      !is.na(padj),
      is.finite(padj),
      padj < padj_cutoff,
      !is.na(log2FoldChange),
      is.finite(log2FoldChange)
    )
  
  if (lfc_direction == "up") {
    
    selected <- deg_clean %>%
      filter(log2FoldChange >= lfc_cutoff)
    
  } else {
    
    selected <- deg_clean %>%
      filter(log2FoldChange <= -lfc_cutoff)
  }
  
  selected <- selected %>%
    arrange(padj) %>%
    distinct(gene_id_clean, .keep_all = TRUE)
  
  message(
    comparison_name,
    " ",
    direction_name,
    ": ",
    nrow(selected),
    " significant genes before mapping."
  )
  
  if (nrow(selected) < minimum_genes) {
    message(
      comparison_name,
      " ",
      direction_name,
      ": fewer than ",
      minimum_genes,
      " significant genes; KEGG skipped."
    )
    return(NULL)
  }
  
  mapping <- AnnotationDbi::select(
    org.Mm.eg.db,
    keys = unique(selected$gene_id_clean),
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
  
  selected_entrez <- intersect(
    unique(as.character(mapping$ENTREZID)),
    as.character(universe_entrez)
  )
  
  message(
    comparison_name,
    " ",
    direction_name,
    ": ",
    length(selected_entrez),
    " mapped Entrez IDs."
  )
  
  if (length(selected_entrez) < minimum_genes) {
    message(
      comparison_name,
      " ",
      direction_name,
      ": fewer than ",
      minimum_genes,
      " mapped genes; KEGG skipped."
    )
    return(NULL)
  }
  
  ekegg <- tryCatch(
    {
      clusterProfiler::enrichKEGG(
        gene = selected_entrez,
        universe = as.character(universe_entrez),
        organism = "mmu",
        keyType = "kegg",
        pvalueCutoff = 0.05,
        pAdjustMethod = "BH",
        qvalueCutoff = 0.20,
        minGSSize = 10,
        maxGSSize = 500
      )
    },
    error = function(e) {
      message(
        comparison_name,
        " ",
        direction_name,
        ": enrichKEGG failed: ",
        conditionMessage(e)
      )
      return(NULL)
    }
  )
  
  if (is.null(ekegg)) {
    return(NULL)
  }
  
  kegg_table <- as.data.frame(ekegg)
  
  message(
    comparison_name,
    " ",
    direction_name,
    ": ",
    nrow(kegg_table),
    " significant KEGG pathways."
  )
  
  if (nrow(kegg_table) == 0) {
    message(
      comparison_name,
      " ",
      direction_name,
      ": no significant KEGG pathways at current cutoffs."
    )
    return(NULL)
  }
  
  ekegg_readable <- tryCatch(
    {
      clusterProfiler::setReadable(
        ekegg,
        OrgDb = org.Mm.eg.db,
        keyType = "ENTREZID"
      )
    },
    error = function(e) {
      message(
        "setReadable failed; using original KEGG result: ",
        conditionMessage(e)
      )
      ekegg
    }
  )
  
  kegg_table_readable <- as.data.frame(ekegg_readable)
  
  csv_file <- file.path(
    output_directory,
    paste0(
      comparison_name,
      "_",
      direction_name,
      "_KEGG.csv"
    )
  )
  
  write.csv(
    kegg_table_readable,
    csv_file,
    row.names = FALSE
  )
  
  p <- enrichplot::dotplot(
    ekegg_readable,
    showCategory = min(
      show_categories,
      nrow(kegg_table_readable)
    ),
    orderBy = "GeneRatio",
    font.size = 11
  ) +
    ggtitle(
      paste(
        gsub("_", " ", comparison_name),
        direction_name,
        "KEGG pathways"
      )
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold"
      ),
      panel.grid.minor = element_blank()
    )
  
  pdf_file <- file.path(
    output_directory,
    paste0(
      comparison_name,
      "_",
      direction_name,
      "_KEGG_dotplot.pdf"
    )
  )
  
  png_file <- file.path(
    output_directory,
    paste0(
      comparison_name,
      "_",
      direction_name,
      "_KEGG_dotplot.png"
    )
  )
  
  ggsave(
    filename = pdf_file,
    plot = p,
    width = 9,
    height = 7
  )
  
  ggsave(
    filename = png_file,
    plot = p,
    width = 9,
    height = 7,
    dpi = 300
  )
  
  message(
    "Saved:\n",
    csv_file,
    "\n",
    pdf_file,
    "\n",
    png_file
  )
  
  invisible(
    list(
      enrichment = ekegg_readable,
      table = kegg_table_readable,
      mapping = mapping,
      plot = p
    )
  )
}