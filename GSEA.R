pathway_output_dir <- file.path(results_dir, "Pathway_Analysis")

kegg_output_dir <- file.path(pathway_output_dir, "KEGG")
reactome_output_dir <- file.path(pathway_output_dir, "Reactome")
gsea_output_dir <- file.path(pathway_output_dir, "GSEA")

dir.create(kegg_output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(reactome_output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(gsea_output_dir, recursive = TRUE, showWarnings = FALSE)

cat(
  "KEGG folder:     ", normalizePath(kegg_output_dir, mustWork = FALSE), "\n",
  "Reactome folder: ", normalizePath(reactome_output_dir, mustWork = FALSE), "\n",
  "GSEA folder:     ", normalizePath(gsea_output_dir, mustWork = FALSE), "\n"
)
#####
library(dplyr)
library(msigdbr)
library(fgsea)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(ggplot2)

hallmark_mouse_df <- msigdbr::msigdbr(
  species = "Mus musculus",
  collection = "H"
)

hallmark_mouse <- split(
  hallmark_mouse_df$ncbi_gene,
  hallmark_mouse_df$gs_name
)

hallmark_mouse <- lapply(
  hallmark_mouse,
  function(x) {
    unique(
      as.character(
        x[!is.na(x)]
      )
    )
  }
)

length(hallmark_mouse)

hallmark_mouse_df <- msigdbr::msigdbr(
  species = "Mus musculus",
  category = "H"
)
run_hallmark_gsea <- function(
    deg_table,
    comparison_name,
    pathways,
    output_directory,
    rank_column = "stat",
    padj_cutoff = 0.05,
    show_categories = 20,
    min_size = 10,
    max_size = 500
) {
  
  dir.create(
    output_directory,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  prefix <- comparison_name
  
  status_file <- file.path(
    output_directory,
    paste0(prefix, "_Hallmark_GSEA_status.txt")
  )
  
  cat(
    paste0(
      "Hallmark GSEA: ",
      comparison_name,
      "\n"
    ),
    file = status_file
  )
  
  write_status <- function(...) {
    text <- paste0(...)
    message(text)
    
    cat(
      text,
      "\n",
      file = status_file,
      append = TRUE
    )
  }
  
  required_columns <- c(
    "gene_id",
    rank_column
  )
  
  missing_columns <- setdiff(
    required_columns,
    colnames(deg_table)
  )
  
  if (length(missing_columns) > 0) {
    write_status(
      "ERROR: missing columns: ",
      paste(missing_columns, collapse = ", ")
    )
    return(NULL)
  }
  
  ranked_table <- deg_table %>%
    mutate(
      gene_id_clean = sub(
        "\\.[0-9]+$",
        "",
        as.character(gene_id)
      ),
      rank_value = .data[[rank_column]]
    ) %>%
    filter(
      !is.na(gene_id_clean),
      gene_id_clean != "",
      !is.na(rank_value),
      is.finite(rank_value)
    )
  
  write_status(
    "Genes with finite ranking values: ",
    nrow(ranked_table)
  )
  
  mapping <- suppressMessages(
    AnnotationDbi::select(
      org.Mm.eg.db,
      keys = unique(ranked_table$gene_id_clean),
      keytype = "ENSEMBL",
      columns = c(
        "ENSEMBL",
        "ENTREZID",
        "SYMBOL"
      )
    )
  ) %>%
    filter(
      !is.na(ENTREZID),
      ENTREZID != ""
    )
  
  ranked_mapped <- ranked_table %>%
    inner_join(
      mapping,
      by = c(
        "gene_id_clean" = "ENSEMBL"
      )
    ) %>%
    arrange(
      desc(abs(rank_value))
    ) %>%
    distinct(
      ENTREZID,
      .keep_all = TRUE
    )
  
  write_status(
    "Genes mapped to unique Entrez IDs: ",
    nrow(ranked_mapped)
  )
  
  write.csv(
    ranked_mapped,
    file.path(
      output_directory,
      paste0(prefix, "_Hallmark_GSEA_ranked_genes.csv")
    ),
    row.names = FALSE
  )
  
  ranks <- ranked_mapped$rank_value
  
  names(ranks) <- as.character(
    ranked_mapped$ENTREZID
  )
  
  ranks <- sort(
    ranks,
    decreasing = TRUE
  )
  
  if (length(ranks) < 100) {
    write_status(
      "SKIPPED: fewer than 100 ranked genes."
    )
    return(NULL)
  }
  
  gsea_result <- tryCatch(
    {
      fgsea::fgsea(
        pathways = pathways,
        stats = ranks,
        minSize = min_size,
        maxSize = max_size,
        eps = 0
      )
    },
    error = function(e) {
      write_status(
        "ERROR from fgsea(): ",
        conditionMessage(e)
      )
      NULL
    }
  )
  
  if (is.null(gsea_result)) {
    write_status(
      "FAILED: fgsea returned NULL."
    )
    return(NULL)
  }
  
  gsea_table <- as.data.frame(
    gsea_result
  ) %>%
    arrange(padj)
  
  if ("leadingEdge" %in% colnames(gsea_table)) {
    gsea_table$leadingEdge <- vapply(
      gsea_table$leadingEdge,
      paste,
      collapse = ";",
      FUN.VALUE = character(1)
    )
  }
  
  write_status(
    "Hallmark pathways tested: ",
    nrow(gsea_table)
  )
  
  write.csv(
    gsea_table,
    file.path(
      output_directory,
      paste0(prefix, "_Hallmark_GSEA_all_results.csv")
    ),
    row.names = FALSE
  )
  
  significant_table <- gsea_table %>%
    filter(
      !is.na(padj),
      padj < padj_cutoff
    ) %>%
    arrange(padj)
  
  write.csv(
    significant_table,
    file.path(
      output_directory,
      paste0(prefix, "_Hallmark_GSEA_significant.csv")
    ),
    row.names = FALSE
  )
  
  write_status(
    "Significant pathways at FDR < ",
    padj_cutoff,
    ": ",
    nrow(significant_table)
  )
  
  if (nrow(significant_table) == 0) {
    write_status(
      "No significant Hallmark pathways."
    )
    
    return(
      list(
        all_table = gsea_table,
        significant_table = significant_table,
        ranks = ranks
      )
    )
  }
  
  plot_table <- significant_table %>%
    filter(NES > 0) %>%
    arrange(padj, desc(NES)) %>%
    slice_head(n = show_categories) %>%
    mutate(
      pathway_label = gsub(
        "^HALLMARK_",
        "",
        pathway
      ),
      pathway_label = gsub(
        "_",
        " ",
        pathway_label
      ),
      pathway_label = factor(
        pathway_label,
        levels = rev(pathway_label)
      )
    )
  
  p <- ggplot(
    plot_table,
    aes(
      x = NES,
      y = pathway_label,
      size = size,
      fill = -log10(padj)
    )
  ) +
    geom_point(
      shape = 21
    ) +
    labs(
      title = paste(
        gsub("_", " ", comparison_name),
        "Positively Enriched Hallmark Pathways"
      ),
      x = "Normalized enrichment score",
      y = NULL,
      size = "Gene set size",
      fill = "-log10 FDR"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold"
      ),
      axis.text.y = element_text(size = 10),
      panel.grid.minor = element_blank()
    )
  
  pdf_file <- file.path(
    output_directory,
    paste0(prefix, "_Hallmark_GSEA_dotplot.pdf")
  )
  
  png_file <- file.path(
    output_directory,
    paste0(prefix, "_Hallmark_GSEA_dotplot.png")
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
  
  write_status(
    "SUCCESS: GSEA tables and plots saved."
  )
  
  invisible(
    list(
      all_table = gsea_table,
      significant_table = significant_table,
      ranks = ranks,
      plot = p
    )
  )
}

gsea_results <- list()

for (comparison_name in names(deg_list)) {
  
  gsea_results[[comparison_name]] <- run_hallmark_gsea(
    deg_table = deg_list[[comparison_name]],
    comparison_name = comparison_name,
    pathways = hallmark_mouse,
    output_directory = gsea_output_dir,
    rank_column = "stat",
    padj_cutoff = 0.05,
    show_categories = 20
  )
}

##
cat("\nKEGG files:\n")
print(
  list.files(
    kegg_output_dir,
    full.names = TRUE
  )
)

cat("\nReactome files:\n")
print(
  list.files(
    reactome_output_dir,
    full.names = TRUE
  )
)

cat("\nGSEA files:\n")
print(
  list.files(
    gsea_output_dir,
    full.names = TRUE
  )
)
