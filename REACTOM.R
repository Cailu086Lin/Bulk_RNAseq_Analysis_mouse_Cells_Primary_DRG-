library(dplyr)
library(ReactomePA)
library(clusterProfiler)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(enrichplot)
library(ggplot2)

run_reactome_debug <- function(
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
  
  dir.create(
    output_directory,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  prefix <- paste0(
    comparison_name,
    "_",
    direction_name
  )
  
  status_file <- file.path(
    output_directory,
    paste0(prefix, "_Reactome_status.txt")
  )
  
  cat(
    paste0(
      "Reactome analysis: ",
      comparison_name,
      " ",
      direction_name,
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
    "padj",
    "log2FoldChange"
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
      !is.na(log2FoldChange),
      is.finite(log2FoldChange)
    )
  
  if (lfc_direction == "up") {
    
    selected <- deg_clean %>%
      filter(
        padj < padj_cutoff,
        log2FoldChange >= lfc_cutoff
      )
    
  } else {
    
    selected <- deg_clean %>%
      filter(
        padj < padj_cutoff,
        log2FoldChange <= -lfc_cutoff
      )
  }
  
  selected <- selected %>%
    arrange(padj, desc(abs(log2FoldChange))) %>%
    distinct(gene_id_clean, .keep_all = TRUE)
  
  write_status(
    "Significant genes before mapping: ",
    nrow(selected)
  )
  
  write.csv(
    selected,
    file.path(
      output_directory,
      paste0(prefix, "_selected_genes.csv")
    ),
    row.names = FALSE
  )
  
  if (nrow(selected) < minimum_genes) {
    write_status(
      "SKIPPED: fewer than ",
      minimum_genes,
      " selected genes."
    )
    return(NULL)
  }
  
  mapping <- suppressMessages(
    AnnotationDbi::select(
      org.Mm.eg.db,
      keys = unique(selected$gene_id_clean),
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
    ) %>%
    distinct(ENTREZID, .keep_all = TRUE)
  
  write.csv(
    mapping,
    file.path(
      output_directory,
      paste0(prefix, "_Entrez_mapping.csv")
    ),
    row.names = FALSE
  )
  
  selected_entrez <- intersect(
    unique(as.character(mapping$ENTREZID)),
    unique(as.character(universe_entrez))
  )
  
  write_status(
    "Mapped Entrez IDs in background: ",
    length(selected_entrez)
  )
  
  write.table(
    selected_entrez,
    file.path(
      output_directory,
      paste0(prefix, "_Entrez_IDs.txt")
    ),
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
  
  if (length(selected_entrez) < minimum_genes) {
    write_status(
      "SKIPPED: fewer than ",
      minimum_genes,
      " mapped genes."
    )
    return(NULL)
  }
  
  reactome_all <- tryCatch(
    {
      ReactomePA::enrichPathway(
        gene = selected_entrez,
        universe = unique(as.character(universe_entrez)),
        organism = "mouse",
        pvalueCutoff = 1,
        qvalueCutoff = 1,
        pAdjustMethod = "BH",
        minGSSize = 10,
        maxGSSize = 500,
        readable = TRUE
      )
    },
    error = function(e) {
      write_status(
        "ERROR from enrichPathway(): ",
        conditionMessage(e)
      )
      NULL
    }
  )
  
  if (is.null(reactome_all)) {
    write_status(
      "FAILED: enrichPathway() returned NULL."
    )
    return(NULL)
  }
  
  all_table <- as.data.frame(
    reactome_all
  )
  
  write_status(
    "Reactome pathways returned: ",
    nrow(all_table)
  )
  
  write.csv(
    all_table,
    file.path(
      output_directory,
      paste0(prefix, "_Reactome_all_results.csv")
    ),
    row.names = FALSE
  )
  
  if (nrow(all_table) == 0) {
    write_status(
      "NO RESULT: Reactome returned no pathways."
    )
    return(NULL)
  }
  
  significant_table <- all_table %>%
    filter(
      !is.na(p.adjust),
      p.adjust < 0.05
    ) %>%
    arrange(p.adjust)
  
  write.csv(
    significant_table,
    file.path(
      output_directory,
      paste0(prefix, "_Reactome_significant.csv")
    ),
    row.names = FALSE
  )
  
  write_status(
    "Significant Reactome pathways: ",
    nrow(significant_table)
  )
  
  if (nrow(significant_table) == 0) {
    write_status(
      "No pathways passed adjusted P < 0.05."
    )
    
    return(
      list(
        enrichment_all = reactome_all,
        all_table = all_table,
        significant_table = significant_table,
        mapping = mapping
      )
    )
  }
  
  reactome_sig <- reactome_all
  
  exclude_terms <- c(
    "Signaling by GPCR",
    "GPCR ligand binding",
    "GPCR downstream signalling",
    "Class A/1 (Rhodopsin-like receptors)",
    "Peptide ligand-binding receptors",
    "G alpha (i) signalling events"
  )
  
  plot_df <- reactome_sig@result %>%
    dplyr::filter(
      !is.na(p.adjust),
      p.adjust < 0.05,
      !Description %in% exclude_terms
    ) %>%
    dplyr::mutate(
      # Convert values such as "12/76" into numeric GeneRatio
      GeneRatio_numeric = vapply(
        strsplit(as.character(GeneRatio), "/"),
        function(x) as.numeric(x[1]) / as.numeric(x[2]),
        numeric(1)
      )
    ) %>%
    dplyr::arrange(dplyr::desc(GeneRatio_numeric)) %>%
    dplyr::slice_head(n = 15) %>%
    dplyr::mutate(
      # First row from descending sort appears at the top
      Description = factor(
        Description,
        levels = rev(Description)
      )
    )
  
  print(
    plot_df %>%
      dplyr::select(
        Description,
        GeneRatio,
        GeneRatio_numeric,
        Count,
        p.adjust
      )
  )
  
  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = GeneRatio_numeric,
      y = Description,
      size = Count,
      color = p.adjust
    )
  ) +
    ggplot2::geom_point() +
    ggplot2::scale_color_gradient(
      low = "#E64B35",
      high = "#4DBBD5",
      trans = "reverse"
    ) +
    ggplot2::labs(
      title = paste(
        gsub("_", " ", comparison_name),
        direction_name,
        "Reactome pathways"
      ),
      x = "GeneRatio",
      y = NULL,
      size = "Count",
      color = "Adjusted P"
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        face = "bold"
      ),
      axis.text.y = ggplot2::element_text(size = 10),
      panel.grid.minor = ggplot2::element_blank()
    )
  
  pdf_file <- file.path(
    output_directory,
    paste0(prefix, "_Reactome_dotplot.pdf")
  )
  
  png_file <- file.path(
    output_directory,
    paste0(prefix, "_Reactome_dotplot.png")
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
    "SUCCESS: Reactome tables and plots saved."
  )
  
  invisible(
    list(
      enrichment = reactome_sig,
      all_table = all_table,
      significant_table = significant_table,
      mapping = mapping,
      plot = p
    )
  )
}


comparisons_to_run <- setdiff(
  names(deg_list),
  "SGE_vs_Control"
)

reactome_results <- list()

for (comparison_name in comparisons_to_run) {
  
  for (direction_name in c("Up", "Down")) {
    
    direction_code <- ifelse(
      direction_name == "Up",
      "up",
      "down"
    )
    
    result_name <- paste0(
      comparison_name,
      "_",
      direction_name
    )
    
    reactome_results[[result_name]] <- run_reactome_debug(
      deg_table = deg_list[[comparison_name]],
      comparison_name = comparison_name,
      direction_name = direction_name,
      lfc_direction = direction_code,
      universe_entrez = background_entrez,
      output_directory = reactome_output_dir,
      padj_cutoff = 0.05,
      lfc_cutoff = 1,
      minimum_genes = 5
    )
  }
}
