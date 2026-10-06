## ============================================================
## Shared ESP/SGE response analysis
## All 9 samples retained
## ============================================================
#1. Load packages and define settings
library(DESeq2)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(ComplexUpset)
library(eulerr)
library(clusterProfiler)
library(ReactomePA)
library(org.Mm.eg.db)
library(msigdbr)
library(fgsea)
library(data.table)

dir.create("Shared_ESP_SGE_analysis", showWarnings = FALSE)
dir.create("Shared_ESP_SGE_analysis/Tables", showWarnings = FALSE)
dir.create("Shared_ESP_SGE_analysis/Figures", showWarnings = FALSE)
dir.create("Shared_ESP_SGE_analysis/Pathways", showWarnings = FALSE)

## Primary thresholds
padj_cutoff <- 0.05
lfc_cutoff  <- 0.50

## More inclusive candidate threshold
nominal_p_cutoff <- 0.05
candidate_lfc_cutoff <- 0.25


###2. Confirm all nine samples and run DESeq2
## Check samples
colnames(dds)
colData(dds)

## Confirm all 9 samples are present
stopifnot(ncol(dds) == 9)

## Make Control the reference
dds$condition <- relevel(factor(dds$condition), ref = "Control")

## Remove unused factor levels
dds$condition <- droplevels(dds$condition)

table(dds$condition)

## Re-run DESeq2 after confirming the design
design(dds) <- ~ condition
dds <- DESeq(dds)

resultsNames(dds)

###3. Extract ESP and SGE results versus Control
res_ESP <- results(
  dds,
  contrast = c("condition", "ESP", "Control"),
  alpha = padj_cutoff
)

res_SGE <- results(
  dds,
  contrast = c("condition", "SGE", "Control"),
  alpha = padj_cutoff
)

res_ESP_df <- as.data.frame(res_ESP) %>%
  rownames_to_column("gene_id") %>%
  rename(
    baseMean_ESP = baseMean,
    log2FC_ESP   = log2FoldChange,
    lfcSE_ESP    = lfcSE,
    stat_ESP     = stat,
    pvalue_ESP   = pvalue,
    padj_ESP     = padj
  )

res_SGE_df <- as.data.frame(res_SGE) %>%
  rownames_to_column("gene_id") %>%
  rename(
    baseMean_SGE = baseMean,
    log2FC_SGE   = log2FoldChange,
    lfcSE_SGE    = lfcSE,
    stat_SGE     = stat,
    pvalue_SGE   = pvalue,
    padj_SGE     = padj
  )


#4. Add mouse gene symbols
shared_results <- full_join(
  res_ESP_df,
  res_SGE_df,
  by = "gene_id"
) %>%
  mutate(
    ensembl_id = sub("\\..*$", "", gene_id)
  )

gene_annotation <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys = unique(shared_results$ensembl_id),
  keytype = "ENSEMBL",
  columns = c("SYMBOL", "ENTREZID", "GENENAME")
) %>%
  distinct(ENSEMBL, .keep_all = TRUE)

shared_results <- shared_results %>%
  left_join(
    gene_annotation,
    by = c("ensembl_id" = "ENSEMBL")
  ) %>%
  mutate(
    gene_symbol = ifelse(
      is.na(SYMBOL) | SYMBOL == "",
      gene_id,
      SYMBOL
    )
  )

#5. Define strict shared upregulated and downregulated genes
shared_results <- shared_results %>%
  mutate(
    strict_ESP_up =
      !is.na(padj_ESP) &
      padj_ESP < padj_cutoff &
      log2FC_ESP >= lfc_cutoff,
    
    strict_SGE_up =
      !is.na(padj_SGE) &
      padj_SGE < padj_cutoff &
      log2FC_SGE >= lfc_cutoff,
    
    strict_ESP_down =
      !is.na(padj_ESP) &
      padj_ESP < padj_cutoff &
      log2FC_ESP <= -lfc_cutoff,
    
    strict_SGE_down =
      !is.na(padj_SGE) &
      padj_SGE < padj_cutoff &
      log2FC_SGE <= -lfc_cutoff,
    
    strict_shared_up = strict_ESP_up & strict_SGE_up,
    strict_shared_down = strict_ESP_down & strict_SGE_down,
    
    strict_shared_direction = case_when(
      strict_shared_up   ~ "Shared up",
      strict_shared_down ~ "Shared down",
      TRUE               ~ "Not strict shared"
    )
  )

strict_shared_up <- shared_results %>%
  filter(strict_shared_up) %>%
  arrange(pmax(padj_ESP, padj_SGE))

strict_shared_down <- shared_results %>%
  filter(strict_shared_down) %>%
  arrange(pmax(padj_ESP, padj_SGE))

strict_shared_all <- bind_rows(
  strict_shared_up,
  strict_shared_down
)

cat("Strict shared upregulated genes:",
    nrow(strict_shared_up), "\n")

cat("Strict shared downregulated genes:",
    nrow(strict_shared_down), "\n")

##Save the results:
write.csv(
strict_shared_up,
"Shared_ESP_SGE_analysis/Tables/strict_shared_upregulated.csv",
row.names = FALSE
)

write.csv(
  strict_shared_down,
  "Shared_ESP_SGE_analysis/Tables/strict_shared_downregulated.csv",
  row.names = FALSE
)

write.csv(
  strict_shared_all,
  "Shared_ESP_SGE_analysis/Tables/strict_shared_all.csv",
  row.names = FALSE
)

##6. Define concordant candidate genes
shared_results <- shared_results %>%
  mutate(
    same_direction =
      sign(log2FC_ESP) == sign(log2FC_SGE) &
      log2FC_ESP != 0 &
      log2FC_SGE != 0,
    
    significant_either =
      (!is.na(padj_ESP) & padj_ESP < padj_cutoff) |
      (!is.na(padj_SGE) & padj_SGE < padj_cutoff),
    
    nominal_both =
      !is.na(pvalue_ESP) &
      !is.na(pvalue_SGE) &
      pvalue_ESP < nominal_p_cutoff &
      pvalue_SGE < nominal_p_cutoff,
    
    adequate_effect_both =
      abs(log2FC_ESP) >= candidate_lfc_cutoff &
      abs(log2FC_SGE) >= candidate_lfc_cutoff,
    
    concordant_candidate =
      same_direction &
      significant_either &
      nominal_both &
      adequate_effect_both,
    
    concordant_direction = case_when(
      concordant_candidate &
        log2FC_ESP > 0 &
        log2FC_SGE > 0 ~ "Concordant up",
      
      concordant_candidate &
        log2FC_ESP < 0 &
        log2FC_SGE < 0 ~ "Concordant down",
      
      TRUE ~ "Other"
    ),
    
    conservative_padj = pmax(
      replace_na(padj_ESP, 1),
      replace_na(padj_SGE, 1)
    ),
    
    mean_log2FC = (log2FC_ESP + log2FC_SGE) / 2,
    
    minimum_abs_log2FC = pmin(
      abs(log2FC_ESP),
      abs(log2FC_SGE)
    )
  )

concordant_candidates <- shared_results %>%
  filter(concordant_candidate) %>%
  arrange(
    conservative_padj,
    desc(minimum_abs_log2FC)
  )

write.csv(
  concordant_candidates,
  "Shared_ESP_SGE_analysis/Tables/concordant_candidate_genes.csv",
  row.names = FALSE
)

table(concordant_candidates$concordant_direction)

##Figures
esp_up_genes <- shared_results %>%
  filter(strict_ESP_up) %>%
  pull(gene_symbol) %>%
  unique()

sge_up_genes <- shared_results %>%
  filter(strict_SGE_up) %>%
  pull(gene_symbol) %>%
  unique()

esp_down_genes <- shared_results %>%
  filter(strict_ESP_down) %>%
  pull(gene_symbol) %>%
  unique()

sge_down_genes <- shared_results %>%
  filter(strict_SGE_down) %>%
  pull(gene_symbol) %>%
  unique()

## Upregulated genes
venn_up <- eulerr::euler(
  list(
    "ESP vs Control up" = esp_up_genes,
    "SGE vs Control up" = sge_up_genes
  )
)

pdf(
  "Shared_ESP_SGE_analysis/Figures/Venn_shared_upregulated.pdf",
  width = 7,
  height = 6
)

plot(
  venn_up,
  quantities = TRUE,
  labels = TRUE,
  edges = TRUE
)

dev.off()


##Downregulated genes
venn_down <- eulerr::euler(
  list(
    "ESP vs Control down" = esp_down_genes,
    "SGE vs Control down" = sge_down_genes
  )
)

pdf(
  "Shared_ESP_SGE_analysis/Figures/Venn_shared_downregulated.pdf",
  width = 7,
  height = 6
)

plot(
  venn_down,
  quantities = TRUE,
  labels = TRUE,
  edges = TRUE
)

dev.off()

##8. Log2-fold-change scatter plot
scatter_df <- shared_results %>%
  filter(
    !is.na(log2FC_ESP),
    !is.na(log2FC_SGE)
  ) %>%
  mutate(
    plot_group = case_when(
      strict_shared_up   ~ "Strict shared up",
      strict_shared_down ~ "Strict shared down",
      
      concordant_direction == "Concordant up" ~
        "Candidate shared up",
      
      concordant_direction == "Concordant down" ~
        "Candidate shared down",
      
      log2FC_ESP > 0 & log2FC_SGE > 0 ~
        "Same direction",
      
      log2FC_ESP < 0 & log2FC_SGE < 0 ~
        "Same direction",
      
      TRUE ~ "Discordant or unchanged"
    )
  )

##Choose labels:
label_genes <- scatter_df %>%
  filter(
    strict_shared_up |
      strict_shared_down |
      concordant_candidate
  ) %>%
  arrange(
    conservative_padj,
    desc(minimum_abs_log2FC)
  ) %>%
  slice_head(n = 25)

##Create plot
p_scatter <- ggplot(
  scatter_df,
  aes(
    x = log2FC_SGE,
    y = log2FC_ESP
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.4
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.4
  ) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dotted",
    linewidth = 0.5
  ) +
  geom_point(
    aes(shape = plot_group),
    alpha = 0.55,
    size = 1.7
  ) +
  geom_point(
    data = scatter_df %>%
      filter(strict_shared_up | strict_shared_down),
    size = 2.7
  ) +
  ggrepel::geom_text_repel(
    data = label_genes,
    aes(label = gene_symbol),
    size = 3,
    max.overlaps = Inf,
    box.padding = 0.4,
    point.padding = 0.25,
    min.segment.length = 0
  ) +
  coord_cartesian(
    xlim = quantile(
      scatter_df$log2FC_SGE,
      c(0.005, 0.995),
      na.rm = TRUE
    ),
    ylim = quantile(
      scatter_df$log2FC_ESP,
      c(0.005, 0.995),
      na.rm = TRUE
    )
  ) +
  labs(
    title = "Shared transcriptional responses to ESP and SGE",
    subtitle = paste0(
      "Strict shared: adjusted P < ", padj_cutoff,
      " and |log2FC| ≥ ", lfc_cutoff,
      " in both comparisons"
    ),
    x = "SGE vs Control log2 fold change",
    y = "ESP vs Control log2 fold change",
    shape = "Gene category"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold")
  )

ggsave(
  "Shared_ESP_SGE_analysis/Figures/shared_log2FC_scatter.pdf",
  p_scatter,
  width = 8,
  height = 7
)

ggsave(
  "Shared_ESP_SGE_analysis/Figures/shared_log2FC_scatter.png",
  p_scatter,
  width = 8,
  height = 7,
  dpi = 400
)

##Calculate the overall fold-change correlation:
fc_cor <- cor.test(
  scatter_df$log2FC_SGE,
  scatter_df$log2FC_ESP,
  method = "spearman",
  use = "complete.obs"
)

capture.output(
  fc_cor,
  file = "Shared_ESP_SGE_analysis/Tables/log2FC_correlation.txt"
)

fc_cor

##9. Heatmap of shared genes
heatmap_genes <- unique(strict_shared_all$gene_id)

if (length(heatmap_genes) < 5) {
  message(
    "Fewer than 5 strict shared genes; using concordant candidates."
  )
  
  heatmap_genes <- unique(concordant_candidates$gene_id)
}

length(heatmap_genes)

#Transform counts:
vsd <- vst(dds, blind = FALSE)
vsd_mat <- assay(vsd)

#Match genes:
heatmap_genes_present <- intersect(
  heatmap_genes,
  rownames(vsd_mat)
)

heatmap_mat <- vsd_mat[
  heatmap_genes_present,
  ,
  drop = FALSE
]

#Add symbols and remove duplicate symbols:
heatmap_annotation <- shared_results %>%
  filter(gene_id %in% heatmap_genes_present) %>%
  select(gene_id, gene_symbol) %>%
  distinct(gene_id, .keep_all = TRUE)

heatmap_symbols <- heatmap_annotation$gene_symbol[
  match(rownames(heatmap_mat), heatmap_annotation$gene_id)
]

heatmap_symbols[
  is.na(heatmap_symbols) |
    heatmap_symbols == ""
] <- rownames(heatmap_mat)[
  is.na(heatmap_symbols) |
    heatmap_symbols == ""
]

rownames(heatmap_mat) <- make.unique(heatmap_symbols)

#Remove zero-variance rows and scale:
row_variance <- apply(heatmap_mat, 1, var)

heatmap_mat <- heatmap_mat[
  is.finite(row_variance) &
    row_variance > 0,
  ,
  drop = FALSE
]

heatmap_z <- t(scale(t(heatmap_mat)))

heatmap_z[!is.finite(heatmap_z)] <- 0

#Sample annotation:
sample_annotation <- as.data.frame(colData(dds)) %>%
  select(condition)

colnames(sample_annotation) <- "Condition"
rownames(sample_annotation) <- colnames(dds)

#Plot
pheatmap::pheatmap(
  heatmap_z,
  annotation_col = sample_annotation,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 7,
  fontsize_col = 10,
  border_color = NA,
  main = "Shared ESP and SGE transcriptional signature",
  filename =
    "Shared_ESP_SGE_analysis/Figures/shared_gene_heatmap.pdf",
  width = 9,
  height = max(6, 3 + 0.16 * nrow(heatmap_z))
)


#A second heatmap using group-average expression can make the pattern clearer:
condition <- factor(
  colData(dds)$condition,
  levels = c("Control", "SGE", "ESP")
)

mean_heatmap <- sapply(
  levels(condition),
  function(group_name) {
    rowMeans(
      heatmap_z[
        ,
        condition == group_name,
        drop = FALSE
      ]
    )
  }
)

colnames(mean_heatmap) <- levels(condition)

pheatmap::pheatmap(
  mean_heatmap,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  border_color = NA,
  main = "Mean shared signature by treatment",
  filename =
    "Shared_ESP_SGE_analysis/Figures/shared_gene_group_means.pdf",
  width = 6,
  height = max(6, 3 + 0.16 * nrow(mean_heatmap))
)

##PATHWAY ANALYSES
#10. Prepare shared upregulated and downregulated Entrez IDs
pathway_up_df <- strict_shared_up
pathway_down_df <- strict_shared_down

if (nrow(pathway_up_df) < 5) {
  pathway_up_df <- concordant_candidates %>%
    filter(concordant_direction == "Concordant up")
  
  message(
    "Using concordant-up candidates for exploratory enrichment."
  )
}

if (nrow(pathway_down_df) < 5) {
  pathway_down_df <- concordant_candidates %>%
    filter(concordant_direction == "Concordant down")
  
  message(
    "Using concordant-down candidates for exploratory enrichment."
  )
}

shared_up_entrez <- pathway_up_df$ENTREZID %>%
  na.omit() %>%
  unique()

shared_down_entrez <- pathway_down_df$ENTREZID %>%
  na.omit() %>%
  unique()

#Define the tested-gene universe:
background_entrez <- shared_results %>%
  filter(
    !is.na(pvalue_ESP) |
      !is.na(pvalue_SGE)
  ) %>%
  pull(ENTREZID) %>%
  na.omit() %>%
  unique()

##11. GO Biological Process enrichment
run_go_shared <- function(
    entrez_genes,
    direction,
    universe_genes
) {
  
  if (length(entrez_genes) < 5) {
    message(
      "Skipping GO ", direction,
      ": fewer than 5 mapped genes."
    )
    return(NULL)
  }
  
  ego <- enrichGO(
    gene = entrez_genes,
    universe = universe_genes,
    OrgDb = org.Mm.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.20,
    readable = TRUE
  )
  
  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
    message("No GO terms for ", direction)
    return(NULL)
  }
  
  ego_simple <- clusterProfiler::simplify(
    ego,
    cutoff = 0.70,
    by = "p.adjust",
    select_fun = min
  )
  
  write.csv(
    as.data.frame(ego_simple),
    paste0(
      "Shared_ESP_SGE_analysis/Pathways/GO_",
      direction,
      ".csv"
    ),
    row.names = FALSE
  )
  
  pdf(
    paste0(
      "Shared_ESP_SGE_analysis/Pathways/GO_",
      direction,
      "_dotplot.pdf"
    ),
    width = 9,
    height = 7
  )
  
  print(
    dotplot(
      ego_simple,
      showCategory = 15
    ) +
      ggtitle(
        paste(
          "Shared",
          direction,
          "genes: GO Biological Process"
        )
      )
  )
  
  dev.off()
  
  ego_simple
}

go_shared_up <- run_go_shared(
  shared_up_entrez,
  "upregulated",
  background_entrez
)

go_shared_down <- run_go_shared(
  shared_down_entrez,
  "downregulated",
  background_entrez
)

##12. Reactome enrichment
run_reactome_shared <- function(
    entrez_genes,
    direction,
    universe_genes
) {
  
  if (length(entrez_genes) < 5) {
    message(
      "Skipping Reactome ", direction,
      ": fewer than 5 mapped genes."
    )
    return(NULL)
  }
  
  reactome_result <- ReactomePA::enrichPathway(
    gene = entrez_genes,
    universe = universe_genes,
    organism = "mouse",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    qvalueCutoff = 0.20,
    readable = TRUE
  )
  
  if (
    is.null(reactome_result) ||
    nrow(as.data.frame(reactome_result)) == 0
  ) {
    message("No Reactome terms for ", direction)
    return(NULL)
  }
  
  write.csv(
    as.data.frame(reactome_result),
    paste0(
      "Shared_ESP_SGE_analysis/Pathways/Reactome_",
      direction,
      ".csv"
    ),
    row.names = FALSE
  )
  
  pdf(
    paste0(
      "Shared_ESP_SGE_analysis/Pathways/Reactome_",
      direction,
      "_dotplot.pdf"
    ),
    width = 9,
    height = 7
  )
  
  print(
    dotplot(
      reactome_result,
      showCategory = 15
    ) +
      ggtitle(
        paste(
          "Shared",
          direction,
          "genes: Reactome"
        )
      )
  )
  
  dev.off()
  
  reactome_result
}

reactome_shared_up <- run_reactome_shared(
  shared_up_entrez,
  "upregulated",
  background_entrez
)

reactome_shared_down <- run_reactome_shared(
  shared_down_entrez,
  "downregulated",
  background_entrez
)

##13. Hallmark over-representation analysis
hallmark_mouse <- msigdbr(
  species = "Mus musculus",
  category = "H"
) %>%
  select(gs_name, entrez_gene) %>%
  distinct() %>%
  filter(!is.na(entrez_gene))

hallmark_term2gene <- hallmark_mouse %>%
  select(
    term = gs_name,
    gene = entrez_gene
  )

#Run enrichment:
run_hallmark_ora <- function(
    entrez_genes,
    direction,
    universe_genes
) {
  
  if (length(entrez_genes) < 5) {
    message(
      "Skipping Hallmark ", direction,
      ": fewer than 5 mapped genes."
    )
    return(NULL)
  }
  
  hallmark_result <- enricher(
    gene = entrez_genes,
    universe = universe_genes,
    TERM2GENE = hallmark_term2gene,
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.20
  )
  
  if (
    is.null(hallmark_result) ||
    nrow(as.data.frame(hallmark_result)) == 0
  ) {
    message("No Hallmark terms for ", direction)
    return(NULL)
  }
  
  write.csv(
    as.data.frame(hallmark_result),
    paste0(
      "Shared_ESP_SGE_analysis/Pathways/Hallmark_ORA_",
      direction,
      ".csv"
    ),
    row.names = FALSE
  )
  
  pdf(
    paste0(
      "Shared_ESP_SGE_analysis/Pathways/Hallmark_ORA_",
      direction,
      "_dotplot.pdf"
    ),
    width = 9,
    height = 6
  )
  
  print(
    dotplot(
      hallmark_result,
      showCategory = 15
    ) +
      ggtitle(
        paste(
          "Shared",
          direction,
          "genes: Hallmark"
        )
      )
  )
  
  dev.off()
  
  hallmark_result
}

hallmark_shared_up <- run_hallmark_ora(
  shared_up_entrez,
  "upregulated",
  background_entrez
)

hallmark_shared_down <- run_hallmark_ora(
  shared_down_entrez,
  "downregulated",
  background_entrez
)

##14. Shared-direction Hallmark GSEA
shared_rank_df <- shared_results %>%
  filter(
    !is.na(stat_ESP),
    !is.na(stat_SGE),
    !is.na(ENTREZID)
  ) %>%
  mutate(
    concordant_sign = case_when(
      stat_ESP > 0 & stat_SGE > 0 ~ 1,
      stat_ESP < 0 & stat_SGE < 0 ~ -1,
      TRUE ~ 0
    ),
    
    shared_rank_score =
      concordant_sign *
      pmin(abs(stat_ESP), abs(stat_SGE))
  ) %>%
  group_by(ENTREZID) %>%
  slice_max(
    order_by = abs(shared_rank_score),
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup()

##Remove zero and tied scores:
shared_rank_df <- shared_rank_df %>%
  filter(shared_rank_score != 0) %>%
  arrange(desc(shared_rank_score))

shared_ranks <- shared_rank_df$shared_rank_score
names(shared_ranks) <- shared_rank_df$ENTREZID

## Tiny deterministic adjustment to avoid exact ties
shared_ranks <- shared_ranks +
  seq_along(shared_ranks) * 1e-10

shared_ranks <- sort(shared_ranks, decreasing = TRUE)

#Create pathway list:
hallmark_pathways <- split(
  hallmark_mouse$entrez_gene,
  hallmark_mouse$gs_name
)

hallmark_pathways <- lapply(
  hallmark_pathways,
  unique
)


#Run FGSEA:
set.seed(123)

fgsea_shared <- fgsea::fgsea(
  pathways = hallmark_pathways,
  stats = shared_ranks,
  minSize = 10,
  maxSize = 500,
  eps = 0
) %>%
  as.data.frame() %>%
  arrange(padj)


write.csv(
  fgsea_shared,
  "Shared_ESP_SGE_analysis/Pathways/Hallmark_shared_direction_GSEA.csv",
  row.names = FALSE
)

#Plot significant pathways:
fgsea_plot_df <- fgsea_shared %>%
  filter(!is.na(padj), padj < 0.05) %>%
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
    )
  ) %>%
  group_by(sign = ifelse(NES > 0, "Shared up", "Shared down")) %>%
  slice_max(
    order_by = abs(NES),
    n = 12,
    with_ties = FALSE
  ) %>%
  ungroup()

if (nrow(fgsea_plot_df) > 0) {
  
  fgsea_plot_df <- fgsea_plot_df %>%
    arrange(NES) %>%
    mutate(
      pathway_label = factor(
        pathway_label,
        levels = pathway_label
      )
    )
  
  p_fgsea <- ggplot(
    fgsea_plot_df,
    aes(
      x = NES,
      y = pathway_label,
      size = -log10(padj)
    )
  ) +
    geom_point() +
    geom_vline(
      xintercept = 0,
      linetype = "dashed"
    ) +
    labs(
      title = "Hallmark pathways shared by ESP and SGE",
      subtitle =
        "Conservative rank based on the weaker concordant Wald statistic",
      x = "Normalized enrichment score",
      y = NULL,
      size = expression(-log[10]("adjusted P"))
    ) +
    theme_classic(base_size = 12)
  
  ggsave(
    "Shared_ESP_SGE_analysis/Pathways/Hallmark_shared_direction_GSEA.pdf",
    p_fgsea,
    width = 9,
    height = 7
  )
}

#Candidate gene prioritization
candidate_table <- shared_results %>%
  filter(
    strict_shared_up |
      strict_shared_down |
      concordant_candidate
  ) %>%
  mutate(
    search_text = toupper(
      paste(
        gene_symbol,
        GENENAME
      )
    ),
    
    candidate_category = case_when(
      grepl(
        "TRPV|TRPA|TRPM|TRPC|TRANSIENT RECEPTOR POTENTIAL",
        search_text
      ) ~ "TRP channels",
      
      grepl(
        "^SCN|SODIUM CHANNEL",
        search_text
      ) ~ "Voltage-gated sodium channels",
      
      grepl(
        "^CACNA|^CACNB|^CACNG|CALCIUM CHANNEL",
        search_text
      ) ~ "Voltage-gated calcium channels",
      
      grepl(
        "^KCNA|^KCNB|^KCNC|^KCND|^KCNH|^KCNJ|^KCNK|^KCNQ|POTASSIUM CHANNEL",
        search_text
      ) ~ "Potassium channels",
      
      grepl(
        "^P2RX|^P2RY|PURINERGIC",
        search_text
      ) ~ "Purinergic signaling",
      
      grepl(
        "^PTGS|^PTGES|^PTGER|^PLA2|PROSTAGLANDIN",
        search_text
      ) ~ "Prostaglandin pathway",
      
      grepl(
        "^CALM|^CAMK|^ITPR|^RYR|^STIM|^ORAI|CALCIUM SIGNAL",
        search_text
      ) ~ "Calcium signaling",
      
      grepl(
        "^NPY|^TAC1|^CALCA|^CALCB|^SST|^PENK|^PDYN|NEUROPEPTIDE",
        search_text
      ) ~ "Neuropeptides",
      
      grepl(
        "^GPR|G PROTEIN.COUpled RECEPTOR|GPCR",
        search_text
      ) ~ "GPCR signaling",
      
      grepl(
        "^CXCL|^CCL|^CXCR|^CCR|CHEMOKINE",
        search_text
      ) ~ "Chemokines and receptors",
      
      grepl(
        "^IL[0-9]|^IL[0-9].*R|^TNF|^IFN|CYTOKINE",
        search_text
      ) ~ "Cytokines and receptors",
      
      grepl(
        "^GABR|^GLRA|^GLRB|^SLC6A1|^SLC6A5|GABA|GLYCINE RECEPTOR",
        search_text
      ) ~ "Inhibitory neurotransmission",
      
      grepl(
        "NEURON|SYNAP|AXON|NOCICEP|PAIN",
        search_text
      ) ~ "Other neuronal function",
      
      TRUE ~ "Other shared response"
    ),
    
    evidence_class = case_when(
      strict_shared_up ~ "Strict shared up",
      strict_shared_down ~ "Strict shared down",
      concordant_direction == "Concordant up" ~
        "Candidate shared up",
      concordant_direction == "Concordant down" ~
        "Candidate shared down",
      TRUE ~ "Other"
    ),
    
    weakest_abs_effect = pmin(
      abs(log2FC_ESP),
      abs(log2FC_SGE)
    ),
    
    weakest_pvalue = pmax(
      replace_na(pvalue_ESP, 1),
      replace_na(pvalue_SGE, 1)
    )
  ) %>%
  select(
    gene_id,
    gene_symbol,
    GENENAME,
    ENTREZID,
    candidate_category,
    evidence_class,
    log2FC_ESP,
    padj_ESP,
    log2FC_SGE,
    padj_SGE,
    weakest_abs_effect,
    weakest_pvalue
  ) %>%
  arrange(
    candidate_category == "Other shared response",
    weakest_pvalue,
    desc(weakest_abs_effect)
  )

##Save all candidates and the neuron-focused subset:
write.csv(
  candidate_table,
  "Shared_ESP_SGE_analysis/Tables/shared_candidate_gene_table_all.csv",
  row.names = FALSE
)

neuron_focused_candidates <- candidate_table %>%
  filter(candidate_category != "Other shared response")

write.csv(
  neuron_focused_candidates,
  "Shared_ESP_SGE_analysis/Tables/shared_candidate_gene_table_neuronal.csv",
  row.names = FALSE
)

neuron_focused_candidates

##16. Create a focused candidate plot
candidate_plot_df <- neuron_focused_candidates %>%
  slice_head(n = 30) %>%
  select(
    gene_symbol,
    candidate_category,
    log2FC_ESP,
    log2FC_SGE
  ) %>%
  pivot_longer(
    cols = c(log2FC_ESP, log2FC_SGE),
    names_to = "Comparison",
    values_to = "log2FC"
  ) %>%
  mutate(
    Comparison = recode(
      Comparison,
      log2FC_ESP = "ESP vs Control",
      log2FC_SGE = "SGE vs Control"
    )
  )

if (nrow(candidate_plot_df) > 0) {
  
  gene_order <- candidate_plot_df %>%
    group_by(gene_symbol) %>%
    summarise(
      mean_effect = mean(log2FC),
      .groups = "drop"
    ) %>%
    arrange(mean_effect) %>%
    pull(gene_symbol)
  
  candidate_plot_df$gene_symbol <- factor(
    candidate_plot_df$gene_symbol,
    levels = gene_order
  )
  
  p_candidates <- ggplot(
    candidate_plot_df,
    aes(
      x = log2FC,
      y = gene_symbol,
      shape = Comparison
    )
  ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed"
    ) +
    geom_point(
      size = 2.8,
      position = position_dodge(width = 0.5)
    ) +
    facet_wrap(
      ~ candidate_category,
      scales = "free_y",
      ncol = 2
    ) +
    labs(
      title = "Shared candidates related to neuronal signaling",
      x = "Log2 fold change versus Control",
      y = NULL
    ) +
    theme_classic(base_size = 11)
  
  ggsave(
    "Shared_ESP_SGE_analysis/Figures/shared_neuronal_candidates.pdf",
    p_candidates,
    width = 11,
    height = 9
  )
}

##17. Final summary table
analysis_summary <- data.frame(
  Gene_set = c(
    "ESP significant up",
    "SGE significant up",
    "Strict shared up",
    "ESP significant down",
    "SGE significant down",
    "Strict shared down",
    "Concordant candidate up",
    "Concordant candidate down"
  ),
  Number_of_genes = c(
    length(esp_up_genes),
    length(sge_up_genes),
    nrow(strict_shared_up),
    length(esp_down_genes),
    length(sge_down_genes),
    nrow(strict_shared_down),
    sum(
      concordant_candidates$concordant_direction ==
        "Concordant up"
    ),
    sum(
      concordant_candidates$concordant_direction ==
        "Concordant down"
    )
  )
)

write.csv(
  analysis_summary,
  "Shared_ESP_SGE_analysis/Tables/shared_analysis_summary.csv",
  row.names = FALSE
)

analysis_summary


