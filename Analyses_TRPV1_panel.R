#1. Load packages and DEG tables
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(ggrepel)
library(forcats)
library(pheatmap)

outdir <- "TRPV1_regulatory_analysis"

dir.create(outdir, showWarnings = FALSE)
dir.create(file.path(outdir, "Tables"), showWarnings = FALSE)
dir.create(file.path(outdir, "Figures"), showWarnings = FALSE)

esp <- read.csv(
  "/mnt/easystore/DBH_26/Stephenson/downstream_GRCm39_decoy/Stephenson_results/DEG/ESP_vs_Control_all_genes.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

sge <- read.csv(
  "/mnt/easystore/DBH_26/Stephenson/downstream_GRCm39_decoy/Stephenson_results/DEG/SGE_vs_Control_all_genes.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

names(esp)
names(sge)

##Your files should contain these key columns:
c(
  "gene_id",
  "log2FoldChange",
  "pvalue",
  "padj",
  "gene_symbol"
)

#2. Define the curated TRPV1 regulatory panel
trpv1_panel <- tribble(
  ~gene_symbol, ~category,
  
  ## TRP channels
  "Trpv1", "TRP channels",
  "Trpv2", "TRP channels",
  "Trpv3", "TRP channels",
  "Trpv4", "TRP channels",
  "Trpa1", "TRP channels",
  "Trpm3", "TRP channels",
  "Trpm8", "TRP channels",
  
  ## Purinergic signaling
  "P2rx2", "Purinergic signaling",
  "P2rx3", "Purinergic signaling",
  "P2rx4", "Purinergic signaling",
  "P2rx7", "Purinergic signaling",
  "P2ry1", "Purinergic signaling",
  "P2ry2", "Purinergic signaling",
  "P2ry12", "Purinergic signaling",
  
  ## Prostaglandin signaling
  "Ptgs1", "Prostaglandin signaling",
  "Ptgs2", "Prostaglandin signaling",
  "Ptges", "Prostaglandin signaling",
  "Ptger1", "Prostaglandin signaling",
  "Ptger2", "Prostaglandin signaling",
  "Ptger3", "Prostaglandin signaling",
  "Ptger4", "Prostaglandin signaling",
  
  ## Cytokines
  "Il1b", "Cytokines",
  "Il6", "Cytokines",
  "Tnf", "Cytokines",
  "Il10", "Cytokines",
  "Tgfb1", "Cytokines",
  
  ## Chemokine signaling
  "Ccl2", "Chemokine signaling",
  "Ccl3", "Chemokine signaling",
  "Ccl4", "Chemokine signaling",
  "Ccl5", "Chemokine signaling",
  "Cxcl1", "Chemokine signaling",
  "Cxcl2", "Chemokine signaling",
  "Cxcl10", "Chemokine signaling",
  "Ccr2", "Chemokine signaling",
  "Cxcr2", "Chemokine signaling",
  
  ## G-protein signaling
  "Gnai1", "G-protein signaling",
  "Gnai2", "G-protein signaling",
  "Gnai3", "G-protein signaling",
  "Gnaq", "G-protein signaling",
  "Gnas", "G-protein signaling",
  "Gpr183", "G-protein signaling",
  "Gpr84", "G-protein signaling",
  
  ## Kinases and MAPK
  "Prkca", "Kinases and MAPK",
  "Prkcb", "Kinases and MAPK",
  "Prkcd", "Kinases and MAPK",
  "Prkce", "Kinases and MAPK",
  "Prkcg", "Kinases and MAPK",
  "Camk2a", "Kinases and MAPK",
  "Camk2d", "Kinases and MAPK",
  "Mapk1", "Kinases and MAPK",
  "Mapk3", "Kinases and MAPK",
  "Mapk8", "Kinases and MAPK",
  "Mapk14", "Kinases and MAPK",
  
  ## Calcium handling
  "Stim1", "Calcium handling",
  "Stim2", "Calcium handling",
  "Orai1", "Calcium handling",
  "Itpr1", "Calcium handling",
  "Itpr2", "Calcium handling",
  "Ryr1", "Calcium handling",
  "Ryr2", "Calcium handling",
  
  ## Voltage-gated calcium channels
  "Cacna1a", "Calcium channels",
  "Cacna1b", "Calcium channels",
  "Cacna1c", "Calcium channels",
  "Cacna2d1", "Calcium channels",
  "Cacna2d2", "Calcium channels",
  
  ## Sodium channels
  "Scn9a", "Sodium channels",
  "Scn10a", "Sodium channels",
  "Scn11a", "Sodium channels",
  
  ## Potassium channels
  "Kcnn1", "Potassium channels",
  "Kcnn2", "Potassium channels",
  "Kcnq2", "Potassium channels",
  "Kcnq3", "Potassium channels",
  
  ## Negative feedback and desensitization
  "Dusp1", "Negative feedback",
  "Dusp5", "Negative feedback",
  "Dusp6", "Negative feedback",
  "Socs1", "Negative feedback",
  "Socs3", "Negative feedback",
  "Nfkbia", "Negative feedback",
  "Tnfaip3", "Negative feedback",
  "Atf3", "Negative feedback",
  "Fos", "Immediate-early response",
  "Jun", "Immediate-early response",
  "Egr1", "Immediate-early response"
) %>%
  distinct(gene_symbol, .keep_all = TRUE)

#3. Standardize and merge the two DEG tables
esp_selected <- esp %>%
  transmute(
    gene_id,
    gene_symbol,
    gene_name,
    baseMean_ESP = baseMean,
    log2FC_ESP = log2FoldChange,
    pvalue_ESP = pvalue,
    padj_ESP = padj
  )

sge_selected <- sge %>%
  transmute(
    gene_id,
    gene_symbol,
    gene_name,
    baseMean_SGE = baseMean,
    log2FC_SGE = log2FoldChange,
    pvalue_SGE = pvalue,
    padj_SGE = padj
  )

trpv1_results <- full_join(
  esp_selected,
  sge_selected,
  by = c("gene_id", "gene_symbol"),
  suffix = c("_ESP", "_SGE")
) %>%
  left_join(
    trpv1_panel,
    by = "gene_symbol"
  ) %>%
  filter(!is.na(category))

#Check which curated genes were not detected:
missing_panel_genes <- setdiff(
  trpv1_panel$gene_symbol,
  trpv1_results$gene_symbol
)

missing_panel_genes

##4. Define evidence tiers and shared priority score
trpv1_results <- trpv1_results %>%
  mutate(
    same_direction =
      !is.na(log2FC_ESP) &
      !is.na(log2FC_SGE) &
      sign(log2FC_ESP) == sign(log2FC_SGE),
    
    shared_up =
      same_direction &
      log2FC_ESP > 0 &
      log2FC_SGE > 0,
    
    shared_down =
      same_direction &
      log2FC_ESP < 0 &
      log2FC_SGE < 0,
    
    significant_ESP =
      !is.na(padj_ESP) &
      padj_ESP < 0.05,
    
    significant_SGE =
      !is.na(padj_SGE) &
      padj_SGE < 0.05,
    
    nominal_ESP =
      !is.na(pvalue_ESP) &
      pvalue_ESP < 0.05,
    
    nominal_SGE =
      !is.na(pvalue_SGE) &
      pvalue_SGE < 0.05,
    
    significant_either =
      significant_ESP |
      significant_SGE,
    
    nominal_both =
      nominal_ESP &
      nominal_SGE,
    
    minimum_abs_log2FC = pmin(
      abs(log2FC_ESP),
      abs(log2FC_SGE),
      na.rm = TRUE
    ),
    
    mean_log2FC = rowMeans(
      cbind(log2FC_ESP, log2FC_SGE),
      na.rm = TRUE
    ),
    
    weakest_pvalue = pmax(
      replace_na(pvalue_ESP, 1),
      replace_na(pvalue_SGE, 1)
    ),
    
    weakest_padj = pmax(
      replace_na(padj_ESP, 1),
      replace_na(padj_SGE, 1)
    ),
    
    evidence_tier = case_when(
      same_direction &
        significant_ESP &
        significant_SGE ~
        "Tier 1: FDR significant in both",
      
      same_direction &
        significant_either &
        nominal_both ~
        "Tier 2: FDR in one, nominal in both",
      
      same_direction &
        nominal_both ~
        "Tier 3: nominal in both",
      
      same_direction &
        significant_either ~
        "Tier 4: FDR in one, same direction",
      
      same_direction ~
        "Tier 5: same direction only",
      
      TRUE ~
        "Discordant or unchanged"
    ),
    
    shared_direction = case_when(
      shared_up ~ "Shared up",
      shared_down ~ "Shared down",
      TRUE ~ "Discordant"
    ),
    
    ## Larger score = stronger shared evidence
    shared_priority_score =
      ifelse(same_direction, 1, 0) *
      minimum_abs_log2FC *
      (
        -log10(weakest_pvalue + 1e-300)
      )
  ) %>%
  arrange(
    desc(shared_priority_score),
    weakest_pvalue
  )
#Save the full table:
top50_trpv1 <- trpv1_results %>%
  filter(same_direction) %>%
  arrange(
    desc(shared_priority_score),
    weakest_pvalue,
    desc(minimum_abs_log2FC)
  ) %>%
  slice_head(n = 50)

#5. Select the top 50 candidates
top50_trpv1 <- trpv1_results %>%
  filter(same_direction) %>%
  arrange(
    desc(shared_priority_score),
    weakest_pvalue,
    desc(minimum_abs_log2FC)
  ) %>%
  slice_head(n = 50)

#Save a concise version:
top50_table <- top50_trpv1 %>%
  select(
    gene_symbol,
    category,
    evidence_tier,
    shared_direction,
    baseMean_ESP,
    log2FC_ESP,
    pvalue_ESP,
    padj_ESP,
    baseMean_SGE,
    log2FC_SGE,
    pvalue_SGE,
    padj_SGE,
    minimum_abs_log2FC,
    shared_priority_score
  )

write.csv(
  top50_table,
  file.path(
    outdir,
    "Tables",
    "Top50_TRPV1_regulatory_candidates.csv"
  ),
  row.names = FALSE
)

top50_table

##Summarize the evidence:
table(top50_trpv1$evidence_tier)
table(top50_trpv1$category)
table(top50_trpv1$shared_direction)

##Plot 1: ESP-versus-SGE fold-change scatter plot
scatter_labels <- trpv1_results %>%
  filter(same_direction) %>%
  arrange(desc(shared_priority_score)) %>%
  slice_head(n = 20)
p_scatter <- ggplot(
  trpv1_results,
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
    intercept = 0,
    slope = 1,
    linetype = "dotted",
    linewidth = 0.5
  ) +
  geom_point(
    aes(
      shape = shared_direction,
      size = -log10(weakest_pvalue)
    ),
    alpha = 0.75
  ) +
  ggrepel::geom_text_repel(
    data = scatter_labels,
    aes(label = gene_symbol),
    size = 3.2,
    max.overlaps = Inf,
    box.padding = 0.4,
    point.padding = 0.25,
    min.segment.length = 0
  ) +
  labs(
    title = "TRPV1-regulatory genes shared by ESP and SGE",
    subtitle = paste0(
      "Point size represents evidence from the weaker comparison;\n",
      "labels show the highest-priority shared candidates"
    ),
    x = "SGE vs Control log2 fold change",
    y = "ESP vs Control log2 fold change",
    shape = "Direction",
    size = expression(-log[10]("weakest P"))
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

ggsave(
  file.path(
    outdir,
    "Figures",
    "TRPV1_regulatory_scatter_ESP_vs_SGE.pdf"
  ),
  p_scatter,
  width = 8.5,
  height = 7.5
)

ggsave(
  file.path(
    outdir,
    "Figures",
    "TRPV1_regulatory_scatter_ESP_vs_SGE.png"
  ),
  p_scatter,
  width = 8.5,
  height = 7.5,
  dpi = 400
)

p_scatter

#Plot 2: paired log2FC plot for the top 50 genes
top50_long <- top50_trpv1 %>%
  select(
    gene_symbol,
    category,
    evidence_tier,
    shared_priority_score,
    log2FC_ESP,
    log2FC_SGE
  ) %>%
  pivot_longer(
    cols = c(log2FC_ESP, log2FC_SGE),
    names_to = "Treatment",
    values_to = "log2FC"
  ) %>%
  mutate(
    Treatment = recode(
      Treatment,
      log2FC_ESP = "ESP vs Control",
      log2FC_SGE = "SGE vs Control"
    )
  )

#Order genes using their average effect:
gene_order <- top50_trpv1 %>%
  arrange(mean_log2FC) %>%
  pull(gene_symbol)

top50_long <- top50_long %>%
  mutate(
    gene_symbol = factor(
      gene_symbol,
      levels = gene_order
    )
  )
#Create the paired plot:
p_paired <- ggplot(
  top50_long,
  aes(
    x = log2FC,
    y = gene_symbol,
    shape = Treatment
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.4
  ) +
  geom_line(
    aes(group = gene_symbol),
    linewidth = 0.45,
    alpha = 0.55
  ) +
  geom_point(
    size = 2.6,
    position = position_dodge(width = 0.25)
  ) +
  labs(
    title = "Top shared TRPV1-regulatory candidates",
    subtitle = "Paired log2 fold changes for ESP and SGE versus Control",
    x = "Log2 fold change",
    y = NULL,
    shape = "Comparison"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.y = element_text(size = 8),
    legend.position = "top"
  )

ggsave(
  file.path(
    outdir,
    "Figures",
    "Top50_TRPV1_regulatory_paired_log2FC.pdf"
  ),
  p_paired,
  width = 8.5,
  height = 12
)

ggsave(
  file.path(
    outdir,
    "Figures",
    "Top50_TRPV1_regulatory_paired_log2FC.png"
  ),
  p_paired,
  width = 8.5,
  height = 12,
  dpi = 400
)

p_paired

##Plot 3: categorized candidate plot
top_candidates_plot <- top50_trpv1 %>%
  slice_head(n = 40) %>%
  select(
    gene_symbol,
    category,
    log2FC_ESP,
    log2FC_SGE
  ) %>%
  pivot_longer(
    cols = c(log2FC_ESP, log2FC_SGE),
    names_to = "Treatment",
    values_to = "log2FC"
  ) %>%
  mutate(
    Treatment = recode(
      Treatment,
      log2FC_ESP = "ESP",
      log2FC_SGE = "SGE"
    )
  )
p_category <- ggplot(
  top_candidates_plot,
  aes(
    x = log2FC,
    y = reorder(gene_symbol, log2FC),
    shape = Treatment
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.4
  ) +
  geom_point(
    size = 2.5,
    position = position_dodge(width = 0.5)
  ) +
  facet_wrap(
    ~category,
    scales = "free_y",
    ncol = 2
  ) +
  labs(
    title = "Shared TRPV1-regulatory candidates by pathway",
    x = "Log2 fold change versus Control",
    y = NULL,
    shape = "Treatment"
  ) +
  theme_classic(base_size = 10) +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.y = element_text(size = 7),
    legend.position = "top"
  )

ggsave(
  file.path(
    outdir,
    "Figures",
    "TRPV1_regulatory_candidates_by_category.pdf"
  ),
  p_category,
  width = 11,
  height = 12
)

p_category

##Plot 4: heatmap-style fold-change plot
heatmap_df <- top50_trpv1 %>%
  select(
    gene_symbol,
    log2FC_ESP,
    log2FC_SGE
  ) %>%
  column_to_rownames("gene_symbol")

heatmap_matrix <- as.matrix(heatmap_df)

colnames(heatmap_matrix) <- c(
  "ESP vs Control",
  "SGE vs Control"
)

pheatmap::pheatmap(
  heatmap_matrix,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  border_color = NA,
  fontsize_row = 8,
  main = "Top TRPV1-regulatory candidates: log2 fold change",
  filename = file.path(
    outdir,
    "Figures",
    "Top50_TRPV1_regulatory_log2FC_heatmap.pdf"
  ),
  width = 7,
  height = 12
)

##6. Make a stricter high-priority candidate table
high_priority_trpv1 <- trpv1_results %>%
  filter(
    same_direction,
    significant_either,
    nominal_both,
    minimum_abs_log2FC >= 0.25
  ) %>%
  arrange(
    desc(shared_priority_score),
    weakest_pvalue
  )

write.csv(
  high_priority_trpv1,
  file.path(
    outdir,
    "Tables",
    "High_priority_TRPV1_regulatory_candidates.csv"
  ),
  row.names = FALSE
)

high_priority_trpv1 %>%
  select(
    gene_symbol,
    category,
    log2FC_ESP,
    padj_ESP,
    log2FC_SGE,
    padj_SGE,
    evidence_tier
  )


