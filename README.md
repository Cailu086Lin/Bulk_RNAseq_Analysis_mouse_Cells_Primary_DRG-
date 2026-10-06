Analysis provided in the manuscript titled "Hookworms attenuate Substance P release from sensory neurons to promote host susceptibility" This repository is meant to provide the source code for the analysis of bulk RNASeq data provided in the above manuscript by
Adriana Stephenson, Ulrich Femoe, Parvathi Annamalai, Fungai Musaigwa, Juan Inclan-Rico, Cailu Lin, Danielle R. Reed, Camila M. Napuri, and De’Broski R. Herbert

Keywords: TRPV1, helminth, neutrophils, Substance P

Abstract15
Whether parasitic nematodes subvert function(s) of skin sensory neurons is unknown. Data show16
that acute infection of mice with the hookworm Nippostrongylus brasiliensis blunts the perception17
of pain. Excretory/secretory products isolated from infectious third stage larvae (iL3) antagonized18
TRPV1-dependent calcium mobilization and neuropeptide release from primary neurons without19
causing cytotoxicity. Indeed, bulk RNA sequencing reveals that heat labile components of ESP20
induce an inflammatory transcriptomic profile in primary dorsal root ganglia cultures, with gene21
upregulation suggestive of primary cells responding to ESP and transitioning into a highly22
regulatory state – perhaps to prevent neuronal cytotoxic excitation. Moreover, preemptive23
activation of TRPV1+ skin afferents reduced the dissemination of N. brasiliensis larvae from skin24
to the pulmonary tract. Conversely, genetic deficiency in the high-affinity Substance P receptor25
(Tacr1) increased parasite dissemination and resulted in greater gamma delt T cell/neutrophil responses as26
compared to littermate controls. These data indicate that parasitic nematodes larvae evolved27
strategies to antagonize neurogenic inflammation to favor host parasitis.

sessionInfo()

R version 4.6.1 (2026-06-24)
Platform: x86_64-pc-linux-gnu
Running under: Ubuntu 24.04.4 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8       
 [4] LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
[10] LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

time zone: America/New_York
tzcode source: system (glibc)

attached base packages:
[1] grid      stats4    stats     graphics  grDevices utils     datasets  methods  
[9] base     

other attached packages:
 [1] EnhancedVolcano_1.30.0      circlize_0.4.18            
 [3] ComplexHeatmap_2.28.0       pheatmap_1.0.13            
 [5] ggrepel_0.9.8               ggplot2_4.0.3              
 [7] readr_2.2.0                 tidyr_1.3.2                
 [9] tibble_3.3.1                dplyr_1.2.1                
[11] clusterProfiler_4.20.0      org.Mm.eg.db_3.23.0        
[13] AnnotationDbi_1.74.0        rtracklayer_1.72.0         
[15] DESeq2_1.52.0               SummarizedExperiment_1.42.0
[17] Biobase_2.72.0              MatrixGenerics_1.24.0      
[19] matrixStats_1.5.0           GenomicRanges_1.64.0       
[21] Seqinfo_1.2.0               IRanges_2.46.0             
[23] S4Vectors_0.50.1            BiocGenerics_0.58.1        
[25] generics_0.1.4              tximport_1.40.0            

loaded via a namespace (and not attached):
  [1] splines_4.6.1            later_1.4.8              BiocIO_1.22.0           
  [4] bitops_1.0-9             ggplotify_0.1.3          polyclip_1.10-7         
  [7] enrichit_0.2.0           XML_3.99-0.23            lifecycle_1.0.5         
 [10] httr2_1.3.0              doParallel_1.0.17        vroom_1.7.1             
 [13] processx_3.9.0           lattice_0.22-9           MASS_7.3-65             
 [16] magrittr_2.0.5           rmarkdown_2.31           yaml_2.3.12             
 [19] httpuv_1.6.17            otel_0.2.0               ggtangle_0.1.2          
 [22] DBI_1.3.0                RColorBrewer_1.1-3       abind_1.4-8             
 [25] purrr_1.2.2              RCurl_1.98-1.19          yulab.utils_0.2.4       
 [28] tweenr_2.0.3             rappdirs_0.3.4           aisdk_1.4.12            
 [31] gdtools_0.5.1            enrichplot_1.32.0        tidytree_0.4.8          
 [34] codetools_0.2-20         DelayedArray_0.38.2      DOSE_4.6.0              
 [37] ggforce_0.5.0            shape_1.4.6.1            tidyselect_1.2.1        
 [40] aplot_0.3.1              farver_2.1.2             GenomicAlignments_1.48.0
 [43] jsonlite_2.0.0           GetoptLong_1.1.1         iterators_1.0.14        
 [46] systemfonts_1.3.2        foreach_1.5.2            tools_4.6.1             
 [49] ggnewscale_0.5.2         ragg_1.5.2               treeio_1.36.1           
 [52] Rcpp_1.1.2               glue_1.8.1               SparseArray_1.12.2      
 [55] xfun_0.60                qvalue_2.44.0            withr_3.0.3             
 [58] fastmap_1.2.0            callr_3.8.0              digest_0.6.39           
 [61] R6_2.6.1                 mime_0.13                gridGraphics_0.5-1      
 [64] textshaping_1.0.5        colorspace_2.1-3         GO.db_3.23.1            
 [67] RSQLite_3.53.3           cigarillo_1.2.1          utf8_1.2.6              
 [70] fontLiberation_0.1.0     data.table_1.18.4        httr_1.4.8              
 [73] htmlwidgets_1.6.4        S4Arrays_1.12.0          scatterpie_0.2.6        
 [76] pkgconfig_2.0.3          gtable_0.3.6             rsconnect_1.10.1        
 [79] blob_1.3.0               S7_0.2.2                 XVector_0.52.0          
 [82] htmltools_0.5.9          fontBitstreamVera_0.1.1  clue_0.3-68             
 [85] scales_1.4.0             png_0.1-9                ggfun_0.2.1             
 [88] knitr_1.51               rstudioapi_0.19.0        tzdb_0.5.0              
 [91] reshape2_1.4.5           rjson_0.2.23             nlme_3.1-169            
 [94] curl_7.1.0               cachem_1.1.0             GlobalOptions_0.1.4     
 [97] stringr_1.6.0            parallel_4.6.1           restfulr_0.0.17         
[100] pillar_1.11.1            vctrs_0.7.3              promises_1.5.0          
[103] tidydr_0.0.6             xtable_1.8-8             cluster_2.1.8.2         
[106] evaluate_1.0.5           cli_3.6.6                locfit_1.5-9.12         
[109] compiler_4.6.1           Rsamtools_2.28.0         rlang_1.3.0             
[112] crayon_1.5.3             labeling_0.4.3           ps_1.9.3                
[115] plyr_1.8.9               fs_2.1.0                 ggiraph_0.9.6           
[118] stringi_1.8.7            BiocParallel_1.46.0      Biostrings_2.80.1       
[121] lazyeval_0.2.3           pacman_0.5.1             GOSemSim_2.38.3         
[124] fontquiver_0.2.1         Matrix_1.7-5             hms_1.1.4               
[127] patchwork_1.3.2          bit64_4.8.2              KEGGREST_1.52.2         
[130] shiny_1.14.0             igraph_2.3.3             memoise_2.0.1           
[133] ggtree_4.3.0             bit_4.6.0                ape_5.8-1               
[136] gson_0.2.0  
