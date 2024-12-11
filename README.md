# Experimental evolution of host range for two isolates of potato virus Y

### Description
R scripts used in Morais et al.

### Overview
To run the analysis, processed data (PVYNb and PVYSt consensus sequences, lofreq output and coverage depth per sample) must be obtained from the [Zenodo repository](https://zenodo.org/) and the contents of the tar.gz file must be extracted here.

#### Contents
##### source_me.R
Loading of the data and main functions.

##### fluctuations.R
Analyses of allele fluctuations and population diversity across the experiments.

The analysis was performed on the following R environment:

```
> sessionInfo()
R version 4.4.1 (2024-06-14)
Platform: aarch64-apple-darwin22.6.0
Running under: macOS Ventura 13.1

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
LAPACK: /opt/homebrew/Cellar/r/4.4.1/lib/R/lib/libRlapack.dylib;  LAPACK version 3.12.0

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: Europe/Madrid
tzcode source: internal

attached base packages:
[1] grid      stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] seqinr_4.2-36     factoextra_1.0.7  cluster_2.1.6     Hmisc_5.1-3       gridExtra_2.3     statcomp_0.1.0    pracma_2.4.4      minpack.lm_1.2-4  igraph_2.0.3     
[10] randtests_1.0.2   ggrepel_0.9.5     viridis_0.6.5     viridisLite_0.4.2 ggpubr_0.6.0      segmented_2.1-1   nlme_3.1-164      MASS_7.3-60.2     readxl_1.4.3     
[19] writexl_1.5.0     ggplot2_3.5.1     tidyr_1.3.1       dplyr_1.1.4      

loaded via a namespace (and not attached):
 [1] gtable_0.3.5      xfun_0.47         htmlwidgets_1.6.4 rstatix_0.7.2     lattice_0.22-6    vctrs_0.6.5       tools_4.4.1       generics_0.1.3    tibble_3.2.1     
[10] fansi_1.0.6       pkgconfig_2.0.3   data.table_1.16.0 checkmate_2.3.2   lifecycle_1.0.4   farver_2.1.2      stringr_1.5.1     compiler_4.4.1    munsell_0.5.1    
[19] carData_3.0-5     htmltools_0.5.8.1 htmlTable_2.4.3   Formula_1.2-5     crayon_1.5.3      pillar_1.9.0      car_3.1-2         rpart_4.1.23      abind_1.4-5      
[28] tidyselect_1.2.1  digest_0.6.37     stringi_1.8.4     purrr_1.0.2       labeling_0.4.3    splines_4.4.1     ade4_1.7-22       cowplot_1.1.3     fastmap_1.2.0    
[37] colorspace_2.1-1  cli_3.6.3         magrittr_2.0.3    base64enc_0.1-3   utf8_1.2.4        broom_1.0.6       foreign_0.8-86    withr_3.0.1       scales_1.3.0     
[46] backports_1.5.0   rmarkdown_2.28    nnet_7.3-19       ggsignif_0.6.4    cellranger_1.1.0  evaluate_0.24.0   knitr_1.48        rlang_1.1.4       Rcpp_1.0.13      
[55] glue_1.7.0        rstudioapi_0.16.0 R6_2.5.1    
```

