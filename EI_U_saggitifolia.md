GC-MS of *U.saggitifolia*
================
Jefferson Pastuña, Pablo Corella
2023-06-14

- [Introduction](#introduction)
- [Before to start](#before-to-start)
- [Notame workflow](#notame-workflow)
- [Preprocessing](#preprocessing)
- [Second PCA and loading plot](#second-pca-and-loading-plot)

## Introduction

The present document aims to record the procedure given for the
statistical analysis of secondary metabolites present in the different
growth stages of *Urospatha saggitifolia*. For each step a brief
explanation, the code and graphics obtained are included.

The workflow used is taken from the paper [“notame”: Workflow for
Non-Targeted LC–MS Metabolic
Profiling](https://doi.org/10.3390/metabo10040135). Which offers a wide
variety of functions to perform metabolomic profile analysis.

## Before to start

The “notame” package accepts as input a feature table that can be
obtained through software such as MZMine, MSDial, among others. In this
case, the table was obtained with the help of MS-DIAL. The (\*.txt) file
was slightly modified to obtain the feature table.

Modifications made to the raw (\*.txt) file can be summarized in adding
and renaming columns. The added columns “Column” and “Ion Mode” allow to
analyze samples with different types of columns and with different
ionization modes respectively. Also, the cells corresponding to mass and
retention time must be renamed so that the package can detect and
process it.

## Notame workflow

As a first step for the analysis, all the necessary libraries were
installed and loaded in Rstudio.

``` r
library(notame)
library(Biobase)
library(BiocGenerics)
library(futile.logger)
library(ggplot2)
library(magrittr)
library(foreach)
library(iterators)
library(parallel)
library(doParallel)
library(pcaMethods)
library(patchwork)
library(cowplot)
library(Rtsne)
library(ggdendro)
library(dplyr)
library(readxl)
library(ggsci)
```

Then, a log system was added to have a record of each process executed.

``` r
init_log(log_file = "Result/GC-MS/log_GC-MS_EI.txt")
```

    ## INFO [2023-07-04 12:35:44] Starting logging

Next, the feature list was imported.

``` r
data <- read_from_excel(file = "Data/EI_feature_list.xlsx", sheet = 1, 
                        corner_row = 4, corner_column = "E", 
                        split_by = c("Column", "Ion Mode"))
```

Once the data was read, the next step was to create a MetaboSet in order
to create a specific R object.

``` r
modes <- construct_metabosets(exprs = data$exprs, 
                              pheno_data = data$pheno_data, 
                              feature_data = data$feature_data,
                              group_col = "Group")
```

Finally, each mode was extracted in a single object.

``` r
mode <- modes$Rxt5_EI
```

As a additional step, we can visualize the raw data in order to inspect
the processing routines.

``` r
Prueba_mode <- modes$Rxt5_EI
POS_raw_sambx <- plot_sample_boxplots(Prueba_mode, order_by = "Group")
POS_raw_pca <- plot_pca(Prueba_mode, center = T)
POS_raw_pca + POS_raw_sambx
```

![](EI_U_saggitifolia_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

In the boxplot we can see how the abundance of metabolites present in
the QCs vary between them, so the drift correction is necessary.

## Preprocessing

The first step of the preprocessing is to change the features with value
equal to 0 to NA.

``` r
mode <- mark_nas(mode, value = 0)
```

Then, features with low detection rate are first flagged and then will
be removed. The notame package employs two criteria to select this
features. First, is the feature presence in a percentage of QC
injections, and then the feature presence in a percentage within a
sample group or class.

``` r
mode <- flag_detection(mode, qc_limit = 0.75, group_limit = 0.9)
```

With these values, features which that were not detected in the 75% of
the QC injections and 90% of sample groups will be discarded.

The next steps for preprocessing correspond to drift correction. The
drift correction can be applied by smoothed cubic spline regression.

``` r
corrected <- correct_drift(mode)
corrected <- flag_quality(corrected)
```

Then we can visualize the correction for QCs.

``` r
EI_corr_sambx <- plot_sample_boxplots(corrected, order_by = "Group")
EI_corr_pca <- plot_pca(corrected, center = T) 
EI_corr_pca + EI_corr_sambx
```

![](EI_U_saggitifolia_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

Contaminant peaks based on the Process Blank were removed with MS-DIAL,
therefore the Process Blank group will be removed from the analysis.

``` r
corrected_no_blank <- corrected[, corrected$Group != "Blank"]
```

Finally, we can plot the PCA.

``` r
EI_PCA_2<-plot_pca(corrected_no_blank)
EI_PCA_2
```

![](EI_U_saggitifolia_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

# Second PCA and loading plot

Droping flagged features

``` r
no_flag <- drop_flagged(corrected_no_blank)

# Extracting feature table (Expression data)
peak_Height <- exprs(no_flag)

# Extracting Phenotipic data
EI_pheno_data <- no_flag@phenoData@data
```

Preparing data and transposing feature table.

``` r
EI_feat_table_pca  <- t(peak_Height)

#Changing NA to 0 
EI_feat_table_pca[is.na(EI_feat_table_pca)]=0

# Centering and Scaling features
EI_pca <- prcomp(EI_feat_table_pca, center = T, scale. = T)
```

Plotting PCA results.

``` r
scores <- EI_pca$x %>%               # Get PC coordinates
  data.frame %>%                            # Convert to data frames
  mutate(Sample_ID = rownames(.)) %>%       # Create a new column with the sample names
  left_join(EI_pheno_data )                  # Adding metadata

ggplot(scores, aes(PC1, PC2, shape = Group, color = Group)) +
  geom_point(size = 3) +
  guides(x=guide_axis(title = "PC 1 (62 %)"), y=guide_axis(title = "PC 2 (15 %)")) +
  theme_classic()
```

![](EI_U_saggitifolia_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

``` r
# Save plot
#ggsave('Result/GC-MS/GC_MS_EI_Score_Plot.png', width = 5, height = 4, device='png', dpi="print")
```

Plotting loading results.

``` r
loadings <- EI_pca$rotation %>%    # Extract loadings
  data.frame(Feature_name = rownames(.))  # New column with feat name
```

Creating an artificial table with Feature name and Compound column.

``` r
EI_feat_name <- readxl::read_excel("Data/EI_Metabolites.xlsx", 2)

# Creating a new small table of the annotated compounds
EI_compouds_all <- left_join(EI_feat_name, loadings)

# Plotting results
ggplot(loadings, aes(PC1, PC2)) + 
  geom_point(alpha = 0.2) +
  theme_classic() + 
  geom_point(data = EI_compouds_all, size = 1) +
  ggrepel::geom_label_repel(data = EI_compouds_all, aes(label = Compound), box.padding = 0.8, label.padding = 0.27, label.r = 0.3, cex = 3) +
  guides(x=guide_axis(title = "PC 1 (62 %)"), y=guide_axis(title = "PC 2 (15 %)")) +
  ggsci::scale_color_aaas()
```

![](EI_U_saggitifolia_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

``` r
# Save plot
#ggsave('Result/GC-MS/GC_MS_EI_Loadings_Plot.png', width = 5, height = 4, device='png', dpi="print")
```

Finish a record.

``` r
finish_log()
```

    ## INFO [2023-07-04 12:35:47] Finished analysis. Tue Jul  4 12:35:47 2023
    ## Session info:
    ## 
    ## INFO [2023-07-04 12:35:47] R version 4.2.0 (2022-04-22)
    ## INFO [2023-07-04 12:35:47] Platform: aarch64-apple-darwin20 (64-bit)
    ## INFO [2023-07-04 12:35:47] Running under: macOS 13.4.1
    ## INFO [2023-07-04 12:35:47] 
    ## INFO [2023-07-04 12:35:47] Matrix products: default
    ## INFO [2023-07-04 12:35:47] BLAS:   /Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/lib/libRblas.0.dylib
    ## INFO [2023-07-04 12:35:47] LAPACK: /Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/lib/libRlapack.dylib
    ## INFO [2023-07-04 12:35:47] 
    ## INFO [2023-07-04 12:35:47] locale:
    ## INFO [2023-07-04 12:35:47] [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## INFO [2023-07-04 12:35:47] 
    ## INFO [2023-07-04 12:35:47] attached base packages:
    ## INFO [2023-07-04 12:35:47] [1] parallel  stats     graphics  grDevices utils     datasets  methods  
    ## INFO [2023-07-04 12:35:47] [8] base     
    ## INFO [2023-07-04 12:35:47] 
    ## INFO [2023-07-04 12:35:47] other attached packages:
    ## INFO [2023-07-04 12:35:47]  [1] ggsci_3.0.0         readxl_1.4.2        dplyr_1.1.2        
    ## INFO [2023-07-04 12:35:47]  [4] ggdendro_0.1.23     Rtsne_0.16          cowplot_1.1.1      
    ## INFO [2023-07-04 12:35:47]  [7] patchwork_1.1.2     pcaMethods_1.90.0   doParallel_1.0.17  
    ## INFO [2023-07-04 12:35:47] [10] iterators_1.0.14    foreach_1.5.2       notame_0.2.1       
    ## INFO [2023-07-04 12:35:47] [13] magrittr_2.0.3      ggplot2_3.4.2       futile.logger_1.4.3
    ## INFO [2023-07-04 12:35:47] [16] Biobase_2.58.0      BiocGenerics_0.44.0
    ## INFO [2023-07-04 12:35:47] 
    ## INFO [2023-07-04 12:35:47] loaded via a namespace (and not attached):
    ## INFO [2023-07-04 12:35:47]  [1] tidyselect_1.2.0     xfun_0.39            purrr_1.0.1         
    ## INFO [2023-07-04 12:35:47]  [4] colorspace_2.1-0     vctrs_0.6.3          generics_0.1.3      
    ## INFO [2023-07-04 12:35:47]  [7] usethis_2.2.1        htmltools_0.5.5      viridisLite_0.4.2   
    ## INFO [2023-07-04 12:35:47] [10] yaml_2.3.7           utf8_1.2.3           rlang_1.1.1         
    ## INFO [2023-07-04 12:35:47] [13] gert_1.9.2           pillar_1.9.0         glue_1.6.2          
    ## INFO [2023-07-04 12:35:47] [16] withr_2.5.0          RColorBrewer_1.1-3   lambda.r_1.2.4      
    ## INFO [2023-07-04 12:35:47] [19] lifecycle_1.0.3      munsell_0.5.0        gtable_0.3.3        
    ## INFO [2023-07-04 12:35:47] [22] cellranger_1.1.0     zip_2.3.0            codetools_0.2-19    
    ## INFO [2023-07-04 12:35:47] [25] evaluate_0.21        labeling_0.4.2       knitr_1.43          
    ## INFO [2023-07-04 12:35:47] [28] fastmap_1.1.1        sys_3.4.2            fansi_1.0.4         
    ## INFO [2023-07-04 12:35:47] [31] highr_0.10           Rcpp_1.0.10          openssl_2.0.6       
    ## INFO [2023-07-04 12:35:47] [34] scales_1.2.1         formatR_1.14         farver_2.1.1        
    ## INFO [2023-07-04 12:35:47] [37] fs_1.6.2             credentials_1.3.2    askpass_1.1         
    ## INFO [2023-07-04 12:35:47] [40] digest_0.6.31        stringi_1.7.12       openxlsx_4.2.5.2    
    ## INFO [2023-07-04 12:35:47] [43] ggrepel_0.9.3        grid_4.2.0           cli_3.6.1           
    ## INFO [2023-07-04 12:35:47] [46] tools_4.2.0          tibble_3.2.1         futile.options_1.0.1
    ## INFO [2023-07-04 12:35:47] [49] tidyr_1.3.0          pkgconfig_2.0.3      MASS_7.3-60         
    ## INFO [2023-07-04 12:35:47] [52] rmarkdown_2.23       rstudioapi_0.14      R6_2.5.1            
    ## INFO [2023-07-04 12:35:47] [55] compiler_4.2.0
