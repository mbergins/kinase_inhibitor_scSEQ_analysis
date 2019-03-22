Kinase Inhibitor Single Cell Analysis
================
Matthew Berginski
February 7, 2019

Read In Data
============

``` r
all_data = rbind(
  #This is a Note
  read.delim(here('raw_data/matrix_DMSO_2.txt.gz')) %>%
    rename(gene_name = X) %>% 
    gather("cell_id","count",-gene_name) %>%
    mutate(treatment = "DMSO"),
  
  read.delim(here('raw_data/matrix_JQ1_2.txt.gz')) %>%
    rename(gene_name = X) %>%
    gather("cell_id","count",-gene_name) %>%
    mutate(treatment = "JQ1"),
  
  read.delim(here('raw_data/matrix_TRAMET_2.txt.gz')) %>%
    rename(gene_name = X) %>%
    gather("cell_id","count",-gene_name) %>%
    mutate(treatment = "GSK"),
  
  read.delim(here('raw_data/matrix_COMBO_2.txt.gz')) %>%
    rename(gene_name = X) %>%
    gather("cell_id","count",-gene_name) %>%
    mutate(treatment = "combo")
)
```

``` r
reads_per_cell = all_data %>%
  group_by(cell_id,treatment) %>%
  summarise(read_count = sum(count)) %>%
  group_by(treatment) %>%
  summarise(treatment_cell_q3 = quantile(read_count,c(0.75)))

quantile_average = mean(reads_per_cell$treatment_cell_q3)

all_data = all_data %>% left_join(reads_per_cell)
```

    ## Joining, by = "treatment"

``` r
all_data = all_data %>%
  mutate(count_normed = count*quantile_average/treatment_cell_q3)
```

Data Cleaning?
==============

![](scAnalysis_files/figure-markdown_github/data_cleaning_filtering-1.png)![](scAnalysis_files/figure-markdown_github/data_cleaning_filtering-2.png)

``` r
reads_per_cell = all_data %>%
  group_by(cell_id,treatment) %>%
  summarise(read_count = sum(count),
            read_count_normed = sum(count_normed))

ggplot(reads_per_cell,aes(x=read_count_normed, color=treatment)) + geom_freqpoly()
```

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

![](scAnalysis_files/figure-markdown_github/unnamed-chunk-2-1.png)

Summary Statistics
==================

``` r
gene_summary = all_data %>%
  group_by(gene_name,treatment) %>%
  summarise(mean_expression = mean(count_normed),
            sd_expression = sd(count_normed),
            percent_cells = 100*(sum(count_normed > 0)/length(unique(cell_id))))
```

``` r
#single gene filter
# temp = gene_summary %>% filter(gene_name == "ACTB")

temp_DMSO = gene_summary %>% 
  filter(treatment == "DMSO") %>% 
  select(-treatment) %>% 
  rename(mean_dmso = mean_expression,
         sd_dmso = sd_expression,
         percent_dmso = percent_cells)

gene_summary = left_join(gene_summary,temp_DMSO)
```

    ## Joining, by = "gene_name"

``` r
gene_summary = gene_summary %>%
  mutate(mean_dmso_diff = mean_expression - mean_dmso,
         sd_dmso_diff = sd_expression - sd_dmso,
         percent_diff = percent_cells - percent_dmso)
```

``` r
gene_summary = gene_summary %>% arrange(desc(sd_dmso_diff))

ggplot(gene_summary, aes(x=sd_dmso_diff)) + geom_histogram(bins=50)
```

![](scAnalysis_files/figure-markdown_github/unnamed-chunk-5-1.png)
