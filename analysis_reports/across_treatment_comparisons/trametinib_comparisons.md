Trametinib Treatment Single Cell Analysis
================
Matthew Berginski
Thu Jul 23 17:15:47 2020

## Data Cleaning/Normalization

The single cell data sets were pre-processed using cellranger and count
matrices produced. Starting with the count matrices, I removed cells
with especially high mitochondrial gene counts:

  - DMSO: 135 removed / 2649 total
  - Trametinib: 42 removed / 678 total

One clear problem is the imbalance in the number of cells captured in
each run. This shouldn’t affect the downstream analysis, but might limit
the differences we can detect.

Let’s take a look at the number of reads assigned to each cell,
splitting the data by the treatment type:

![](trametinib_comparisons_files/figure-gfm/read%20counts%20by%20cell-1.png)<!-- -->

Nothing too surprising here, the fewer Trametinib cells has resulted in
a higher average sequencing depth per cell. This might prove to be a
problem in dealing with very low count genes (like the kinases). There
aren’t any clear outliers in terms of the cells present, so we won’t do
any other cell filtering. Now I’ll normalize each cell to the 75
percentile of the reads per cell to allow comparisons between the counts
in each cell.

## Kinases Present in Normalized Data

Just a quick look at the kinases that have some evidence of expression.
I’m setting this level at seeing at least a normalized count of 1 or
greater in a specific cell. Every cell in the DMSO and Trametinib
treatments have a kinase with at least 1 read. Nearly all of the kinases
in my master list are present in at least one cell with at least one
count (564 total, 159 dark kinases). Here’s the distribution of kinase
counts on a per cell basis:

![](trametinib_comparisons_files/figure-gfm/number%20of%20kinases-1.png)<!-- -->

Overall the distributions look pretty similar for the kinases between
DMSO and Trametinib. Let’s move onto trying to detect different
expression patterns.

## Genome Wide Expression Changes

It’s sort of the wild-west when it comes to calling gene expression
changes from single cell sequencing data, so I pieced together some
techniques that made the most sense. First off, all the comparisons are
made with the per-cell normalized count values. Otherwise, I used the
wilcox test to look for genes with differential expression levels. I
also calculated the log 2 fold change in average expression, which
raises a bit of a complication for single cell sequencing. Since the
sequencing depth per cell isn’t very high, low prevalence transcripts
are often missed. In practice, this means I’ve got lots of mean zero
counts in both the DMSO and Trametinib, but the resulting infinite fold
change values aren’t all that interesting. To combat these infinite
values, I decided to toss out the all the fold changes when either
treatment was zero, except when the other treatment was showing up in at
least 5 percent of the cell. I choose 5 percent somewhat arbitrarily,
but I figure if there is no evidence of expression while the gene shows
up in 5 percent of cells in the other treatment, that gene shouldn’t be
removed.

### Number of Diff Expressed Genes

Overall, there are 5672 genes deferentially expressed. Of these, 154 are
kinases and 45 are understudied. To get an idea of what these
differences actually look like I’ll produce a few plots from some of the
dark kinases with larger changes.

### Sample Dark Kinases that Change

PKMYT1 happens to one of the biggest drops after trametinib treatment:

![](trametinib_comparisons_files/figure-gfm/PKMYT1%20expression-1.png)<!-- -->

RIOK1 also drops:

![](trametinib_comparisons_files/figure-gfm/RIOK1%20expression-1.png)<!-- -->

ALPK3 goes up:

![](trametinib_comparisons_files/figure-gfm/ALPK3%20expression-1.png)<!-- -->

ULK4 also goes up:

![](trametinib_comparisons_files/figure-gfm/ULK4%20expression-1.png)<!-- -->

### All Dark Kinase Plots

Here’s a small multiples plot with all the differentially expression
dark kinases for reference:

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

![](trametinib_comparisons_files/figure-gfm/all%20dark%20kinases-1.png)<!-- -->

## Bulk Comparisons

Now I’ll try to compare these results to the DESeq2 based results. In
general the single cell analysis is finding more differences at the gene
level, with only a few genes being only picked up by DESeq2.

### Overall Fold Change Comparison

The overall correlation is 0.78. The red line indicates the one to one
line:

![](trametinib_comparisons_files/figure-gfm/overall%20fold%20change-1.png)<!-- -->

### Overlap with All Differentially Expressed Genes

![](trametinib_comparisons_files/figure-gfm/all%20genes%20upset-1.png)<!-- -->

### Overlap with Kinase Differentially Expressed Genes

![](trametinib_comparisons_files/figure-gfm/just%20kinases%20upset-1.png)<!-- -->

### Overlap with Dark Kinase Differentially Expressed Genes

![](trametinib_comparisons_files/figure-gfm/just%20dark%20kinases%20upset-1.png)<!-- -->
