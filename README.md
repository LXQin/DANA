
# DANA

<!-- badges: start -->
<!-- badges: end -->

DANA is an approach for assessing the performance of normalization for microRNA-Seq data based on biology-motivated and data-driven metrics.
Our approach takes advantage of well-known biological features of microRNAs for their expression pattern and polycistronic clustering to assess (1) how effectively normalization removes handling effects and (2) how normalization biases true biological signals.
DANA is implemented in R and can be used for assessing any normalization method (under minimal assumptions) for any microRNA-Seq data set and only requires additional information on polycistronic clustering or marker location on the genome, which is typically readily available.

In general, you can apply the DANA assessment for any normalization method as long as the normalized counts are non-negative.
For simplicity, we provide eight commonly used normalization methods:

 - Total Count (TC)
 - Upper Quartile (UQ)
 - Median (median)
 - Trimmed Median of Means (TMM)
 - DESeq
 - Quantile Normalization (QN)
 - Remove Unwanted Variation (with subtypes: RUVg, RUVr, and RUVs)
 - PoissonSeq
 
For each normalization method under study, DANA computes two assessment metrics: `cc` and `mscr`:
 
 - `cc` measures the preservation of biological signals before versus after normalization. A high value indicates a high preservation of biological signals (`cc` <= 1).
 - `mscr` measures the relative reduction of handling before versus after normalization. A high `mscr` indicates higher removal of handling effects.

An optimal normalization maximally removes handling effects (high `mscr`) while keeping biological signals intact (`cc` close to 1).
However, in most cases, there is no clear "best" method with maximal `mscr` _and_ maximal `cc`.
Therefore, one should aim for the best possible trade-off between the proposed statistics for negative and positive controls with an emphasis on keeping biological signals intact (`cc` close to 1).
The two metrics can be easily assessed by plotting the metrics in a scatter plot (see function `plotDANA`) for each normalization method under study, where a preferable method should be located towards the top-right quadrant of the plot.



## Installation

You can install the released version of DANA directly from GitHub using devtools by:

```R
# install.packages("devtools")
devtools::install_github("LXQin/DANA")
```

The R package `PoissonSeq` for PoissonSeq normalization was removed from CRAN, but you can install the archived version from github using:

```R
devtools::install_github("cran/PoissonSeq")
```

## Usage

Please refer to the [documentation page](https://lxqin.github.io/DANA/) for detailed function references.

Generally, users need to provide a raw count matrix from RNA-Seq studies (without any normalization or transformation) and information about polycistronic clustering of miRNAs.
Clustering information can be readily computed using the `defineClusters` function if miRNA chromosome and location information (base-pair/nucleotide location) is available.
For data using the notation from miRBase (v22), such as data from TCGA, DANA includes all necessary information.
For data using other notation, the user must provide the chromosome and nucleotide location on each chromosome for each miRNA in the data as inputs `chr` and `pos`, respectively, to the `defineClusters` function.

Assume that `raw.counts` hold the raw un-normalized count data where row names correspond to genes and column names correspond to sample names.
`groups` is a vector of sample groups (e.g. sample sub-types).

First, normalize the raw count data using the provided normalization routine:

```R
normalized <- applyNormalization(
  raw.counts,
  groups,
  method = c("TC", "UQ", "median", "TMM", "DESeq", "PoissonSeq", "QN", "RUV"))
```

Next, define polycistronic clusters for all miRNAs.
For miRNA data not using miRBase (v22) notation, `chr` and `pos` provide information about the chromosome and location of each marker on the chromosome in `raw.counts`. 

```R
# miRBase (v22) notation
clusters <- defineClusters(rownames(raw.counts))

# other notation -> provide "chr" and "pos"
clusters <- defineClusters(rownames(raw.counts), chr, pos)
```
Define and positive and negative control markers using the `defineControls` function. 
Positive control markers are well-expressed (in [tWell, inf)), clustered markers representing biological effects.
Negative control markers are poorly-expressed (in [tZero, tPoor]) and primarily represent handling effects.
Adjust the bounds for positive and negative controls to your data set, e.g. using the helper plot functions `plotMeanSD` and `plotCountHist`.

```R
controls <- defineControls(
  raw.counts, 
  tZero = 2,     # lower count bound for negative controls
  tPoor = 5,     # upper count bound for negative controls
  tWell = 100,   # lower count bound for positive controls
  clusters)
```

Use DANA to assess the normalization for your data set:

```R
res <- assessNormalization(
  raw=raw.counts,
  normalized  = normalized,
  negControls = controls$negControls,
  posControls = controls$posControls,
  clusters    = clusters)
```

Finally, assess normalization as stated above and use the `plotDANA` function to compare the result metrics `cc` and `mscr` for each normalization method.

```R
plotDANA(res)
```







