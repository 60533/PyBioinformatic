**<a id="link1">User guide of RNA-seq</a>**<br />

- [1. RNA-seq analysis](#1-rna-seq-analysis)
  - [1.1 Standardize the expression levels of genes](#11-standardize-the-expression-levels-of-genes)
  - [1.2 Differential expression analysis](#12-differential-expression-analysis)
  - [1.3 Allele expression pattern clustering](#13-allele-expression-pattern-clustering)

# 1. RNA-seq analysis
## 1.1 Standardize the expression levels of genes
```python
from pandas import read_csv
from pybioinformatic import get_FPKM, get_TPM

# Data reading and preprocessing
reads_count = read_csv(
    'test/RNA-seq/reads_count.txt',
    sep='\t',
    index_col=0
).iloc[:, 4:]

# calculate FPKM and TPM
FPKM = get_FPKM(reads_count, min_value=0)
print(FPKM.head())

TPM = get_TPM(reads_count, min_value=0)
print(TPM.head())
```
```text
# PFKM
                                 leaf  ...         xylem
Geneid                                 ...              
Posim.A01G000100.1.v1.0   9947.798049  ...  20539.057139
Posim.A01G000100.2.v1.0   8448.458192  ...  16840.578670
Posim.A01G000100.3.v1.0   8758.543619  ...  17445.378219
Posim.A01G000100.4.v1.0   7497.544554  ...  14347.661828
Posim.A01G000200.1.v1.0  36334.326673  ...  48958.187516

[5 rows x 4 columns]
```
```text
# TPM
                                 leaf  ...          xylem
Geneid                                 ...               
Posim.A01G000100.1.v1.0  23745.228638  ...   53609.408586
Posim.A01G000100.2.v1.0  20166.329315  ...   43955.935106
Posim.A01G000100.3.v1.0  20906.498076  ...   45534.534646
Posim.A01G000100.4.v1.0  17896.514263  ...   37449.122421
Posim.A01G000200.1.v1.0  86729.433995  ...  127786.755760

[5 rows x 4 columns]
```

## 1.2 Differential expression analysis
```python
from pybioinformatic import GeneExpressionAnalysis


gea = GeneExpressionAnalysis('test/RNA-seq/reads.count.xls')
gea.DESeq2(
    metadata='test/RNA-seq/DEG/metadata.xls',
    min_reads_count=3,
    padj=0.05,
    log2fc=1,
    figure_size=(10, 10),
    out_path='test/RNA-seq/DEG'
)
```
![image](https://github.com/60533/PyBioinformatic/blob/main/test/RNA-seq/DEG/OE_vs_CK.volcano.png)

## 1.3 Allele expression pattern clustering
```python
from pybioinformatic import AllelicExpressionAnalyzer

aea = AllelicExpressionAnalyzer(
    allele_exp_file='FPKM.xls',
    allele_pairs_file='allele_pairs.xls',  # three columns: "locus_id", "allele1_id", and "allele2_id"
    num_processing=10
)

# Filter out low-expression alleles and construct an allele pair feature matrix.
aea.preprocess_data(
    min_expression=0.5,  # Minimum expression value threshold
    min_samples=1  # How many samples at least need to be expressed
)

# Alleles expression patterns clustering.
aea.clustering(
    method='kmeans',   # Clustering method: {'kmeans', 'hierarchical', 'dbscan', 'agglomerative'}
    n_clusters=4,  # Number of clusters (it is only applicable to methods that require specifying the number of clusters)
    features_to_use=None  # A list of features for clustering: {'mean_ratio', 'fold_change', 'pearson_corr', 't_statistic', 'dominance'}
)

# Visualize allele clustering results.
aea.visualize_clusters(
    method='pca',
    figsize=(12, 10),
    out_file='./clustering_results.png'
)

# Identify the difference patterns between two clusters.
aea.identify_differential_patterns(
    cluster_id1=0,
    cluster_id2=1,
    alpha=0.05
)
aea.identify_differential_patterns(
    cluster_id1=0,
    cluster_id2=2,
    alpha=0.05
)
aea.identify_differential_patterns(
    cluster_id1=0,
    cluster_id2=3,
    alpha=0.05
)

aea.save_clustering_results(output_prefix='./allele')
```
**Output**
- `<output_prefix>_pairwise_features.xls`: It contains the quantitative characteristics and clustering labels of each allele pair.
  - `loci_id`: loci id
  - `mean_ratio`: mean(allele1) / mean(allele2)
  - `fold_change`: log2(mean(allele1)/mean(allele2))
  - `pearson_corr`: pearson correlation coefficient
  - `t_statistic`: Welch's t-test statistic
  - `ttest_pvalue`: t-test p value
  - `dominance`: allele1 dominance "mean(allele1 > allele2)"
  - `expression_balance`: expression the degree of imbalance "abs(fold_change)"
  - `[method]_cluster`: cluster label
- `<output_prefix>_clustering_report.txt`: The process and effect evaluation of cluster analysis were recorded.
- `<output_prefix>_differential_patterns.xlsx`: This is generated after executing the function `identify_differential_patterns`, and it compares the feature differences between the two clusters.
```python
# Obtain all alleles in the specified cluster.
alleles = aea.get_cluster_alleles(cluster_id=1)
print(alleles)
```
```text
['Posim.A01G000100.v1.0', 'Posim.A01G002800.v1.0', ...]
```
