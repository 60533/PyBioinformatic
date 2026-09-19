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

reads_count = read_csv(
    'test/RNA-seq/reads.count.fc.xls',
    sep='\t',
    index_col=0
)

# calculate FPKM and TPM
FPKM = get_FPKM(reads_count, min_value=0)
print(FPKM)

TPM = get_TPM(reads_count, min_value=0)
print(TPM)
```
```text
# PFKM
                    CK_1        CK_2  ...       KO_2       KO_3
Geneid                                ...                      
Pe.001G000100   0.000000    0.000000  ...   0.000000   0.000000
Pe.001G000200   0.000000    0.000000  ...   0.000000   0.000000
Pe.001G000300  35.229230   36.196582  ...  46.798874  34.981561
Pe.001G000400  31.993198   44.485188  ...  31.840140  40.098900
Pe.001G000500  30.067137   25.076738  ...  33.889950  17.274440
...                  ...         ...  ...        ...        ...
Pe.019G115500   0.244846    0.433143  ...   0.049441   0.096276
Pe.019G115600  94.832151  149.217301  ...  65.517412  64.933021
Pe.019G115700   0.000000    0.000000  ...   0.000000   0.000000
Pe.019G115800   0.092146    0.073355  ...   0.156298   0.043479
Pe.019G115900   0.000000    0.000000  ...   0.000000   0.000000

[36436 rows x 9 columns]
```
```text
# TPM
                     CK_1        CK_2  ...        KO_2        KO_3
Geneid                                 ...                        
Pe.001G000100    0.000000    0.000000  ...    0.000000    0.000000
Pe.001G000200    0.000000    0.000000  ...    0.000000    0.000000
Pe.001G000300   68.737634   69.933088  ...   96.429703   69.664905
Pe.001G000400   62.423638   85.946972  ...   65.607031   79.855958
Pe.001G000500   58.665598   48.449153  ...   69.830694   34.401615
...                   ...         ...  ...         ...         ...
Pe.019G115500    0.477732    0.836849  ...    0.101874    0.191731
Pe.019G115600  185.032078  288.293152  ...  134.999498  129.312489
Pe.019G115700    0.000000    0.000000  ...    0.000000    0.000000
Pe.019G115800    0.179792    0.141724  ...    0.322054    0.086588
Pe.019G115900    0.000000    0.000000  ...    0.000000    0.000000

[36436 rows x 9 columns]
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

## 1.3 Differential expression enrichment analysis
```python
from pybioinformatic import GeneExpressionAnalysis

gea = GeneExpressionAnalysis(exp_matrix='test/RNA-seq/reads.count.xls')
ret = gea.enrichment_analysis(
    foreground_gene_set='test/RNA-seq/enrichment/DEGs.lst',
    background_gene_set='test/RNA-seq/enrichment/KEGG_anno.xls',
    min_exp=3,
    n_top=15,
    title='OE_vs_CK',
    color_map='RdYlBu',
    wrap_width=35,
    figure_size=(6, 8),
    out_path='test/RNA-seq/enrichment'
)
```
![image](https://github.com/60533/PyBioinformatic/blob/main/test/RNA-seq/enrichment/dotplot.png)

## 1.4 Allele expression pattern clustering
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
