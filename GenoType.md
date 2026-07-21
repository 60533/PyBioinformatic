**<a id="link1">User guide of GenoType class</a>**<br />
- [1. Genotype file format](#1-genotype-file-format)
- [2. Instantiate a GenoType object and convert it to the pandas.DataFrame class](#2-instantiate-a-genotype-object-and-convert-it-to-the-pandasdataframe-class)
- [3. Data preprocessing](#3-data-preprocessing)
- [4. Merge data by position](#4-merge-data-by-position)
- [5. Filter mendelian errors](#5-filter-mendelian-errors)
- [6. Calculate the missing rate, heterozygosity rate and minor allele frequency (MAF) of SNP](#6-calculate-the-missing-rate-heterozygosity-rate-and-minor-allele-frequency-maf-of-snp)
- [7. Identify the SNPs that can distinguish different sample groups](#7-identify-the-snps-that-can-distinguish-different-sample-groups)
- [8. Find the smallest set of SNP loci that can distinguish all samples by greedy forward choice](#8-find-the-smallest-set-of-snp-loci-that-can-distinguish-all-samples-by-greedy-forward-choice)
- [9. Calculate genotype similarity](#9-calculate-genotype-similarity)
  - [9.1. Compare with the genotype matrix of the database. It can be used for cultivar identification](#91-compare-with-the-genotype-matrix-of-the-database-it-can-be-used-for-cultivar-identification)
  - [9.2. Compare the genotype consistency among different replicate samples](#92-compare-the-genotype-consistency-among-different-replicate-samples)
- [10. Draw genotype heatmap](#10-draw-genotype-heatmap)
- [11. Build finger-print](#11-build-finger-print)
  - [11.1. Select core SNPs](#111-select-core-snps)
  - [11.2. Visualized fingerprint](#112-visualized-fingerprint)
- [12. Analysis of phenotypic differences by genotype](#12-analysis-of-phenotypic-differences-by-genotype)

# 1. Genotype file format
**A standard genotype file should at least contain 5 columns, with each column separated by a tab. The column names are `ID`, `Chrom`, `Position`, `Ref` and `Sample_name`.**
```text
ID     Chrom  Position  Ref  Sample01  Sample02  ...
SNP01  Chr01  33103326  G    GG        GG        ...
SNP02  Chr01  44570920  G    GT        GT        ...
SNP03  Chr02  8789570   T    CC        CC        ...
SNP04  Chr03  436874    T    TC        TT        ...
SNP05  Chr04  13430686  A    AG        AG        ...
SNP06  Chr15  14366854  G    GA        GA        ...
SNP07  Chr16  10504815  T    TC        TC        ...
SNP08  Chr17  3911259   G    GA        AA        ...
SNP09  Chr18  3552443   C    CT        CC        ...
SNP10  Chr19  10936856  C    CA        AA        ...
...    ...    ...       ...  ...       ...       ...
```
- column1: SNP ID
- column2: chromosome name
- column3: chromosome position
- column4: the base located in the second and third columns of the reference genome
- column5-n: sample name

# 2. Instantiate a GenoType object and convert it to the pandas.DataFrame class
```python
from pybioinformatic import GenoType

with GenoType('test/genotype/genotype.txt') as gt:
    df = gt.to_dataframe(index_col=0, sort_allele=True)
    print(df)
```
```text
       Chrom  Position Ref Sample01  ... Sample25 Sample26 Sample27 Sample28
ID                                   ...                                    
SNP01  Chr01  33103326   G       GG  ...       GT       GT       TT       TT
SNP02  Chr01  44570920   G       GT  ...       GG       GG       GT       GT
SNP03  Chr02   8789570   T       CC  ...       TT       TT       TT       TT
SNP04  Chr03    436874   T       CT  ...       CC       CC       CT       TT
SNP05  Chr04  13430686   A       AG  ...       AA       AA       AG       AG
SNP06  Chr15  14366854   G       AG  ...       GG       AG       AG       AG
SNP07  Chr16  10504815   T       CT  ...       CT       TT       CT       CT
SNP08  Chr17   3911259   G       AG  ...       GG       AG       AG       AA
SNP09  Chr18   3552443   C       CT  ...       TT       CT       CC       CT
SNP10  Chr19  10936856   C       AC  ...       AC       AC       AC       CC

[10 rows x 31 columns]
```
**Then, you can call the methods of the pandas.DataFrame class to conduct personalized analysis.**

# 3. Data preprocessing
**Many of the Genotype class methods are more robust when conducting biallelic SNP analysis. Therefore, we suggest that you process your genotype file before conducting the data analysis. Assuming that your genotype file is a multiple allele SNP, then the following code can be used to select the biallelic SNP.**
```python
from pybioinformatic import GenoType

with GenoType('test/genotype/raw_genotype.txt') as gt:
    df = gt.filter_allele(min_alleles=2, max_alleles=2)
    df.to_csv('test/genotype/biallelic_genotype.txt', sep='\t', na_rep='NA')
```
**The multiple allele SNP "SNP05" will be removed.**

# 4. Merge data by position
```python
from pybioinformatic import GenoType

with GenoType('test/genotype/genotype.txt') as gt1, \
        GenoType('test/genotype/genotype2.txt') as gt2:
    outer = gt1.merge(other=gt2, how='outer')  # union
    print(outer)
    inner = gt1.merge(other=gt2, how='inner')  # intersection
    print(inner)
```
```text
# outer
                  Chr       Pos Ref  ... Sample28 Sample29 Sample30
ID                                   ...                           
Chr01_33103326  Chr01  33103326   G  ...       TT       GG       GG
Chr01_44570920  Chr01  44570920   G  ...       GT       NA       NA
Chr02_8789570   Chr02   8789570   T  ...       TT       CC       CC
Chr03_436874    Chr03    436874   T  ...       TT       CT       TT
Chr04_13430686  Chr04  13430686   A  ...       AG       AG       AG
Chr15_14366854  Chr15  14366854   G  ...       AG       AG       AG
Chr16_10504815  Chr16  10504815   T  ...       CT       CT       CT
Chr17_3911259   Chr17   3911259   G  ...       AA       AG       AA
Chr18_3552443   Chr18   3552443   C  ...       CT       CT       CC
Chr19_10936856  Chr19  10936856   C  ...       CC       AC       AA
Chr19_25831234  Chr19  25831234   A  ...       NA       AC       AA

[11 rows x 33 columns]
```
```text
# inner
                  Chr       Pos Ref  ... Sample28 Sample29 Sample30
ID                                   ...                           
Chr01_33103326  Chr01  33103326   G  ...       TT       GG       GG
Chr02_8789570   Chr02   8789570   T  ...       TT       CC       CC
Chr03_436874    Chr03    436874   T  ...       TT       CT       TT
Chr04_13430686  Chr04  13430686   A  ...       AG       AG       AG
Chr15_14366854  Chr15  14366854   G  ...       AG       AG       AG
Chr16_10504815  Chr16  10504815   T  ...       CT       CT       CT
Chr17_3911259   Chr17   3911259   G  ...       AA       AG       AA
Chr18_3552443   Chr18   3552443   C  ...       CT       CT       CC
Chr19_10936856  Chr19  10936856   C  ...       CC       AC       AA

[9 rows x 33 columns]
```

# 5. Filter mendelian errors
```python
from pybioinformatic import GenoType

with GenoType('test/genotype/filter_mendelian_errors/raw_genotype.txt') as gt:
    pass_genotype, failed_genotype = gt.filter_mendelian_errors(
        parents=['Sample01', 'Sample02'],
        offsprings=['Sample03', 'Sample04'],
        num_processing=1
    )
    print(pass_genotype)
    print(failed_genotype)
```
```text
# pass_genotype
      ID  Chrom  Position Ref Sample01 Sample02 Sample03 Sample04
0  SNP02  Chr01  44570920   G       GT       GT       TT       GT
1  SNP05  Chr04  13430686   A       AG       AG       AG       AG
2  SNP06  Chr15  14366854   G       AG       AG       AA       AA
3  SNP07  Chr16  10504815   T       CT       CT       CT       CT
4  SNP09  Chr18   3552443   C       CT       CC       CT       CT
5  SNP10  Chr19  10936856   C       AC       AA       AA       AC
```
```text
# failed_genotype
      ID  Chrom  Position Ref Sample01 Sample02 Sample03 Sample04
0  SNP01  Chr01  33103326   G       GG       GG       GG       GT
1  SNP03  Chr02   8789570   T       CC       CC       CC       TT
2  SNP04  Chr03    436874   T       CT       TT       TT       CC
3  SNP08  Chr17   3911259   G       AG       AA       AA       GG
```

# 6. Calculate the missing rate, heterozygosity rate and minor allele frequency (MAF) of SNP
```python
from pybioinformatic import GenoType

with GenoType('test/genotype/genotype.txt') as gt:
    df = gt.parallel_stat_MHM(num_processing=1)
    print(df)
```
```text
      ID  Chrom  Position Ref  MissRate(%)  HetRate(%)       MAF
0  SNP01  Chr01  33103326   G          0.0   50.000000  0.464286
1  SNP02  Chr01  44570920   G          0.0   28.571429  0.428571
2  SNP03  Chr02   8789570   T          0.0   46.428571  0.410714
3  SNP04  Chr03    436874   T          0.0   35.714286  0.464286
4  SNP05  Chr04  13430686   A          0.0   82.142857  0.446429
5  SNP06  Chr15  14366854   G          0.0   60.714286  0.482143
6  SNP07  Chr16  10504815   T          0.0   60.714286  0.482143
7  SNP08  Chr17   3911259   G          0.0   42.857143  0.464286
8  SNP09  Chr18   3552443   C          0.0   46.428571  0.410714
9  SNP10  Chr19  10936856   C          0.0   57.142857  0.428571
```

# 7. Identify the SNPs that can distinguish different sample groups
```python
from pybioinformatic import GenoType

with GenoType('test/genotype/genotype.txt') as gt:
    df = gt.filter_multi_group_genotypes(
      groups=[
        ['Sample01', 'Sample03'],  # group 1
        ['Sample05', 'Sample07']  # group 2
      ]
    )
    print(df)
```
```text
               Chrom  Position Ref Sample01  ... Sample25 Sample26 Sample27 Sample28
ID                                           ...                                    
Chr02_8789570  Chr02   8789570   T       CC  ...       TT       TT       TT       TT

[1 rows x 31 columns]
```

# 8. Find the smallest set of SNP loci that can distinguish all samples by greedy forward choice
```python
from pybioinformatic import GenoType

with GenoType('test/genotype/genotype.txt') as gt:
    df = gt.find_min_distinguishing_snps()
    print(df)
```
```text
       Chrom  Position Ref Sample01  ... Sample25 Sample26 Sample27 Sample28
ID                                   ...                                    
SNP02  Chr01  44570920   G       GT  ...       GG       GG       GT       GT
SNP03  Chr02   8789570   T       CC  ...       TT       TT       TT       TT
SNP04  Chr03    436874   T       CT  ...       CC       CC       CT       TT
SNP08  Chr17   3911259   G       AG  ...       GG       AG       AG       AA
SNP09  Chr18   3552443   C       CT  ...       TT       CT       CC       CT

[5 rows x 31 columns]
```
**These 5 SNP genotype combinations can distinguish among the 28 samples.**

# 9. Calculate genotype similarity
## 9.1. Compare with the genotype matrix of the database. It can be used for cultivar identification
```python
from pybioinformatic import GenoType

with GenoType('test/genotype/genotype.txt') as gt1, \
        GenoType('test/genotype/database_genotype.txt') as gt2:  # gt2 is database
  gt1.compare(
    other=gt2,
    cmap='crest',  # color map of heatmap, 
    # more colors see https://matplotlib.org/stable/users/explain/colors/colormaps.html
    output_path='test/genotype/gs'
  )
```
**Compare using the common SNP set from the two genotype matrices. The common set of SNP is determined by the SNP ID. [`Compare.GT.xls`](https://github.com/wenlinXu-njfu/Leo/blob/main/test/genotype/gs/Compare.GT.xls) lists the SNPs which are used to compare with the database samples.**
```text
# Compare.GT.xls
ID     Ref  Sample01  Sample02  Sample03  Sample04  Sample05  ......
SNP01  G    GG        GG        GG        GT        TT        ......
SNP02  G    GT        GT        TT        GT        TT        ......
SNP03  T    CC        CC        CC        TT        CT        ......
SNP04  T    CT        TT        TT        CC        CC        ......
SNP05  A    NA        AA        AG        AG        AG        ......
SNP06  G    AG        AG        AA        AA        AG        ......
SNP07  T    CT        CT        CT        CT        TT        ......
SNP08  G    AG        AA        AA        GG        GG        ......
SNP09  C    CT        CC        CT        CT        CT        ......
SNP10  C    AC        AA        AA        AC        CC        ......
```
**Output the comparison results of genotypes in two formats.**<br />
1. **[`Sample.consistency.fmt1.xls`](https://github.com/wenlinXu-njfu/Leo/blob/main/test/genotype/gs/Sample.consistency.fmt1.xls) is a pairwise comparison between the query sample and the database samples, presenting detailed results of genotype similarity.**
```text
# Sample.consistency.fmt1.xls
Sample1   Sample2   IdenticalCount  NaCount  TotalCount  GS(%)
Sample01  Sample01  9               1        9           100.00
Sample01  Sample27  6               1        9           66.67
Sample01  Sample07  6               0        10          60.00
Sample01  Sample08  6               0        10          60.00
Sample01  Sample15  6               0        10          60.00
Sample01  Sample22  6               0        10          60.00
Sample01  Sample02  5               0        10          50.00
Sample01  Sample03  5               0        10          50.00
Sample01  Sample04  5               0        10          50.00
......    ......    ......          ......   ......      ......
```
- column1 and column2: Query and database samples for comparison<a id="link2"></a>
- column3: The number of SNPs with the same genotype between the two compared samples after excluding SNPs with missing genotype
- column4: The number of SNPs with the missing genotype between the two compared samples
- column5: The total number of SNPs without missing genotype that used for comparing genotypes
- column6: The ratio of the column3 to the column5 (genetic similarity)
2. **[`Sample.consistency.fmt2.xls`](https://github.com/wenlinXu-njfu/Leo/blob/main/test/genotype/gs/Sample.consistency.fmt2.xls) is the genotype similarity matrix**
```text
# Sample.consistency.fmt2.xls
          Sample01  Sample02  Sample03  Sample04  Sample05  ......
Sample01                                          
Sample02  50.00                                   
Sample03  50.00     70.00                         
Sample04  50.00     30.00     40.00               
Sample05  30.00     20.00     30.00     40.00     
Sample06  20.00     30.00     30.00     30.00     30.00
Sample07  60.00     60.00     40.00     30.00     30.00     ......
Sample08  60.00     20.00     50.00     60.00     40.00     ......
Sample09  22.22     22.22     11.11     11.11     33.33     ......
......    ......    ......    ......    ......    ......    ......
```
**Also outputs the visualized results of genotype similarity**
![image](https://github.com/wenlinXu-njfu/Leo/blob/main/test/genotype/gs/GS_heatmap.png)
![image](https://github.com/wenlinXu-njfu/Leo/blob/main/test/genotype/gs/GS_cluster_heatmap.png)
## 9.2. Compare the genotype consistency among different replicate samples
```python
from pybioinformatic import GenoType

with GenoType('test/genotype/genotype.txt') as gt1, \
  GenoType('test/genotype/genotype_repeat2.txt') as gt2:
  gt1.self_compare(
      other=gt2,
      output_path='test/genotype/gs'
  )
```
**Output two files: [`Sample.consistency.xls`](https://github.com/wenlinXu-njfu/Leo/blob/main/test/genotype/gs/Sample.consistency.xls) and [`Interval.stat.xls`](https://github.com/wenlinXu-njfu/Leo/blob/main/test/genotype/gs/Interval.stat.xls)**
```text
# Sample.consistency.xls
SampleName  IdenticalCount  NaCount  TotalCount  GS(%)
Sample01    9               1        9           100.0
Sample02    9               0        10          90.0
Sample03    10              0        10          100.0
Sample04    10              0        10          100.0
Sample05    10              0        10          100.0
Sample06    8               2        8           100.0
Sample07    10              0        10          100.0
Sample08    9               1        9           100.0
Sample09    10              0        10          100.0
```
**[`Sample.consistency.xls`](https://github.com/wenlinXu-njfu/Leo/blob/main/test/genotype/gs/Sample.consistency.xls) format instructions can be reference in the previous [section](#link2).**<br />
```text
# Interval.stat.xls
[0, 5]     0.0
(5, 10]    0.0
(10, 15]   0.0
(15, 20]   0.0
(20, 25]   0.0
(25, 30]   0.0
(30, 35]   0.0
(35, 40]   0.0
(40, 45]   0.0
(45, 50]   0.0
(50, 55]   0.0
(55, 60]   0.0
(60, 65]   0.0
(65, 70]   0.0
(70, 75]   0.0
(75, 80]   0.0
(80, 85]   0.0
(85, 90]   1.0
(90, 95]   0.0
(95, 100]  27.0
total      28.0
```
**[`Interval.stat.xls`](https://github.com/wenlinXu-njfu/Leo/blob/main/test/genotype/gs/Interval.stat.xls) count the number of samples in each similarity interval.**

# 10. Draw genotype heatmap
```python
from pybioinformatic import GenoType

with GenoType('test/genotype/genotype.txt') as gt:
    gt.draw_genotype_heatmap(
        genotype_colors=None,  # use default color ['#2c7bb6', '#ffffbf', '#d7191c']
        figure_size=(10, 8),
        font_size=10,
        out_file='test/genotype/genotype_heatmap/genotype_heatmap_with_default_colors.png'
    )
```
![image](https://github.com/wenlinXu-njfu/Leo/blob/main/test/genotype/genotype_heatmap/genotype_heatmap_with_default_colors.png)
**Or use customized genotype colors**
```python
from pybioinformatic import GenoType

with GenoType('test/genotype/genotype.txt') as gt:
    gt.draw_genotype_heatmap(
        genotype_colors=['#AAC026', '#F09B39', '#F24532'],
        figure_size=(10, 8),
        font_size=10,
        out_file='test/genotype/genotype_heatmap/genotype_heatmap_with_customized_colors.png'
    )
```
![image](https://github.com/wenlinXu-njfu/Leo/blob/main/test/genotype/genotype_heatmap/genotype_heatmap_with_customized_colors.png)

# 11. Build finger-print
## 11.1. Select core SNPs
**Select a minimal core set of SNP loci that can fully discriminate all samples using a genetic algorithm (GA) with optional uniform chromosomal distribution constraints.**<br />
![image](https://github.com/wenlinXu-njfu/Leo/blob/main/test/genotype/finger_print/genetic_algorithm_flowchart.png)
```python
from pybioinformatic import GenoType

with GenoType('test/genotype/finger_print/candidate_SNP.GT.txt') as gt:
    df = gt.select_core_snps_ga(
        min_per_chrom=6,
        pop_size=100,
        max_gen=200,
        cross_prob=0.8,
        mut_prob=0.05,
        elite_count=2,
        seed=42,
        no_elimination=False,
        output_curve=True,
        curve_prefix='test/genotype/finger_print/identification_rate_curve'
    )
    df.to_csv('test/genotype/finger_print/coreSNP.txt', sep='\t', index=False)
```
```text
# coreSNP.txt
ID       Chrom  Position  Ref  sample_01  sample_02  sample_03  sample_04  ...
SNP0163  Chr01  8039838   T    TT         CC         CC         CC         ...
SNP0346  Chr01  17283882  C    CC         CC         TT         CT         ...
SNP0503  Chr01  25819482  C    TT         TT         CC         CT         ...
SNP0688  Chr01  34640052  C    GG         CC         CG         CC         ...
SNP0907  Chr01  43266259  C    TT         TT         TT         CT         ...
SNP1120  Chr01  51505719  C    TT         CC         CC         CT         ...
SNP1159  Chr02  4408344   T    AT         AA         TT         TT         ...
SNP1226  Chr02  7850760   C    TT         TT         TT         TT         ...
SNP1291  Chr02  12057445  A    CC         CC         CC         AA         ...
...      ...    ...       ...  ...        ...        ...        ...        ...
```
![image](https://github.com/wenlinXu-njfu/Leo/blob/main/test/genotype/finger_print/identification_rate_curve.png)<br />
**More detail information use `help(gt.select_core_snps_ga)`**
## 11.2. Visualized fingerprint
```python
from pybioinformatic import GenoType

with GenoType('test/genotype/finger_print/coreSNP.txt') as gt:
    gt.draw_finger_print(output_file='test/genotype/finger_print/fingerprint.png')
```
![image](https://github.com/wenlinXu-njfu/Leo/blob/main/test/genotype/finger_print/fingerprint.png)

# 12. Analysis of phenotypic differences by genotype
**Perform phenotypic association analysis (ANOVA) for the genotype of each SNP in the input genotype matrix and visualize the results.**
```python
from pybioinformatic import GenoType

with GenoType('test/genotype/genotype.txt') as gt:
    gt.genotype_phenotype_anova_analysis(
        trait_file='test/genotype/genotype_phenotype_anova/trait.txt',
        trait_name='Trait3',
        out_path='test/genotype/genotype_phenotype_anova',
        mark_sample_list=['Sample07', 'Sample20']
    )
```
![image](https://github.com/wenlinXu-njfu/Leo/blob/main/test/genotype/genotype_phenotype_anova/Chr03_436874_Trait3.png)

[Back to top](#link1)