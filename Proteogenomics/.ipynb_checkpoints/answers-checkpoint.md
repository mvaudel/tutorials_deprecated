---
layout: page
title: Answers-Proteogenomics
---

# Answers - Proteogenomics

## 1. Novel Peptides

##### [:thought_balloon:](novel_peptides.md#thought_balloon-based-on-you-knowledge-of-peptide-and-protein-identification-can-you-anticipate-challenges-posed-by-these-proteogenomic-databases) _Based on your knowledge of peptide and protein identification, can you anticipate challenges posed by these proteogenomic databases?_

The first challenge posed by proteogenomic databases is that they are very large, and therefore require much more computational power. The inflation in the collection of possible sequences, called search space, also increases the probability that two peptide sequences are similar, hence reducing the discrimination power of search engines. Thus reduced discrimination power results in a reduced identification rate at fixed error rate [(1)](#references). The estimation of error rates in proteomics usually relies on the target-decoy approach [(2)](#references), which is very reliable to track random matches, but not the errors due to partial matches [(3)](#references). Consequently, error rate estimation is very challenging in proteogenomic databases and requires spectrum-level inspection of the matches [(4)](#references). Finally, the increased number of similar sequences reduces the probability to find peptides unique to a protein, hence complexifying protein inference [(5)](#references).

##### [:thought_balloon:](novel_peptides.md#thought_balloon-what-do-the-different-columns-in-the-table-represent) _What do the different columns in the table represent?_

This document contains two tables, a list of novel peptides and a list of peptides containing single amino acid variants (SAAV).

- Novel Peptides

| Column Name | Description |
| ----------- | ----------- |
| Peptide Sequence | The amino acid sequence of the identified peptide. |
| MSGF+ SpecEva | The score (e-value) produced by the search engine (ms-gf+) when evaluating the match between the peptide and the fragmentation spectrum. The lower, the better. |
| # PSMs | The number of spectra where this peptide was matched. The more, the better. |
| Annotation - GRCh37 - hg19 | The annotation of the genetic locus coding the peptide. _GRCh37/hg19_ corresponds to the build of the genome used for the analysis. |
| Chromosome, start, end, strand, locus number | Genomic coordinates of the locus coding the peptide in the GRCh37 build. |
| Category and sequence similarity to known proteins | Quality control report on possible mismatch with other proteins. Note that this is provided at two different steps of the bioinformatic pipeline. See _Proteogenomics search and Class-specific FDR_ in the [Supplementary Information](../resources/Johansson_et_al_breast_cancer_quantitative_proteome_and_proteogenomic_landscape/supplementary_information.pdf). |
| Peptide explained by nsSNP in CanProVar 2.0 | Known variation that could explain the peptide. |
| hg38_coordinates | Genomic coordinates of the locus coding the peptide in the _GRCh38/hg38_ build. |
| (closest) matched protein | Protein best matched with this peptide. |
| Nterm.seq.3aa. - Aligned sequence - Cterm.seq.3aa. | Sequence, upstream, and downstream amino acids in the mached protein. |
| Identity | Identity score between the peptide and its matched counterpart. The higher the better. |
| Peptide length - Alignment length - # mismatches - # gaps | Summary information on the identified and matched peptides. |
| top 0.2% MHCflurry | class I peptide/MHC binding affinity prediction: see _NSearch neoantigen candidates in draft proteomes data and MHC binding prediction_ in the [Supplementary Information](../resources/Johansson_et_al_breast_cancer_quantitative_proteome_and_proteogenomic_landscape/supplementary_information.pdf). |
| Identified in Kim et al draft proteome | Indicates whether the peptide was identified in [(6)](#references). |
| Class | Class of the locus coding this peptide. |
| Associated gene (closest genes in the genome) | Nearest gene for the given locus. |
| Category | Category of the nearest gene. |
| TMT quantification | Abundance of the peptide in the different tumors screened. |
| Orthogonal data | Additional data on the peptide: see _Orthogonal validation data of proteogenomics searches_ in the [Supplementary Information](../resources/Johansson_et_al_breast_cancer_quantitative_proteome_and_proteogenomic_landscape/supplementary_information.pdf). |

- SAAV

| Column Name | Description |
| ----------- | ----------- |
| Peptide Sequence | The amino acid sequence of the identified peptide. |
| Position of amino acid substitution | The position of the SAAV on the peptide sequence. |
| Matched protein(s) | The proteins matched. |
| Peptide sequence with modifications | The peptide sequence annotated with the masses of the modifications. |
| Peptide length | The sequence length in number of amino acids. |
| Fragmentation Method | The method used to fragment the peptide during mass spectrometry analysis. |
| Precursor m/z | The mass over charge ratio of the precursor selected for fragmentation. |
| IsotopeError | The isotopic difference between the measured m/z and the m/z of the peptide. |
| Precursor Error (ppm) | The relative difference in parts per million (ppm) between the measured m/z and the m/z of the peptide after isotopic correction. |
| Charge | The charge of the precursor. |
| DeNovoScore, MSGFScore | Score of the match between the fragmentation spectrum and the peptide. The higher, the better. |
| SpecEValue, EValue | E-value of the score of the match between the fragmentation spectrum and the peptide. The lower, the better. |
| SAAV Only - Peptide FDR | FDR estimation at the given score for the SAAV peptides only. See _Proteogenomics search and Class-specific FDR_ in the [Supplementary Information](../resources/Johansson_et_al_breast_cancer_quantitative_proteome_and_proteogenomic_landscape/supplementary_information.pdf). |
| SpectraFile, TMT set, ScanNum | File, TMT set, and scan number for the spectrum where the peptide was identified. |
| SpectrumAI verification | Summary information of the spectrum-level validation of the SAAV using flanking fragment ions [(4)](#references). |
| TMT intensity | Intensity of the fragment ions from the different TMT channels. |


##### [:thought_balloon:](novel_peptides.md#thought_balloon-can-you-speculate-on-how-these-different-classes-of-loci-can-yield-novel-peptides) _Can you speculate on how these different classes of loci can yield novel peptides?_

Mutations can yield new coding regions through alteration of start, stop, and splicing information. Exons can hence be extended by upstream and downstream regions, extending canonical sequences and disrupting splicing. Alternatively, completely new peptides or proteins can be created from intergenic regions.

##### [:thought_balloon:](novel_peptides.md#thought_balloon-what-do-these-categories-represent) _What do these categories represent?_

| Category | # Peptides | Description |
| -------- | ---------- | ----------- |
| pseudogene | 93 | Pseudogenes are segments of DNA related to real genes but have lost functionality. |
| 5UTR | 43 | 5′ untranslated region (5′ UTR) is directly upstream of the translated region. |
| intronic | 33 | Introns are nucleotide sequence removed by RNA splicing. |
| exonic.Alt.ORF | 18 | Exons from alternative open reading frames. |
| exon_extension | 14 | Extension of exon. |
| ncRNA | 12 | Non-coding RNA |
| intergenic | 2 | Intergenic section |
| 3UTR | 1 | 3′ untranslated region (3′ UTR) is directly downstream of the translated region. |


## 2. Variation Analysis

##### [:thought_balloon:](variation_analysis.md#thought_balloon-if-we-assume-a-linear-relationship-between-number-of-alleles-and-peptide-abundance-what-should-be-the-peptide-distribution-for-each-genotype) _If we assume a linear relationship between number of alleles and peptide abundance, what should be the peptide distribution for each genotype?_

![peptide_per_allele](resources/images/peptideAllelePlot.png?raw=true "Peptide abundance vs genotype")

Under this hypothesis: (1) if the genotype is homozygous for the reference allele, _AA_, one would expect only peptides carrying the reference amino acid; (2) if the genotype is heterozygous, _AB_, one would expect a 50-50 distribution between peptides carrying the reference and alternative amino acids; (3) if the genotype is homozygous for the alternative allele, _BB_, one would expect only peptides carrying the alternative amino acid.

##### [:thought_balloon:](answers.md#thought_balloon-how-do-we-need-to-transform-the-tables-to-compare-saav-peptide-level-intensities-to-gene-level-intensities) _How do we need to transform the tables to compare SAAV peptide-level intensities to gene-level intensities?_

The data in the SAAV peptides table are not normalized. Thus, we need to scale all intensities using the reference channel (131). Each line represents a PSM, peptides might therefore be present multiple times, we need to aggregate these different measurements per peptide in order to have a peptide-level table. Then, the different TMT batches need to be aligned and linked to the different tumors. Finally, the gene-level intensities need to be extracted for each peptide. Note that one peptide might map to multiple genes.

##### [:thought_balloon:](answers.md#thought_balloon-why-can-there-be-multiple-peptide-per-gene-and-gene-per-peptide-is-it-correct-to-represent-peptides-by-their-sequence) _Why can there be multiple peptide per gene, and gene per peptide? Is it correct to represent peptides by their sequence?_

A gene can contain multiple peptides with different SAAV. Conversely, a peptide carrying an SAAV can be shared between different genes. Consequently genes and peptides can appear multiple times in the table and their quantification values as well.

If a protein carries a post-translational modification (PTM), its abundance and function will differ from the non-modified counterpart. Therefore, modifications ough to be considered when comparing peptides. In the data set of this tutorial, modifications were not accounted for.

##### [:thought_balloon:](answers.md#thought_balloon-why-are-all-intensities-at-the-bottom-how-can-we-transform-the-data-to-better-visualize-these-distributions) _Why are all intensities at the bottom? How can we transform the data to better visualize these distributions?_

On a proteome-wide scale, measured peptide and protein abundances generally distribute log-normally. This means that they span several orders of magnitude, with most abundances are low, with a long tail towards high abundances, and a short tail towards lower abundances. Note that the detectability and quantification performance strongly influence this distribution, and it can be that the distribution of protein abundances in the sample differ from what we measure.

In order to visualize such distributions, you can transform the intensities logarithmically, or use a log-scale on the x-axis, then their distribution should look bell-shaped and symmetrical.

![abundance_distribution](resources/images/peptideProteinDistribution.png?raw=true "Peptide and Protein abundance distributions")

##### [:thought_balloon:](answers.md#thought_balloon-why-does-the-curve-have-this-shape) _Why does the curve have this shape?_

The model used for normalization is based on a log-Normal distribution, see the `family` argument provided to _GAMLSS_:
```
model <- gamlss(
        formula = as.formula(paste0(column, " ~ x")),
        family = LOGNO,
        data = trainingDF,
        
    )
```

This accounts for the log-normal distribution of the intensities observed in the previous question, and explains why the curve has this logarithmic shape. Using a log-scale on the x axis therefore aligns the points.

![zScore_log](resources/images/zLogIntensity.png?raw=true "Z-score vs log")


## 3. CNA-protein

##### [:thought_balloon:](cna-protein.md#thought_balloon-what-do-the-columns-represent-what-is-the-difference-between-pearson-and-spearman-correlations) What do the columns represent? What is the difference between Pearson and Spearman correlations?

| Column Name | Description |
| ----------- | ----------- |
| gene, chromosome, start, end, band, position | The chromosomic coordinates of the CNA. |
| mRNA_ANOVA_AdjPval, protein_ANOVA_AdjPval | The significance that the RNA or protein, respectively, are differentially abundant between the tumors using an ANOVA analysis. [-log10] |
| mRNA_Pearson_correlation, protein_Pearson_correlation | The Pearson correlation coefficient (R) for the RNA or the protein abundance, respectively, with the CNA. |
| mRNA_Pearson_pval, protein_Pearson_pval | The Pearson correlation significance for the RNA or the protein abundance, respectively, with the CNA. |
| mRNA_Spearman_correlation, protein_Spearman_correlation | The Spearman correlation coefficient (R) for the RNA or the protein abundance, respectively, with the CNA. |
| mRNA_Spearman_pval, protein_Spearman_pval | The Spearman correlation significance for the RNA or the protein abundance, respectively, with the CNA. |

The Spearman correlation is a Pearson correlation of the ranks of the values. Working on the ranks is more robust, especially towards outliers. See the R help on _cor_ for more details: `?cor`.

##### [:thought_balloon:](cna-protein.md#thought_balloon-how-many-gaussian-distributions-were-suggested-by-the-model-what-do-the-mixing-probabilities-means-and-variances-represent) What do the _Mixing probabilities_, _Means_, and _Variances_ represent?

The Model suggested using two gaussian distributions. This is based on the [Bayesian Information Criterion (BIC)](https://en.wikipedia.org/wiki/Bayesian_information_criterion), which evaluates the quality of the fitting relatively to the model complexity: an overly complex model will fit very well training data but is likely to overfit, while an overly simple model is unlikely to overfit but is likely to underperform. Optimizing the BIC allows finding a trade-off between the two. Below, the BIC is plotted against the number of components.

![gmm_n_bic](resources/images/gmm_n.png?raw=true "BIC vs. n")

For each Gaussian distribution, mclust returns the following paramters:
- Mixing probabilities: the contribution of the distribution to the overall density. Their sum should be 1.
- Mean: the mean of the distribution.
- Variance: the variance of the distribution.

##### [:thought_balloon:](cna-protein.md#thought_balloon-based-on-this-how-many-cnas-are-considered-attenuated) Based on this, how many CNAs are considered attenuated?

The mixing probability of the distribution with lowest mean (_i.e._ most attenuated) is approx 45.6%, considering 9533 CNAs this makes 4343 CNAs attenuated.

##### [:thought_balloon:](cna-protein.md#thought_balloon-which-component-represents-the-attenuated-distribution-how-can-we-classify-cnas-based-on-these-distributions) Which component represents the attenuated distribution? How can we classify CNAs based on these distributions?

The distribution to the left, i.e. lowest mean, is the one representing the attenuated population. 

It is possible to use the ratio of a distribution to the sum to estimate the probability at a given attenuation score to belong from either category. Alternatively, it is possible to use the cumulative distribution function to evaluate the appartenance to a population.

##### [:thought_balloon:](cna-protein.md#thought_balloon-using-this-threshold-what-is-the-share-of-cnas-considered-attenuated-that-would-come-from-the-component-2-distribution-what-is-the-share-of-cnas-considered-not-attenuated-that-would-come-from-the-component-1-how-will-this-influence-the-analyses) Using this threshold: What is the share of CNAs considered attenuated that would come from the component 2 distribution? What is the share of CNAs considered not attenuated that would come from the component 1? How will this influence the analyses?

| Category | Component 1 | Component 2 |
| -------- | ----------- | ----------- |
| Attenuated | 1383 | 255 |
| Not attenuated | 3049 | 4846 |

Approximately 16% of the CNAs considered attenuated are derived from the second component. Approximately 39% of the CNAs considered not attenuated are derived from the first component. This means that our classification has many false positives/negatives, which will reduce our ability to find association with biological features.

##### [:thought_balloon:](cna-protein.md#thought_balloon-what-is-the-number-of-proteins-per-gene-reported-in-the-cna-table-how-is-this-going-to-influence-the-analysis) What is the number of proteins per gene reported in the CNA table? How is this going to influence the analysis?

According to the UniProt mapping, the number of proteins per gene ranges from 1 to approx. 8000. The genes mapping to most proteins are listed in the table below.

| Gene | # accessions mapped |
| ---- | ------------------- |
| HLA-B | 7983 |
| HLA-A | 6874 |
| HLA-C | 6516 |
| HLA-DRB1 | 3095 |
| HLA-DQB1 | 2145 |
| HLA-DPB1 | 1761 |
| HLA-DQA1 | 207 |
| ATXN3 | 204 |
| HLA-DPA1 | 127 |
| BRCA1 | 125 |
| TP53 | 120 |
| APOL1 | 107 |

The result for each gene can be obtained using the following script in R.

```
geneOccurenceDF <- as.data.frame(table(accessionsMapping$gene))
occurrenceFrequencyDF <- as.data.frame(table(geneOccurenceDF[[2]]))
```

When mapping the CNA results to functional databases, these genes are going to map everywhere in the proteome, reducing our ability to find meaningful information.

## References

(1) [Anatomy and evolution of database search engines-a central component of mass spectrometry based proteomic workflows](https://www.ncbi.nlm.nih.gov/pubmed/28902424)

(2) [Target‐decoy search strategy for increased confidence in large‐scale protein identifications by mass spectrometyr](https://www.ncbi.nlm.nih.gov/pubmed/17327847)

(3) [Analysis of the resolution limitations of peptide identification algorithms](https://www.ncbi.nlm.nih.gov/pubmed/21995378)

(4) [Discovery of coding regions in the human genome by integrated proteogenomics analysis workflow](https://www.ncbi.nlm.nih.gov/pubmed/29500430)

(5) [Interpretation of shotgun proteomic data: the protein inference problem](https://www.ncbi.nlm.nih.gov/pubmed/16009968)

(6) [A draft map of the human proteome](https://www.ncbi.nlm.nih.gov/pubmed/24870542)
