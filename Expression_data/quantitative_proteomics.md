Analysis of Quantitative Proteomics Data
================

The result of a proteomic quantitative analysis is a list of peptide and protein abundances for every protein in different samples, or abundance ratios between the samples. The downstream interpretation of these data varies according to the experimental design. In this chapter we will describe different generic methods for the interpretation of quantitative datasets.

The dataset used here for illustrative purposes is freely available through the ProteomeXchange [(1)](#references) consortium via the PRIDE [(2)](#references) partner repository under the accession number PXD000441, and available in the *resources* folder under *Results\_MQ\_5cell-line-mix*. It consists of five cell lines derived from acute myeloid leukemia (AML) patients measured with a spiked-in internal stanard (IS) obtained from the combination of the same five AML cell lines that have been metabolically labeled with heavy isotopes [(3)](#references) and analyzed using MaxQuant [(4)](#references) version 1.4.1.2. See [(5)](#references) for details.

This chapter introduces the basic methods, and is by no means aiming at covering all possibilities or present a reference workflow. Please continue exploring the data on your own and critically adapt the interpretation workflow to your experiment - there is no one workflow fits all!

Installation
------------

This tutorial uses [R](https://www.r-project.org), a language that allows the simple manipulation of large datasets. We will use R from the open source [RSudio](www.rstudio.com) environment. Please make sure to have RStudio installed on your computer.

File import
-----------

The MaxQuant report files are tab separated text files that can readily be imported in R as data frame. A data frame allows the convenient manipulation of large tables.

``` r
proteinGroupsInput <- read.table(file = "proteinGroups_5cell-line-mix.txt", header = T, stringsAsFactors = F, sep = "\t")
```

You should see the proteinGroupsInput in your *Environment* panel. You can get the size of the table using the *length* function.

``` r
nColumns <- length(proteinGroupsInput)
nLines <- length(proteinGroupsInput$Protein.IDs)
paste("Number of columns: ", nColumns, ", number of protein groups: ", nLines, sep="")
```

    ## [1] "Number of columns: 191, number of protein groups: 3982"

Filtering of the protein groups
-------------------------------

MaxQuant indicates contaminants, decoy sequences, and proteins only identified by site by a '+' in their columns. You can access the number of every category by using the *table* function.

``` r
nContaminants <- table(proteinGroupsInput$Contaminant)
paste("Number of Contaminants: ", nContaminants[2], ", Others: ", nContaminants[1], sep="")
```

    ## [1] "Number of Contaminants: 104, Others: 3878"

``` r
nDecoys <- table(proteinGroupsInput$Reverse)
paste("Number of decoys: ", nDecoys[2], ", Others: ", nDecoys[1], sep="")
```

    ## [1] "Number of decoys: 76, Others: 3906"

``` r
nOnlyIdentifiedBySite <- table(proteinGroupsInput$Only.identified.by.site)
paste("Number of only identified by site: ", nOnlyIdentifiedBySite[2], ", Others: ", nOnlyIdentifiedBySite[1], sep="")
```

    ## [1] "Number of only identified by site: 109, Others: 3873"

Filter these lines and save the result in a new table.

``` r
proteinGroups <- proteinGroupsInput[proteinGroupsInput$Contaminant != '+' & proteinGroupsInput$Reverse != '+' & proteinGroupsInput$Only.identified.by.site != '+',]
nLines <- length(proteinGroups$Protein.IDs)
nFiltered <- length(proteinGroupsInput$Protein.IDs) - length(proteinGroups$Protein.IDs)
paste("New number of protein groups: ", nLines, ", Number of protein groups removed: ", nFiltered, sep="")
```

    ## [1] "New number of protein groups: 3732, Number of protein groups removed: 250"

References
----------

1.  [Vizcaino, J.A. et al., *ProteomeXchange provides globally coordinated proteomics data submission and dissemination*, Nature Biotechnology, 2014](https://www.ncbi.nlm.nih.gov/pubmed/24727771)
2.  [Martens, L. et al., *PRIDE: the proteomics identifications database*, Proteomics, 2005](https://www.ncbi.nlm.nih.gov/pubmed/16041671)
3.  [Geiger, T. et al., *Super-SILAC mix for quantitative proteomics of human tumor tissue*, Nature Methods, 2010](https://www.ncbi.nlm.nih.gov/pubmed/20364148)
4.  [Cox, J. and Mann, M, *MaxQuant enables high peptide identification rates, individualized p.p.b. range mass accuracies and proteome-wide protein quantification*, Nature Biotechnology, 2008](https://www.ncbi.nlm.nih.gov/pubmed/19029910)

GitHub Documents
----------------

This is an R Markdown format used for publishing markdown documents to GitHub. When you click the **Knit** button all R code chunks are run and a markdown file (.md) suitable for publishing to GitHub is generated.

Including Code
--------------

You can include R code in the document as follows:

``` r
summary(cars)
```

    ##      speed           dist       
    ##  Min.   : 4.0   Min.   :  2.00  
    ##  1st Qu.:12.0   1st Qu.: 26.00  
    ##  Median :15.0   Median : 36.00  
    ##  Mean   :15.4   Mean   : 42.98  
    ##  3rd Qu.:19.0   3rd Qu.: 56.00  
    ##  Max.   :25.0   Max.   :120.00

Including Plots
---------------

You can also embed plots, for example:

![](quantitative_proteomics_files/figure-markdown_github/pressure-1.png)

Note that the `echo = FALSE` parameter was added to the code chunk to prevent printing of the R code that generated the plot.
