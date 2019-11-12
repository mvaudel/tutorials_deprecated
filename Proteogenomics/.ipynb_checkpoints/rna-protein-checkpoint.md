RNA-Protein
================
Marc Vaudel
2019-10-28

# RNA-Protein

In the [*novel peptides*](novel_peptides.md) and [*variation
analysis*](variation_analysis.md) chapters, novel and variant peptides
were identified by Johansson *et al.* [(1)](#references) using a
proteogenomic search strategy, where genomic data are used to inform the
proteomics search strategy. It is possible to use RNA sequencing data
complementarily or instead of the genetic sequence to identify
non-canonical protein products.

##### [:thought\_balloon:](answers.md#thought_balloon-what-are-the-advantages-and-shortcomings-of-using-transcript-sequences-instead-of-or-in-addition-to-genomic-data) *What are the advantages and shortcomings of using transcript sequences instead of or in addition to genomic data?*

## Libraries

We will need the following libraries, please make sure that they are
installed.

``` r
library(conflicted)
library(tidyr)
library(dplyr)
library(ggplot2)
library(scico)
library(gamlss)
library(igraph)

theme_set(theme_bw(base_size = 11))

conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
```

## Transcript - protein abundance relationship

In addition to the identification of non-canonical genetic products,
studying the relationship between transcript and protein abundance
levels can inform on biological processes ongoing in the studied
samples. Johansson *et al.* [(1)](#references) provide the correlation
between transcript and protein abundance levels across tumors per gene
in [Supplementary Table
1](../resources/Johansson_et_al_breast_cancer_quantitative_proteome_and_proteogenomic_landscape).
The table was extracted to an R-friendly text format for this tutorial,
and is available in
[resources/data/rna-protein.gz](resources/data/rna-protein.gz). Like in
the [*variation analysis chapter*](variation_analysis.md), we will
generate random correlations to compare these data to.

##### :pencil2: Load the data in R as in the code below, Z-transform the protein count, and plot a histogram of the correlation.

``` r
rnaProteinCorDF <- read.table(
    file = "resources/data/rna-protein.gz",
    header = T,
    sep = "\t",
    comment.char = "",
    quote = "",
    stringsAsFactors = F
) %>%
    filter(
        !is.na(mRNA_protein_correlation)
    )


# Generate random correlations

rnaProteinCorDF$random_correlation <- NA

for (i in 1:nrow(rnaProteinCorDF)) {
    
    x <- rnorm(45)
    y <- rnorm(45)
    
    rnaProteinCorDF$random_correlation[i] <- cor(
        x = x,
        y = y,
        method = "spearman"
    )
    
}

# Normalize the protein abundance, set missing values to the median

trainingDF <- rnaProteinCorDF %>%
    filter(
        !is.na(protein_copy_number)
    )

trainingDF$x <- 0

model <- gamlss(
    formula = as.formula("protein_copy_number ~ x"),
    family = LOGNO,
    data = trainingDF,
    
)
```

    ## GAMLSS-RS iteration 1: Global Deviance = 295213.6 
    ## GAMLSS-RS iteration 2: Global Deviance = 295213.6

``` r
trainingDF$z_protein_copy_number <- centiles.pred(
    obj = model, 
    xname = "x", 
    xvalues = trainingDF$x, 
    yval = trainingDF$protein_copy_number, 
    type = "z-scores"
)

rnaProteinCorDF <- trainingDF %>%
    select(
        gene, z_protein_copy_number
    ) %>%
    right_join(
        rnaProteinCorDF,
        by = "gene"
    ) %>%
    mutate(
        z_protein_copy_number = ifelse(is.na(z_protein_copy_number), 0, z_protein_copy_number)
    )


# Transform the data frame from wide to long format

rnaProteinCorLongDF <- rnaProteinCorDF %>%
    gather(
        mRNA_protein_correlation, random_correlation,
        key = "data",
        value = "correlation"
    ) %>%
    mutate(
        data = ifelse(data == "mRNA_protein_correlation", "Data", "Random")
    ) %>%
    arrange(
        desc(data)
    )


# Build plot

ggplot(
    data = rnaProteinCorLongDF
) +
    geom_histogram(
        mapping = aes(
            x = correlation,
            fill = data
        ),
        alpha = 0.8,
        bins = 20,
        position = "dodge"
    ) +
    scale_x_continuous(
        name = "Correlation [Spearman]"
    ) +
    scale_y_continuous(
        name = "# SAAV Peptide - Protein"
    ) +
    scale_fill_manual(
        values = c("darkblue", "black")
    ) +
    theme(
        legend.position = "top",
        legend.title = element_blank()
    )
```

![](rna-protein_files/figure-gfm/import-1.png)<!-- -->

##### [:thought\_balloon:](answers.md#thought_balloon-what-are-the-advantages-and-shortcomings-of-using-transcript-sequences-instead-of-or-in-addition-to-genomic-data) *How does the correlation evolve relatively to the abundance of the protein? Should it be used as covariate? What other protein characteristics could influence the correlation?*

##### :pencil2: Compute Z-scores of the correlation using protein abundance as covariate, plot one against the other.

``` r
# Normalize the correlation

rnaProteinCorDF$x <- 0

model <- gamlss(
    formula = as.formula("mRNA_protein_correlation ~ x"),
    family = NO,
    data = rnaProteinCorDF[, c("mRNA_protein_correlation", "x")]
)
```

    ## GAMLSS-RS iteration 1: Global Deviance = -602.3349 
    ## GAMLSS-RS iteration 2: Global Deviance = -602.3349

``` r
rnaProteinCorDF$z_mRNA_protein_correlation <- centiles.pred(
    obj = model, 
    xname = "x", 
    xvalues = rnaProteinCorDF$x, 
    yval = rnaProteinCorDF$mRNA_protein_correlation, 
    type = "z-scores"
)

model <- gamlss(
    formula = as.formula("random_correlation ~ x"),
    family = NO,
    data = rnaProteinCorDF[, c("random_correlation", "x")]
)
```

    ## GAMLSS-RS iteration 1: Global Deviance = -9177.734 
    ## GAMLSS-RS iteration 2: Global Deviance = -9177.734

``` r
rnaProteinCorDF$z_random_correlation <- centiles.pred(
    obj = model, 
    xname = "x", 
    xvalues = rnaProteinCorDF$x, 
    yval = rnaProteinCorDF$random_correlation, 
    type = "z-scores"
)


# Transform the data frame from wide to long format

z_rnaProteinCorLongDF <- rnaProteinCorDF %>%
    select(
        gene, z_mRNA_protein_correlation, z_random_correlation
    ) %>%
    gather(
        z_mRNA_protein_correlation, z_random_correlation,
        key = "data",
        value = "z_correlation"
    ) %>%
    mutate(
        data = ifelse(data == "z_mRNA_protein_correlation", "Data", "Random")
    ) %>%
    arrange(
        desc(data)
    )

rnaProteinCorLongDF <- rnaProteinCorLongDF %>%
    left_join(
        z_rnaProteinCorLongDF,
        by = c("gene", "data")
    )


# Build plot

ggplot(
    data = rnaProteinCorLongDF
) +
    geom_point(
        mapping = aes(
            x = correlation,
            y = z_correlation,
            col = data
        ),
        alpha = 0.1
    ) +
    scale_x_continuous(
        name = "Correlation"
    ) +
    scale_y_continuous(
        name = "Correlation [Z-score]"
    ) +
    scale_color_manual(
        values = c("darkblue", "black")
    ) +
    theme(
        legend.position = "top",
        legend.title = element_blank()
    )
```

![](rna-protein_files/figure-gfm/normalization-1.png)<!-- -->

## Reduced correlation through participation in complexes and biochemical reactions

Johansson *et al.* [(1)](#references) show that the correlation between
transcript and protein abundances is influenced by the involvement of
proteins in complexes and protein-protein interactions (PPIs). Like in
the [*CNA-protein*](cna_protein.md) chapter, we are now going to
investigate the relationship between RNA-protein correlation and protein
participation in complexes, biochemical reactions, and PPIs.

##### :pencil2: Load the identifiers mapping as done in the [*CNA-protein*](cna_protein.md) chapter.

``` r
# Uniprot accession mapping

accessionsMapping <- read.table(
    file = "resources/data/HUMAN_9606_idmapping.dat.gz", 
    header = F, 
    sep = "\t", 
    quote = "", 
    comment.char = "", 
    stringsAsFactors = F
)
names(accessionsMapping) <- c("accession", "source", "gene")
accessionsMapping <- accessionsMapping %>%
    filter(
        source == "Gene_Name"
    ) %>%
    select (
        gene, accession
    ) %>%
    filter(
        gene %in% rnaProteinCorDF$gene
    ) %>%
    distinct()

# Exclude gene mapping with more than 10 proteins per gene

geneOccurenceDF <- as.data.frame(
    x = table(accessionsMapping$gene), 
    stringsAsFactors = F
)
names(geneOccurenceDF) <- c("gene", "nProteins")

excludedGenes <- geneOccurenceDF$gene[geneOccurenceDF$nProteins >= 100]
```

##### :pencil2: Load complex data, plot the RNA-protein correlation against the number of complexes a protein is involved in.

``` r
# Complexes

complexesDF <- read.table(
    file = "resources/data/complexes.03.11.19.tsv.gz", 
    header = T, 
    sep = "\t", 
    stringsAsFactors = F, 
    quote = "", 
    comment.char = ""
)


# Extract the protein to complex link, extract protein pairs involved in the same complex

complexesParticipantsList <- list()
complexesEdgesList <- list()

for (i in 1:nrow(complexesDF)) {
    
    complex <- complexesDF$complex_accession[i]
    participants <- complexesDF$participants_stoichiometry[i]
    
    proteinsSplit <- strsplit(
        x = participants,
        split = "\\|"
    )[[1]]
    participantsList <- strsplit(
        x = proteinsSplit,
        split = "[\\(\\)]"
    )
    
    participantsDF <- as.data.frame(
        t(
            unname(
                as.data.frame(
                    participantsList,
                    stringsAsFactors = F
                )
            )
        ),
        stringsAsFactors = F
    )
    names(participantsDF) <- c("protein", "stoichiometry")
    participantsDF$complex = complex
    
    complexesParticipantsList[[i]] <- participantsDF
    
    edgesDF <- data.frame(
        from = rep(participantsDF$protein, times = nrow(participantsDF)),
        to = rep(participantsDF$protein, each = nrow(participantsDF)),
        stringsAsFactors = F
    ) %>%
        filter(
            from != to
        ) %>%
        mutate(
            complex = complex
        )
    
    complexesEdgesList[[i]] <- edgesDF
    
}

complexesParticipantsDF <- do.call("rbind", complexesParticipantsList)
edgesComplexes <- do.call("rbind", complexesEdgesList)


# Get the number of complexes per gene

nComplexesPerProteinDF <- complexesParticipantsDF %>%
    group_by(
        protein
    ) %>%
    summarise(
        n = n()
    ) %>%
    rename(
        accession = protein
    )

rnaComplexesDF <- accessionsMapping %>%
    filter(
        !gene %in% excludedGenes
    ) %>% 
    left_join(
        nComplexesPerProteinDF,
        by = "accession"
    ) %>%
    mutate(
        n = ifelse(is.na(n), 0, n)
    ) %>%
    select(
        -accession
    ) %>%
    distinct() %>%
    group_by(
        gene
    ) %>%
    summarise(
        n = sum(n)
    ) %>% 
    ungroup() %>%
    inner_join(
        rnaProteinCorLongDF,
        by = "gene"
    ) 


# Build plot

rnaComplexesDF$nFactor <- factor(
    x = ifelse(rnaComplexesDF$n > 8, ">8", 
               ifelse(rnaComplexesDF$n > 4, "5-8", as.character(rnaComplexesDF$n))),
    levels = as.character(c(0:9, "5-8", ">8"))
)

levels <- levels(rnaComplexesDF$nFactor)

for (i in 1:length(levels)) {
    
    level <- levels[i]
    
    n <- sum(rnaComplexesDF$nFactor == level)
    
    level <- paste0(level, "\n(n = ", n, ")")
    
    levels[i] <- level
    
}

levels(rnaComplexesDF$nFactor) <- levels

ggplot(
    data = rnaComplexesDF
) +
    geom_violin(
        mapping = aes(
            x = nFactor,
            y = z_correlation,
            col = data
        ), 
        alpha = 0.1,
        width = 0.5
    ) +
    geom_boxplot(
        mapping = aes(
            x = nFactor,
            y = z_correlation,
            col = data
        ),
        width = 0.5
    ) +
    scale_x_discrete(
        name = "Number of Complexes"
    ) +
    scale_y_continuous(
        name = "RNA-Protein correlation"
    ) +
    scale_color_manual(
        values = c("darkblue", "black")
    ) +
    theme(
        legend.position = "top",
        legend.title = element_blank()
    )
```

![](rna-protein_files/figure-gfm/complexes-1.png)<!-- -->

##### :speech\_balloon: According to this analysis, does the participation in complexes seem to be associated with CNA attenuation? How does this compare to CNA attenuation?

##### :pencil2: Load pathway data, plot the RNA-protein correlation against on the number of reactions a protein is involved in.

``` r
# Protein reactions from Reactome

reactomeSearchDF <- read.table(
    file = "resources/data/search.tsv.gz", 
    header = T, 
    sep = "\t", 
    stringsAsFactors = F, 
    quote = "", 
    comment.char = ""
)


# Get the number of reactions per gene

nReactionsPerProteinDF <- reactomeSearchDF %>%
    select(
        UNIPROT, REACTION_STID
    ) %>%
    rename(
        protein = UNIPROT,
        reaction = REACTION_STID
    ) %>%
    group_by(
        protein
    ) %>%
    summarise(
        n = n()
    ) %>%
    rename(
        accession = protein
    )

rnaReactionsDF <- accessionsMapping %>%
    filter(
        !gene %in% excludedGenes
    ) %>% 
    left_join(
        nReactionsPerProteinDF,
        by = "accession"
    ) %>%
    mutate(
        n = ifelse(is.na(n), 0, n)
    ) %>%
    select(
        -accession
    ) %>%
    distinct() %>%
    group_by(
        gene
    ) %>%
    summarise(
        n = sum(n)
    ) %>% 
    ungroup() %>%
    inner_join(
        rnaProteinCorLongDF,
        by = "gene"
    )


# Build plot

ggplot(
    data = rnaReactionsDF
) +
    geom_point(
        mapping = aes(
            x = n,
            y = z_correlation,
            col = data
        ), 
        alpha = 0.1
    ) +
    geom_smooth(
        mapping = aes(
            x = n,
            y = z_correlation,
            col = data
        ),
        method = "loess"
    ) +
    scale_x_log10(
        name = "Number of Reactions"
    ) +
    scale_y_continuous(
        name = "Correlation [Z-score]"
    ) +
    scale_color_manual(
        values = c("darkblue", "black")
    ) +
    theme(
        legend.position = "top",
        legend.title = element_blank()
    )
```

    ## Warning: Transformation introduced infinite values in continuous x-axis
    
    ## Warning: Transformation introduced infinite values in continuous x-axis

    ## Warning: Removed 7608 rows containing non-finite values (stat_smooth).

![](rna-protein_files/figure-gfm/pathways-1.png)<!-- -->

##### :speech\_balloon: According to this analysis, does the participation in biochemical reactions seem to be associated with RNA correlation? How does this compare to CNA attenuation?

##### :pencil2: Load PPI data, plot the RNA-protein correlation against on the number of interactors a protein has.

``` r
# Protein interactions from String

edgesString <- read.table(
    file = "resources/data/string_v10.5_high_edges.gz", 
    header = T, 
    sep = " ", 
    stringsAsFactors = F, 
    quote = "", 
    comment.char = ""
)


# Build and clean graph

graphString <- graph_from_data_frame(edgesString)
graphString <- simplify(graphString, remove.multiple = T, remove.loops = T, edge.attr.comb = "first")


# Build degree data frame

degreeDF <- data.frame(
    accession = V(graphString)$name,
    degree = degree(graphString),
    stringsAsFactors = F
)


# Get the number of interactors per gene

rnaDegreeDF <- accessionsMapping %>%
    filter(
        !gene %in% excludedGenes
    ) %>% 
    left_join(
        degreeDF,
        by = "accession"
    ) %>%
    mutate(
        degree = ifelse(is.na(degree), 0, degree)
    ) %>%
    select(
        -accession
    ) %>%
    distinct() %>%
    group_by(
        gene
    ) %>%
    summarise(
        degree = sum(degree)
    ) %>% 
    ungroup() %>%
    inner_join(
        rnaProteinCorLongDF,
        by = "gene"
    )


# Build plot

ggplot(
    data = rnaDegreeDF
) +
    geom_point(
        mapping = aes(
            x = log10(degree),
            y = z_correlation,
            col = data
        ), 
        alpha = 0.1
    ) +
    geom_smooth(
        mapping = aes(
            x = log10(degree),
            y = z_correlation,
            col = data
        ),
        method = "loess"
    ) +
    scale_x_continuous(
        name = "Number of Interactors"
    ) +
    scale_y_continuous(
        name = "Correlation [Z-score]"
    ) +
    scale_color_manual(
        values = c("darkblue", "black")
    ) +
    theme(
        legend.position = "top",
        legend.title = element_blank()
    )
```

    ## Warning: Removed 3190 rows containing non-finite values (stat_smooth).

![](rna-protein_files/figure-gfm/string-1.png)<!-- -->

##### :speech\_balloon: According to this analysis, does the number of interactors seem to be associated with RNA correlation? How does this compare to reactions? To CNA attenuation?

## References

1)  [Breast cancer quantitative proteome and proteogenomic
    landscape](https://www.ncbi.nlm.nih.gov/pubmed/30962452)
