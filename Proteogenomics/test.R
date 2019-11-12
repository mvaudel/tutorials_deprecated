
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


rnaProteinCorDF$x <- 0

model <- gamlss(
    formula = as.formula("mRNA_protein_correlation ~ x"),
    family = NO,
    data = rnaProteinCorDF[, c("mRNA_protein_correlation", "x")]
)

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

# Exclude gene mapping with more than 100 proteins per gene

geneOccurenceDF <- as.data.frame(
    x = table(accessionsMapping$gene), 
    stringsAsFactors = F
)
names(geneOccurenceDF) <- c("gene", "nProteins")

excludedGenes <- geneOccurenceDF$gene[geneOccurenceDF$nProteins >= 100]


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



reactomeSearchDF <- read.table(
    file = "resources/data/search.tsv.gz", 
    header = T, 
    sep = "\t", 
    stringsAsFactors = F, 
    quote = "", 
    comment.char = ""
)

pathwaysDF <- reactomeSearchDF %>%
    select(
        PATHWAY_STID, PATHWAY_DISPLAY_NAME
    ) %>%
    distinct() %>%
    mutate(
        z_mRNA_protein_correlation = NA,
        z_random_correlation = NA,
        p = NA
    )

for (i in 1:nrow(pathwaysDF)) {
    
    pathwayId <- pathwaysDF$PATHWAY_STID[i]
    proteins <- reactomeSearchDF$UNIPROT[pathwaysDF$PATHWAY_STID == pathwayId]
    genes <- accessionsMapping$gene[accessionsMapping$accession %in% proteins]
    
    z_mRNA_protein_correlation <- rnaProteinCorDF$z_mRNA_protein_correlation[rnaProteinCorDF$gene %in% genes]
    z_random_correlation <- rnaProteinCorDF$z_random_correlation[rnaProteinCorDF$gene %in% genes]
    
    tTest <- t.test(
        x = z_random_correlation,
        y = z_mRNA_protein_correlation,
        alternative = "two.sided"
    )
    
    p <- tTest$p.value
    
    pathwaysDF$z_mRNA_protein_correlation[i] <- mean(z_mRNA_protein_correlation)
    pathwaysDF$z_random_correlation[i] <- mean(z_random_correlation)
    pathwaysDF$p[i] <- p
    
}

pathwaysDF <- pathwaysDF %>%
    arrange(
        desc(p)
    ) %>%
    mutate(
        bhFDR <- p.adjust(p, method = "BH")
    )

ggplot(
    data = pathwaysDF
) + 
    geom_point(
        mapping = aes(
            x = z_mRNA_protein_correlation,
            y = -log10(p),
            col = bhFDR
        ),
        alpha = 0.5
    ) +
    scale_x_continuous(
        name = "Correlation [Z-score]"
    ) +
    scale_y_continuous(
        name = "p-value [-log10]"
    ) + 
    scale_color_scico(
        name = "BH FDR",
        palette = "turku",
        direction = -1,
        end = 0.8
    ) +
    theme(
        legend.position = "top"
    )

View(pathwaysDF)


expectedPvalues <- -log10(sort(runif(n = nrow(pathwaysDF))))
observedPvalues <- -log10(sort(pathwaysDF$p))

maxP <- max(c(expectedPvalues, observedPvalues))

ggplot() +
    geom_point(
        mapping = aes(
            x = expectedPvalues,
            y = observedPvalues
        )
    )

