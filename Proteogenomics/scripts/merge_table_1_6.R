##
#
# This script merges tables 1 and 6 using gene names.
# 
##

# Libraries

library(tidyr)
library(dplyr)

# Load data

tumorDF <- read.table(
    file = "pages/proteogenomics/resources/data/tumor.gz",
    header = T,
    sep = "\t",
    comment.char = "",
    quote = "",
    stringsAsFactors = F
)

genesDF <- read.table(
    file = "pages/proteogenomics/resources/data/genes.gz",
    header = T,
    sep = "\t",
    comment.char = "",
    quote = "",
    stringsAsFactors = F
) %>%
    select(
        -contains("nspsm")
    )

saavDF <- read.table(
    file = "pages/proteogenomics/resources/data/saav.gz",
    header = T,
    sep = "\t",
    comment.char = "",
    quote = "",
    stringsAsFactors = F
)


# Normalize the intensities of the SAAV peptides using the reference channel (131).

for (tmtColumn in c("tmt126", "tmt127N", "tmt127C", "tmt128N", "tmt128C", "tmt129N", "tmt129C", "tmt130N", "tmt130C", "tmt131")) {
    
    saavDF <- saavDF %>% 
        mutate(
            !!tmtColumn := !!sym(tmtColumn) / tmt131
        )
}

saavDF <- saavDF %>%
    select(
        -tmt131
    )


# Organize the SAAV peptides table using kits as column, keep one quantification value per peptide.

saavPeptidesProteinDF <- saavDF %>% 
    select(
        peptide_sequence,
        protein
    ) %>%
    distinct()

for (tmtSet in 1:5) {
    
    for (tmtColumn in c("tmt126", "tmt127N", "tmt127C", "tmt128N", "tmt128C", "tmt129N", "tmt129C", "tmt130N", "tmt130C")) {
        
        newColumn <- paste0(tmtColumn, "_set", tmtSet)
        saavPeptidesProteinDF[[newColumn]] <- NA
        
    }
    
    for (i in 1:nrow(saavPeptidesProteinDF)) {
        
        js <- saavDF$tmt_set == tmtSet & saavDF$peptide_sequence == saavPeptidesProteinDF$peptide_sequence[i]
        
        if (sum(js) > 0) {
            
            for (tmtColumn in c("tmt126", "tmt127N", "tmt127C", "tmt128N", "tmt128C", "tmt129N", "tmt129C", "tmt130N", "tmt130C")) {
                
                newColumn <- paste0(tmtColumn, "_set", tmtSet)
                
                values <- saavDF[js, tmtColumn]
                
                if (sum(!is.na(values)) > 0) {
                    
                    value <- median(
                        x = values, 
                        na.rm = T
                    )
                    
                    saavPeptidesProteinDF[i, newColumn] <- value
                    
                }
            }
        }
    }
}


# Rename columns using the tumor ids

tumorDF <- tumorDF %>%
    mutate(
        tmtColumn = paste0("tmt", tmtTag, "_set", tmtSet)
    )

for (i in 1:nrow(tumorDF)) {
    
    oldColumn <- tumorDF$tmtColumn[i]
    newColumn <- tumorDF$tumor_id[i]
    
    names(saavPeptidesProteinDF)[names(saavPeptidesProteinDF) == oldColumn] <- newColumn
    names(genesDF)[names(genesDF) == oldColumn] <- newColumn
    
}


# Extract gene names and merge with the protein table.

table16DF <- NULL

genesSplit <- strsplit(
    x = saavPeptidesProteinDF$protein,
    split = ";"
)

for (i in 1:nrow(saavPeptidesProteinDF)) {
    
    geneLines <- genesSplit[[i]]
    
    geneLines <- geneLines[
        startsWith(
            x = geneLines,
            "CanProVar"
        ) |
            startsWith(
                x = geneLines,
                "COSMIC"
            ) |
            startsWith(
                x = geneLines,
                "RNAediting"
            )
        ]
    
    peptideDF <- data.frame(
        gene = character(length(geneLines)),
        peptide_sequence = character(length(geneLines)),
        gene_variant = character(length(geneLines)),
        protein_variant = character(length(geneLines)),
        stringsAsFactors = F
    )
    
    processedGenes <- c()
    processedPeptides <- c()
    
    for (j in 1:length(geneLines)) {
        
        peptide_sequence <- saavPeptidesProteinDF$peptide_sequence[i]
        
        geneLine <- geneLines[j]
        
        gene <- NA
        geneVariant <- NA
        proteinVariant <- NA
        
        if (
            startsWith(
                x = geneLine,
                "CanProVar"
            )) {
            
            tempSplit <- strsplit(
                x = geneLine,
                split = "_"
            )[[1]]
            
            gene <- tempSplit[3]
            geneVariant <- tempSplit[2]
            proteinVariant <- tempSplit[4]
            
        } else if (
            startsWith(
                x = geneLine,
                "COSMIC"
            )) {
            
            tempSplit <- strsplit(
                x = geneLine,
                split = ":"
            )[[1]]
            
            gene <- tempSplit[2]
            geneVariant <- strsplit(
                x = tempSplit[4],
                split = "\\."
            )[[1]][2]
            proteinVariant <- strsplit(
                x = tempSplit[5],
                split = "\\."
            )[[1]][2]
            
            tempSplit <- strsplit(
                x = gene,
                split = "_"
            )[[1]]
            
            if (length(tempSplit) > 1) {
                
                gene <- tempSplit[1]
                
            }
        } else if (
            startsWith(
                x = geneLine,
                "RNAediting"
            )) {
            
            tempSplit <- strsplit(
                x = geneLine,
                split = "_"
            )[[1]]
            
            tempSplit2 <- strsplit(
                x = tempSplit[7],
                split = "\\."
            )[[1]]
            
            gene <- tempSplit[3]
            geneVariant <- paste(tempSplit[8], tempSplit[9], tempSplit[2], sep = ".")
            proteinVariant <- paste(tempSplit2[1], tempSplit[6], tempSplit2[2], sep = ".")
            
        }
        
        if (!gene %in% processedGenes | !peptide_sequence %in% processedPeptides) {
        
        peptideDF$peptide_sequence[j] <- peptide_sequence
        peptideDF$gene[j] <- gene
        peptideDF$gene_variant[j] <- geneVariant
        peptideDF$protein_variant[j] <- proteinVariant
        
        for (tumor in tumorDF$tumor_id) {
            
            peptideColumn <- paste0(tumor, "_saavPeptide")
            peptideValue <- saavPeptidesProteinDF[i, tumor]
            peptideDF[j, peptideColumn] <- peptideValue
            
            proteinColumn <- paste0(tumor, "_gene")
            
            k <- genesDF$gene_symbol == peptideDF$gene[j]
            
            if (sum(k) == 1) {
                
                proteinValue <- genesDF[k, tumor]
                
            } else {
                
                proteinValue <- NA
                
            }
            
            peptideDF[j, proteinColumn] <- proteinValue
            
        }
        
        processedGenes[length(processedGenes) + 1] <- gene
        processedPeptides[length(processedPeptides) + 1] <- peptide_sequence
        
        }
    }
    
    peptideDF <- peptideDF %>%
        filter(
            gene != ""
        ) %>%
        distinct()
    
    if (is.null(table16DF)) {
        
        table16DF <- peptideDF
        
    } else {
        
        table16DF <- rbind(table16DF, peptideDF)
        
    }
}

# Save to file

write.table(
    x = table16DF,
    file = gzfile("pages/proteogenomics/resources/data/table16.gz"),
    sep = "\t",
    col.names = T,
    row.names = F,
    quote = F
)
