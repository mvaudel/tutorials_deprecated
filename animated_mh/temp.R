

# Load packages
library(digest, lib.loc = lib)
library(reshape2, lib.loc = lib)
library(labeling, lib.loc = lib)
library(withr, lib.loc = lib)
library(ggplot2, lib.loc = "~/R")
library(scico, lib.loc = "~/R")
library(gganimate, lib.loc = "~/R")

# Set ggplot theme once and for all
theme_set(theme_bw(base_size = 11))




#' Returns a gganimate object with the MH for the given association data frame.
#' 
#' @param associationDF data frame containing the results of the association
#' @param pathwaysDF pathway matching data frame
#' @param orderedPathways an ordered vector of pathways
#' @param chrColumn name of the chromosome column
#' @param bpColumn name of the chromosomic coordinates column
#' @param idColumn name of the ids column
#' @param pColumn name of the log-transformed p-values column
#' @param pathwayColumn name of the pathway column
#' @param snpColors vector of two colors to use for markers alternatively on chromosomes
#' @param annotationColor color to use for annotated markers
#' @param maxP the maximal p-value used to scale the y axis
#' @param threshold the p-value threshold to display
#' @param thresholdColor the color to use for the threshold
#' 
#' @return a ggplot object of the MH
getAnimatedMh <- function(
  associationDF, 
  pathwaysDF,
  orderedPathways,
  chrColumn = "chr", 
  bpColumn = "bp", 
  idColumn = "id", 
  pColumn = "p", 
  pathwayColumn = "topPathway",
  snpColors = scico(n = 2, begin = 0.2, end = 0.4, palette = "grayC"), 
  annotationColor = scico(n = 1, begin = 0.8, end = 0.8, palette = "grayC"), 
  maxP = max(associationDF[[pColumn]]),
  threshold = -log10(5e-8),
  thresholdColor = "green4",
  pathwayPalette = "batlow") {
  
  # Chromosome coordinates 
  chromosomeStart <- cumsum(chromosomeLength) - chromosomeLength
  chromosomeMiddle <- chromosomeStart + chromosomeLength / 2
  
  # Make data frame for plotting
  baseDF <- data.frame(chr = associationDF[[chrColumn]],
                       bp = associationDF[[bpColumn]],
                       id = associationDF[[idColumn]],
                       p = associationDF[[pColumn]],
                       stringsAsFactors = F)
  
  # Get genomic coordinates
  baseDF$x <- chromosomeStart[baseDF$chr] + baseDF$bp
  
  # Set base colors
  baseDF$color <- ifelse(as.numeric(baseDF$chr) %% 2 == 0, snpColors[1], snpColors[2])
  
  # Get pathway colors
  pathwayColors <- scico(n = length(orderedPathways), palette = pathwayPalette)
  
  # Make one data frame per pathway with colored markers and merge them
  pathwayDFs <- list()
  for (i in 1:length(orderedPathways)) {
    
    pathway <- orderedPathways[i]
    pathwayColor <- pathwayColors[i]
    
    pathwayDF <- baseDF
    pathwayDF$pathway <- pathway
    
    matchedIds <- unique(pathwaysDF$id[pathwaysDF[[pathwayColumn]] == pathway])
    pathwayDF$color[pathwayDF$id %in% matchedIds] <- pathwayColor
    
    pathwayDFs[[i]] <- pathwayDF
    
  }
  plotDF <- do.call("rbind", pathwayDFs)
  
  # Factor the colors and pathways, and sort the data frame with the annotated SNPs in the front
  plotDF$pathway <- factor(plotDF$pathway, levels = orderedPathways)
  plotDF$color <- factor(plotDF$color, levels = c(snpColors, pathwayColors))
  plotDF <- plotDF[order(plotDF$pathway, plotDF$color, plotDF$x), ]
  
  # X-axis labels
  xLabels <- chromosomes
  xLabelsI <- 1:length(xLabels)
  xLabels[xLabelsI %% 2 == 0 & xLabelsI > 17] <- ""
  
  # Threshold DF
  thresholdDF <- data.frame(threshold = rep(threshold, length(orderedPathways)), pathway = orderedPathways, stringsAsFactors = F)
  thresholdDF$pathway <- factor(thresholdDF$pathway, levels = orderedPathways)
  
  # Make plot
  mhPlot <- ggplot() + 
    
    geom_hline(data = thresholdDF, mapping = aes(yintercept = threshold), col = thresholdColor, size = 0.3) + 
    geom_point(data = plotDF, mapping = aes(x = x, y = p, col = color), alpha = 0.8, size = 1, stroke = 0, shape = 16) + 
    
    scale_y_continuous(name = "p-value [-log10]", expand = c(0, 0), limits = c(0, 1.05 * maxP)) + 
    scale_x_continuous(name = "Chromosome", breaks = chromosomeMiddle, labels = xLabels, limits = c(0, genomeLength), expand = c(0.01, 0.01)) + 
    scale_color_manual(values = c(snpColors, pathwayColors)) + 
    
    ggtitle("{closest_state}") +
    
    theme(legend.position = "none",
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(size = 0.3),
          plot.title = element_text(hjust = 1)) +
    
    transition_states(states = pathway, transition_length = 3, state_length = 5) +
    enter_fade() +
    exit_fade()
  
  return(mhPlot)
  
}

print(paste0(Sys.time(), " Loading Giant"))

# Read table
giantData <- read.table(file = "resources/bmi_p-values.gz", header = T, stringsAsFactors = F, sep = "\t")

# Select markers and columns of interest
giantData <- giantData[giantData$SNPNAME != '-' & !is.na(giantData$Pvalue), c("CHR", "POS", "SNPNAME", "Pvalue")]

# Formatting
names(giantData) <- c("chr", "bp", "id", "p")
giantData$p <- -log10(giantData$p)

# Chromosome lengths in GRCh37.p13 (hg19) from Ensembl
chromosomes <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X")
chromosomeLength <- c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747,	135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566, 	155270560)
genomeLength <- sum(chromosomeLength)

# Factor and sort markers
giantData$chr <- factor(giantData$chr, levels = chromosomes)
giantData <- giantData[order(giantData$chr, giantData$bp), ]


print(paste0(Sys.time(), " Loading Pathways"))

# Read table
pathwaysDF <- read.table(file = "resources/search.tsv.gz", header = T, stringsAsFactors = F, sep = "\t", quote = "", comment.char = "")

# Select columns of interest
pathwaysDF <- pathwaysDF[, c("RSID", "PATHWAY_DISPLAY_NAME", "TOP_LEVEL_PATHWAY_DISPLAY_NAME")]
names(pathwaysDF) <- c("id", "pathway", "topPathway")

# Get pathways
pathways <- unique(pathwaysDF$pathway)
topLevelPathways <- unique(pathwaysDF$topPathway)

# Extract number of pathways
nPathways <- length(pathways)
nTopLevelPathways <- length(topLevelPathways)
print(paste("Number of top level pathways: ", nTopLevelPathways))
print(paste("Number of pathways: ", nPathways))


print(paste0(Sys.time(), " Clustering pathways"))

# Get the number of markers mapped to each pathway
pathwaySnps <- list()

for (i in 1:nTopLevelPathways) {
  
  pathwayI <- topLevelPathways[i]
  pathwaySnps[[i]] <- unique(pathwaysDF$id[pathwaysDF$topPathway == pathwayI])
  
}

# Make a matrix of number of shared markers
pathwayOverlap <- matrix(nrow = nTopLevelPathways, ncol = nTopLevelPathways, data = 0, dimnames = list(topLevelPathways, topLevelPathways))
diag(pathwayOverlap) <- 1

for (i in 1:(nTopLevelPathways-1)) {
  for (j in (i+1):nTopLevelPathways) {
    
    pathwayI <- topLevelPathways[i]
    pathwayJ <- topLevelPathways[j]
    
    pathwaySnpsI <- pathwaySnps[[i]]
    pathwaySnpsJ <- pathwaySnps[[j]]
    
    nI <- length(pathwaySnpsI)
    nJ <- length(pathwaySnpsJ)
    
    nIJ <- sum(pathwaySnpsI %in% pathwaySnpsJ)
    
    pathwayOverlap[i, j] <- nIJ / (nI + nJ - nIJ)
    pathwayOverlap[j, i] <- pathwayOverlap[i, j]
    
  }
}

# Make distance matrix
pathwayDist <- 1 - pathwayOverlap

# Cluster
pathwayClust <- hclust(as.dist(pathwayDist))
orderedPathways <- topLevelPathways[pathwayClust$order]


print(paste0(Sys.time(), " Animate"))

# Animate
animatedMH <- getAnimatedMh(associationDF = giantData, pathwaysDF = pathwaysDF, orderedPathways = orderedPathways)

animate(animatedMH, nframes = 500, height = 900, width = 1600)
anim_save(filename = "pathwayMh.gif", path = "plots")