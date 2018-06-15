# 
# This script builds MH and QQ plots.
#
startTimeAll <- proc.time()


## Libraries

library(ggplot2)


## Functions

mh <- function(x, snp='SNPNAME', chr='CHR', bp='POS', p='Pvalue', thresholdLow = 5, thresholdHigh = -log10(5e-8), thresholdLowColor = "darkgreen", thresholdHighColor = "green3"){
  
  # constants
  
  chromosomeLength <- c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747,	135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566, 155270560) # bp length from Ensembl GRCh37.p13
  genomeLength <- sum(chromosomeLength)
  
  # Make a data frame for the plot
  
  manhattanData <- data.frame(snp = x[[snp]], chr = x[[chr]], bp = x[[bp]], p = x[[p]], stringsAsFactors = F)
  manhattanData$chr[manhattanData$chr == 'X'] <- 23
  manhattanData$chr <- as.numeric(manhattanData$chr)
  manhattanData <- manhattanData[!is.na(manhattanData$p) & manhattanData$p > 0, ]
  manhattanData$logP <- -log10(manhattanData$p)
  
  
  # Map chromosomic coordinates to x axis
  
  xValues <- c() # The position of every SNP on the x axis
  xBreak <- c() # The center of every chromosome
  xBreakLabels <- c() # The labels to use for every chromosome
  chrStart <- c() # The start of every second chromosome
  chrEnd <- c() # The end of every second chromosome
  
  manhattanData <- manhattanData[order(manhattanData$chr, manhattanData$bp), ]
  chromosomes <- sort(unique(manhattanData$chr))
  
  xOffset <- 0
  start <- T
  xMin <- -1
  xMax <- -1
  
  for (chromosomeNumber in chromosomes) {
    
    tempLength <- chromosomeLength[chromosomeNumber]
    
    bpTemp <- manhattanData$bp[manhattanData$chr == chromosomeNumber]
    
    if (length(bpTemp) > 0) {
      
      if (xMin == -1) {
        
        xMin <- xOffset
        
      }
      
      xMax <- xOffset + tempLength
      
      xTemp <- bpTemp + xOffset
      xValues <- c(xValues, xTemp)
      
    }
    
    if (length(bpTemp) > 0 || xMin != -1) {
      
      breakValue <- xOffset + tempLength / 2
      xBreak <- c(xBreak, breakValue)
      
      if (chromosomeNumber < 12 || chromosomeNumber %% 2 == 0) {
        
        xBreakLabels <- c(xBreakLabels, chromosomeNumber)
        
      } else if (chromosomeNumber == 23) {
        
        xBreakLabels <- c(xBreakLabels, "X")
        
      } else {
        
        xBreakLabels <- c(xBreakLabels, "")
        
      }
    }
    
    xOffset <- xOffset + tempLength
    
    if (start) {
      
      chrStart <- c(chrStart, xOffset)
      
    } else {
      
      chrEnd <- c(chrEnd, xOffset)
      
    }
    
    start <- !start
    
  }
  
  manhattanData$xValues <- xValues
  
  # y axis scale
  
  yMax <- max(round(max(manhattanData$logP)+1), thresholdHigh)
  
  
  # Color according to chromosome
  
  manhattanData$color <- manhattanData$chr %% 2
  
  manhattanData <- manhattanData[order(manhattanData$color), ]
  
  manhattanData$color <- factor(manhattanData$color)
  
  
  # Make a ggplot object
  
  manhattanPlot <- ggplot()
  
  
  # Add background rectangle for every second chromosome
  
  # manhattanPlot <- manhattanPlot + geom_rect(aes(xmin = chrStart, xmax = chrEnd, ymin = 0, ymax = yMax), alpha = 0.2)
  
  
  # Plot all markers
  
  manhattanPlot <- manhattanPlot + geom_point(data = manhattanData, aes(x = xValues, y = logP, col = color), size = 1.5)
  
  
  # Plot thresholds
  
  manhattanPlot <- manhattanPlot + geom_hline(aes(yintercept = thresholdLow), col = thresholdLowColor)
  manhattanPlot <- manhattanPlot + geom_hline(aes(yintercept = thresholdHigh), col = thresholdHighColor)
  
  
  # Set axes labels and limits
  
  manhattanPlot <- manhattanPlot + scale_y_continuous(name = "-log10(p)", limits = c(0, yMax), expand = c(0, 0))
  manhattanPlot <- manhattanPlot + scale_x_continuous(name = NULL, breaks = xBreak, labels = xBreakLabels, expand = c(0.01, 0), limits = c(0, genomeLength))
  
  
  # Set colors
  
  manhattanPlot <- manhattanPlot + scale_color_manual(name = NULL, values = c("grey20", "grey40"))
  
  
  # Format background and grid
  
  manhattanPlot <- manhattanPlot + theme_bw(base_size = 14)
  manhattanPlot <- manhattanPlot + theme(panel.background = element_rect(fill = NA, colour = "grey50"),
                                         panel.grid.major.y = element_line(colour = "grey50", linetype = "dotted"),
                                         panel.grid.minor.y = element_blank(),
                                         panel.grid.major.x = element_blank(),
                                         panel.grid.minor.x = element_blank(),
                                         legend.position = "none")
  
  
  # Return the plot
  
  return(manhattanPlot)
  
}


qq <- function(x, snp='SNPNAME', p='Pvalue'){
  
  # Prepare data in a data frame
  
  qqData <- data.frame(snp = x[[snp]], p = x[[p]], stringsAsFactors = F)
  qqData <- qqData[!is.na(qqData$p) & qqData$p > 0, ]
  qqData$logP <- -log10(qqData$p)
  qqData$expectedLogP <- NA
  
  
  # Estimate expected p-value
  
  qqData <- qqData[order(qqData$p), ]
  qqData$expectedLogP <- -log10(ppoints(n = nrow(qqData)))
  
  
  # axis scale and labels
  
  yMax <- max(round(max(qqData$logP)+1))
  yMaxFloored <- floor(yMax)
  
  xMax <- floor(max(qqData$expectedLogP)) + 1
  
  
  # Make a ggplot object
  
  qqPlot <- ggplot()
  
  
  # Plot all markers
  
  qqPlot <- qqPlot + geom_point(data = qqData, aes(x = expectedLogP, y = logP), size = 1.5)
  
  
  # Plot thresholds
  
  qqPlot <- qqPlot + geom_abline(intercept = 0, slope = 1, linetype="dotted", color = "darkgray")
  
  
  # Plot pathway legend
  
  # qqPlot <- qqPlot + geom_text_repel(data = qqData, aes(x = expectedLogP, y = logP, label = annotation), nudge_x = 10, nudge_y = 0, segment.color = NA, na.rm = T)
  
  
  # Set axes labels and limits
  
  qqPlot <- qqPlot + scale_y_continuous(name = "-log10(p) Observed", limits = c(0, yMax), expand = c(0, 0))
  qqPlot <- qqPlot + scale_x_continuous(name = "-log10(p) Expected", limits = c(0, xMax), expand = c(0, 0))
  
  
  # Format plot
  
  qqPlot <- qqPlot + theme_bw(base_size = 14)
  qqPlot <- qqPlot + theme(legend.position = "none")
  
  return(qqPlot)
  
}


## Script

# Load data

giantDF <- read.table("resources/bmi_p-values.gz", header = T, stringsAsFactors = F)

mhPlot <- mh(giantDF)

png("mh_giant.png", height = 600, width = 800)
plot(mhPlot)
dummy <- dev.off()


qqPlot <- qq(giantDF)

png("qq_giant.png", height = 600, width = 800)
plot(qqPlot)
dummy <- dev.off()

