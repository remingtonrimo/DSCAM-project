library(ape)
library(gplots)
#setting the working directory for the correct time treees
setwd("C:/Users/remin/Desktop/School/Programs/BEAST.v.2.7.6/DSCAM/chronogram")

####reading in dscam & mitochondria trees
  dscam<-read.tree("dscamdrop02.new")
  plot(dscam,cex=0.6)
  mito<-read.tree("mitodrop02.new")
  plot(mito,cex=0.6)
  nodelabels()

####gathering branch lengths for DSCAM and mitochondria
  dscamdist<-cophenetic(dscam)
  mitodist<-cophenetic(mito)
#dividing distances by 2
  dscam_dist2<-dscamdist/2
  mito_dist2<-mitodist/2
####writing out as cv
  write.csv(dscam_dist2, "dscam_bl.csv")
  write.csv(mito_dist2, "mito_bl.csv")
####reading in csv file
  diff<-read.csv("dscam-mito.csv")
  diff[upper.tri(diff)] <- NA
####makinng a heatmap with mito-dscam
  m<-as.matrix(diff[,-1])
  rownames <- diff$Species
  
# Assign row names to the matrix
  rownames(m) <- rownames
  
  heatmap((m),
    col = colorRampPalette(c("#05041A","#440154","#21918c", "limegreen","#fde725"))(100),
    symm = TRUE,
    main = "Cophenetic Matrix Heatmap",
    cexRow = 0.6,
    cexCol=0.6,
    Colv = NA,
    Rowv = NA)

####Heat map with a continuous legend  
  #The heat map I used
#Color palettes were from https://waldyrious.net/viridis-palette-generator/
  viridis <- colorRampPalette(c("#440154","#383167","#3b528b","#21918c","#A8CF2D"))(100)
  inferno <- colorRampPalette(c("#000004","#280b53","#57106e","#9f2a63","#bc3754","#f98e09"))(100)
  magma <- colorRampPalette(c("#000004","#51127c","#b73779","#fc8961","#fcfdbf"))(100)
  plasma <- colorRampPalette(c("#0d0887", "#7e03a8","#cc4778", "#f89540","#f0f921"))(100)
  bluered <- colorRampPalette(c("#0000CC","#400099","#800066","#bf0033","#FF0000"))(100)
  bloodymary <- colorRampPalette(c("#000004","#440154","#68157B","#DD2476","#EE3B53","#FF512F"))(100)
  custom<-colorRampPalette(c("#69000C","#7D0025","#DA3500","#F39300"))(100)
  
  heatmap.2(m,
            col = custom,
            symm = TRUE,
            main = "Cophenetic Matrix Heatmap", #Title of graph
            cexRow = 0.6,
            cexCol = 0.6,
            Colv = NA,
            Rowv = NA,
            key = TRUE,  # Display the legend
            key.title = "T(d)-T(m)",  # Title of the legend
            key.xlab = "Branch length difference",  # X-label
            key.ylab = "Distribution",# Y-label
            trace = "none",
            mar = c(9, 9)
  )

