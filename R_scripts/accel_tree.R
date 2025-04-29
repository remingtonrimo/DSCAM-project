colorChronogram <- function(phy, phy2, node.vals, col = c("green", "orange"), rotate = TRUE, assoc = NULL, 
                            edge.width=2, cex=1, ...){
  
  # Get tip labels associations if not provided
  if (is.null(assoc)) {
    assoc <- intersect(phy$tip.label, phy2$tip.label)
    assoc <- if (length(assoc) > 0) cbind(assoc, assoc) else NULL
    if (is.null(assoc)) stop("No shared tip labels to associate.")
  }
  
  # Optimize rotation if requested
  if (rotate) {
    cat("Rotating nodes to optimize matching...\n")
    assoc_map <- setNames(sapply(assoc[, 2], match, table = phy2$tip.label), assoc[, 1])
    phy <- tipRotate(phy, assoc_map * Ntip(phy) / Ntip(phy2), right = phy2, assoc = assoc, ...)
    
    assoc_map <- setNames(sapply(assoc[, 1], match, table = phy$tip.label), assoc[, 2])
    phy2 <- tipRotate(phy2, assoc_map * Ntip(phy2) / Ntip(phy), left = phy, assoc = assoc, ...)
    cat("Done with rotation optimization.\n")
  }
  
  #Convert node map to edge map by averaging rootward and tipward node values:
  edge.vector <- c()
  for(edge.index in 1:dim(phy$edge)[1]){
    edge.vector[edge.index] <- (node.vals[phy$edge[edge.index,1]] + node.vals[phy$edge[edge.index,2]])/2
  }
  
  #Get limits and get us a buffer:
  rate.lims <- range(edge.vector)
  lims.percentage.correction <- 0.001
  rate.lims[1] <- rate.lims[1] - lims.percentage.correction*abs(rate.lims[1])
  rate.lims[2] <- rate.lims[2] + lims.percentage.correction*abs(rate.lims[2])
  rate_normalized <- round((rate.lims[2] - edge.vector)/(rate.lims[2] - rate.lims[1]), 3)*1000
  rate.cols <- rev(colorRampPalette(col, space="Lab")(1001))
  phy$edge.colors <- rate.cols[rate_normalized]
  
  plot.phylo(x=phy2, edge.color = "black", edge.width=edge.width, cex = cex, show.tip.label = FALSE, ...)
  
  par(new = TRUE)
  plot.phylo(x=phy, edge.color = phy$edge.colors, edge.width=edge.width, cex = cex, show.tip.label = FALSE, ...)
  
}

# Plot t1 with color mapping and overlay t2 in black
colorChronogram(phy = exon4_time, phy2 = pruned_dscam1, node.vals = values_vector, col = c("#fde725","#21918c","#440154"), 
                rotate = TRUE, edge.width = 2, cex = 1)



node.height(exon4_time)
nodeheight(pruned_dscam1, 100)
a <-nodeheight(exon4_time, 100)
b <-nodeheight(pruned_dscam1, 100)
a-b
a
b

colorChronogram <- function(phy, phy2, col = c("green", "orange"), rotate = TRUE, assoc = NULL, 
                            edge.width = 2, cex = 1, ...) {
  
  # Get tip labels associations if not provided
  if (is.null(assoc)) {
    assoc <- intersect(phy$tip.label, phy2$tip.label)
    assoc <- if (length(assoc) > 0) cbind(assoc, assoc) else NULL
    if (is.null(assoc)) stop("No shared tip labels to associate.")
  }
  
  # Optimize rotation if requested
  if (rotate) {
    cat("Rotating nodes to optimize matching...\n")
    assoc_map <- setNames(sapply(assoc[, 2], match, table = phy2$tip.label), assoc[, 1])
    phy <- tipRotate(phy, assoc_map * Ntip(phy) / Ntip(phy2), right = phy2, assoc = assoc, ...)
    
    assoc_map <- setNames(sapply(assoc[, 1], match, table = phy$tip.label), assoc[, 2])
    phy2 <- tipRotate(phy2, assoc_map * Ntip(phy2) / Ntip(phy), left = phy, assoc = assoc, ...)
    cat("Done with rotation optimization.\n")
  }
  
  # creating a vector a branch lengths differences
  a <-length(phy$tip.label) + phy$Nnode
  node.vals <- numeric(a) 
  
  # Loop through each node using `i`
  for (i in 1:374) {
    b <- nodeheight(exon4_time, i)
    c <- nodeheight(pruned_dscam1, i)
    node.vals[i] <- b - c
  }
  
  # Convert node differences to edge values
  edge.vector <- numeric(nrow(phy$edge))
  for (edge.index in 1:nrow(phy$edge)) {
    edge.vector[edge.index] <- (node.vals[phy$edge[edge.index, 1]] + node.vals[phy$edge[edge.index, 2]]) / 2
  }
  
  # Set explicit symmetric range
  adjusted.lims <- c(-50, 50)
  # Clamp values to the adjusted range
  edge.vector.clamped <- pmax(pmin(edge.vector, adjusted.lims[2]), adjusted.lims[1])
  # Normalize values to the range [0, 1] for coloring
  rate_normalized <- (edge.vector.clamped - adjusted.lims[1]) / (adjusted.lims[2] - adjusted.lims[1])
  # Generate colors
  rate.cols <- rev(colorRampPalette(col, space="Lab")(1001))  # Create a color gradient
  color_indices <- round(rate_normalized * 1000) + 1  # Map normalized values to color indices
  phy$edge.colors <- rate.cols[color_indices]  # Assign colors
  
  # Plot trees
  plot.phylo(x = phy2, edge.color = "black", edge.width = edge.width, cex = cex, show.tip.label = FALSE, ...)
  par(new = TRUE)
  plot.phylo(x = phy, edge.color = phy$edge.colors, edge.width = edge.width, cex = cex, show.tip.label = FALSE, ...)
}



colorChronogram(phy = exon4_time, phy2 = pruned_dscam1, col = c("#fde725","#21918c","#440154"), 
                rotate = TRUE, edge.width = 2, cex = 1)



length(exon4_time$Nnode)
# creating a vector a branch lengths differences
a <-length(exon4_time$tip.label) + exon4_time$Nnode
node.vals <- numeric(a) 






getDescendants(exon4_time, 220)


