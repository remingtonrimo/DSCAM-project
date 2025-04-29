# First we need to drop tips from the original phylogeny that aren't included in the exon alignment files 
# This tree will be included as the prior and constraint when constructing the phylogenies of the exons
plot(dscam1)
alignment <- read.nexus.data("variable_exon9.nex")

# Extract species names
tree_species <- tree$tip.label
alignment_species <- names(alignment)

# Check for species in the tree but not in the alignment
missing_in_alignment <- setdiff(tree_species, alignment_species)
if (length(missing_in_alignment) > 0) {
  cat("Species in the tree but not in the alignment:\n")
  print(missing_in_alignment)
} else {
  cat("All species in the tree are present in the alignment.\n")
}

# Check for species in the alignment but not in the tree
missing_in_tree <- setdiff(alignment_species, tree_species)
if (length(missing_in_tree) > 0) {
  cat("Species in the alignment but not in the tree:\n")
  print(missing_in_tree)
} else {
  cat("All species in the alignment are present in the tree.\n")
}


# list of species that disrupt the tree topology
notpresent <- c("Anopheles_coustani","Anopheles_bellator","Diorhabda_carinulata")

# We then drop the disruptive species
exon9guide <- drop.tip(dscam1, notpresent)

write.tree(exon9guide, file ="exon9guide.new")

library(ape)
library(phytools)
library(phangorn)
# change the directory to where all of the trees are
setwd("~/dscam")
  csv <- read.csv("dscam.csv")
  dscam1 <- read.newick("iqtree_run/dscam1.new")
setwd("~/dscam/variable_exons")

# These trees need to be dated, using dating_phylogenies.R script
exon4_tree <- read.nexus("exon4/exon4-con.nex")
  plot(exon4_tree, show.tip.label = TRUE, cex=0.4, edge.width = 2)
exon6_tree <- read.nexus("exon6/exon6-con.nex")
  plot(exon6_tree)
exon9_tree <- read.nexus("exon9/exon9-con.nex")
  plot(exon9_tree)

# adding underscores
  csv$Species <- gsub(" ", "_", csv$Species)
# Identify the first occurrence of each family
  first_occurrence <- csv[!duplicated(csv$Family), ]
# Get a list of all species
  all_species <- csv$Species
# Get a list of species to keep (the first occurrence for each family)
  species_to_keep <- first_occurrence$Species
# Identify species to remove
  species_to_remove <- setdiff(all_species, species_to_keep)
# Drop the species from the tree
  exon4_fam <- drop.tip(exon4_time, species_to_remove)
  exon6_fam <- drop.tip(exon6_time, species_to_remove)
  exon9_fam <- drop.tip(exon9_time, species_to_remove)
  dscam1_fam <- drop.tip(dscam1, species_to_remove)
# ceating a vector to rename the tips of the phylogeny
  rename_map <- setNames(first_occurrence$Family, first_occurrence$Species)
# Rename the tips in the tree
  exon4_fam$tip.label <- rename_map[exon4_fam$tip.label]
  exon6_fam$tip.label <- rename_map[exon6_fam$tip.label]
  exon9_fam$tip.label <- rename_map[exon9_fam$tip.label]
  dscam1_fam$tip.label <- rename_map[dscam1_fam$tip.label]
  
  tip_name <- exon9_fam$tip.label[8]  # Get the name of the 8th tip
  exon9_fam <- drop.tip(exon9_fam, tip_name)
  
# from here I need to take the difference between dscam1 branch lengths with the exon trees
# Then once I do that, I can incorporate node coloring into my function
dscamdist<-cophenetic(dscam1_fam)
exon4dist<-cophenetic(exon4_fam)
exon6dist<-cophenetic(exon6_fam)
exon9dist<-cophenetic(exon9_fam)

#dividing distances by 2
dscamdist<-dscamdist/2
exon4dist<-exon4dist/2
exon6dist<-exon6dist/2
exon9dist<-exon9dist/2

# making the matrix order
dscam1_fam$tip.label
matrix_order <- c("Coenagrionidae", "Acrididae","Archotermopsidae", "Bacillidae", "Thripidae", "Cicadellidae",
                  "Aleyrodidae", "Aphididae", "Adelgidae", "Phylloxeridae", "Pseudococcidae", "Braconidae",
                  "Bethylidae", "Vespidae", "Cephidae","Formicidae", "Megachilidae", "Colletidae","Halictidae",
                  "Apidae","Diprionidae","Agaonidae","Eulophidae","Figitidae","Chrysopidae","Scarabaeidae","Lampyridae",
                  "Silphidae","Staphylinidae","Coccinellidae","Chrysomelidae","Curculionidae","Brentidae",
                  "Tenebrionidae","Bombycidae","Papilionidae","Tortricidae","Cosmopterigidae","Gelechiidae",
                  "Noctuidae","Pieridae","Nymphalidae","Pyralidae","Crambidae","Pulicidae","Psychodidae",
                  "Sciaridae","Cecidomyiidae","Stratiomyidae","Dolichopodidae","Syrphidae","Glossinidae",
                  "Diopsidae","Tephritidae","Drosophilidae","Ceratopogonidae","Culicidae")
# Create an empty matrix for the results
diff4 <- matrix(0, nrow = length(matrix_order), ncol = length(matrix_order), dimnames = list(matrix_order, matrix_order))
diff6 <- matrix(0, nrow = length(matrix_order), ncol = length(matrix_order), dimnames = list(matrix_order, matrix_order))
diff9 <- matrix(0, nrow = length(matrix_order), ncol = length(matrix_order), dimnames = list(matrix_order, matrix_order))

# Loop through the rows and columns to perform the subtraction
# Loop through the rows and columns to perform the subtraction for all exon matrices
for (row_name in rownames(dscamdist)) {
  for (col_name in colnames(dscamdist)) {
    # Check for exon4dist
    if (row_name %in% rownames(exon4dist) && col_name %in% colnames(exon4dist)) {
      diff4[row_name, col_name] <- dscamdist[row_name, col_name] - exon4dist[row_name, col_name]
    }
    # Check for exon6dist
    if (row_name %in% rownames(exon6dist) && col_name %in% colnames(exon6dist)) {
      diff6[row_name, col_name] <- dscamdist[row_name, col_name] - exon6dist[row_name, col_name]
    }
    # Check for exon9dist
    if (row_name %in% rownames(exon9dist) && col_name %in% colnames(exon9dist)) {
      diff9[row_name, col_name] <- dscamdist[row_name, col_name] - exon9dist[row_name, col_name]
    }
  }
}

phylo <- plot.phylo(dscam1_fam, type= "radial", cex = 0.8, edge.width = 2)

diff4[upper.tri(diff4)] <- NA
diff6[upper.tri(diff6)] <- NA
diff9[upper.tri(diff9)] <- NA

col <- colorRampPalette(c("#fbb61a","#ed6925","#bc3754","#781c6d","#320a5e","#000004"))(100)
image(1:100, 1, as.matrix(1:100), col = col, axes = FALSE, xlab = "Value", ylab = "", 
      main = "Branch Differences")

breaks <- seq(-100, 200, length.out = 101)  # Scale from -100 to 100
heatmap.2(diff4,
          col = col,
          symm = TRUE,
          main = "Exon 4", #Title of graph
          cexRow = 0.6,
          cexCol = 0.6,
          dendrogram = "none",
          Colv = NA,
          Rowv = NA,
          key = T,# Display the legend
          key.title = "Δ=T(gene)-T(exon)",  # Title of the legend
          key.xlab = "Branch Differences",  # X-label
          key.ylab = "Distribution",# Y-label
          trace = "none",
          mar = c(9, 9),
          breaks = breaks
)
heatmap.2(diff6,
          col = col,
          symm = TRUE,
          main = "Exon 6", #Title of graph
          cexRow = 0.6,
          cexCol = 0.6,
          dendrogram = "none",
          Colv = NA,
          Rowv = NA,
          key = T,  # Display the legend
          key.title = "Δ=T(gene)-T(exon)",  # Title of the legend
          key.xlab = "Branch Differences",  # X-label
          key.ylab = "Distribution",# Y-label
          trace = "none",
          mar = c(9, 9),
          breaks = breaks
)
heatmap.2(diff9,
          col = col,
          symm = TRUE,
          main = "Exon 9", #Title of graph
          cexRow = 0.6,
          cexCol = 0.6,
          dendrogram = "none",
          Colv = NA,
          Rowv = NA,
          key = T,  # Display the legend
          key.title = "Δ=T(gene)-T(exon)",  # Title of the legend
          key.xlab = "Branch Differences",  # X-label
          key.ylab = "Distribution",# Y-label
          trace = "none",
          mar = c(9, 9),
          breaks = breaks
)

compare_trees <- list(dscam1_fam, exon4_fam,exon6_fam,exon9_fam)
kronoviz(compare_trees, horiz=FALSE, show.tip.label = TRUE)

# Overall, exon 4,6 and 9 all seem to have longer branch lengths in hymenoptera
# for the family level phylogenies

# extra code I was working with
# Function to rename nodes using unique descendant tips
rename_nodes_with_edge <- function(tree) {
  # Initialize a vector for new node labels
  new_node_labels <- rep(NA, tree$Nnode)
  
  # Iterate over internal nodes
  for (node_index in 1:tree$Nnode) {
    # Calculate the real node ID (internal nodes start after tips)
    real_node <- Ntip(tree) + node_index
    
    # Get all descendants of the current node
    descendants <- Descendants(tree, real_node, type = "tips")
    
    # Get tip labels for descendants
    if (length(descendants[[1]]) >= 2) {
      descendant_labels <- tree$tip.label[descendants[[1]]]
      # Use the first two unique tips to name the node
      new_node_labels[node_index] <- paste(descendant_labels[1], descendant_labels[2], sep = "_")
    } else if (length(descendants[[1]]) == 1) {
      # If only one descendant, name using it
      new_node_labels[node_index] <- tree$tip.label[descendants[[1]][1]]
    }
  }
  
  # Assign the new labels to the tree
  tree$node.label <- new_node_labels
  return(tree)
}

# Apply the function to rename nodes
exon4_time <- rename_nodes_with_edge(exon4_time)

# Check renamed nodes
exon4_time$node.label


# Create an empty vector to store the values
values_vector <- numeric(length(exon4_time$tip.label) + Nnode(exon4_time))
names(values_vector) <- c(exon4_time$tip.label, as.character(1:Nnode(exon4_time) + Ntip(exon4_time)))

# Assign 0 to all tip labels in values_vector
values_vector[1:Ntip(exon4_time)] <- 0

# Create a set (character vector) to track already assigned node labels
assigned_labels <- character()

# Loop over internal nodes to assign values from the difference matrix
for (node in (Ntip(exon4_time) + 1):(Ntip(exon4_time) + Nnode(exon4_time))) {
  # Get descendant tip indices for this node using Descendants
  descendants <- Descendants(exon4_time, node)[[1]]  # Use Descendants without "type" argument
  
  # Retrieve the species names corresponding to these descendant tips
  descendant_labels <- exon4_time$tip.label[descendants]
  
  # Create a unique node label for this internal node (e.g., "species1_species2")
  if (length(descendant_labels) == 2) {
    # Sort the species names to avoid label mismatches (e.g., Bombus_a_Bombus_b == Bombus_b_Bombus_a)
    node_label <- paste(sort(descendant_labels), collapse = "_")
    
    # If this label has not been assigned yet, assign the value from the difference matrix
    if (!(node_label %in% assigned_labels)) {
      # Check if the node label exists in the difference matrix
      if (node_label %in% rownames(difference_matrix)) {
        # Assign the value from the difference matrix to this node
        values_vector[node] <- difference_matrix[node_label, node_label]
        # Mark the label as assigned
        assigned_labels <- c(assigned_labels, node_label)
      }
    }
    
    # Assign the unique label to the node in the tree
    exon4_time$node.label[node - Ntip(exon4_time)] <- node_label
  }
  else if (length(descendant_labels) > 2) {
    # If there are more than two descendants, use the first and last species to form the node label
    node_label <- paste(descendant_labels[1], descendant_labels[length(descendant_labels)], sep = "_")
    
    # Assign the label to the node in the tree
    exon4_time$node.label[node - Ntip(exon4_time)] <- node_label
  }
}


colorTree <- function(phy, node.vals, col = c("blue", "red"), edge.width=2, cex=1, ...){
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
  
  plot.phylo(x=phy, edge.color = phy$edge.colors, edge.width=edge.width, cex = cex, ...)
}

test <- colorTree(phy = exon4_time, node.vals = values_vector)









# Plot the result
plot(cont_map_result, main = "Phylogenetic Tree with Node Differences Mapped")

plot(exon4_time, cex = 0.6)
nodelabels()


# Example of assigning simple labels to internal nodes
for (node in (Ntip(exon4_time) + 1):(Ntip(exon4_time) + Nnode(exon4_time))) {
  exon4_time$node.label[node - Ntip(exon4_time)] <- paste("Node", node)
}

# Check current node labels
exon4_time$node.label

# Find the index of 'Node 189'
index <- which(exon4_time$node.label == "Node 190")

# Rename it to 'Ischnura_elegans_Achroia_grisella'
exon4_time$node.label[index] <- "Zootermopsis_nevadensis_Achroia_grisella"







plot.phylo(dscam1, type= "fan", cex = 0.6, edge.width = 2, show.tip.label = FALSE)
