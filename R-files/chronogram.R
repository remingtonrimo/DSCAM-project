# We will be making a cophlo to compare the species tree (TimeTree) and the DSCAM1 tree
# This is so we can compare the branch heights to each other
# First we load in the neccessary packages
library(phytools)
library(ape)
library(gplots)

# We read in the trees
setwd("~/dscam/iqtree_run")
dscam1 <- read.newick("dscam1.new")
setwd("~/dscam/control_tree")
control <- read.newick("control_tree.nwk")

# dscam vs. control tree
# Extract tip labels
tips_tree3 <- dscam$tip.label
tips_tree2 <- control$tip.label
# Find tips unique to each tree
unique_to_tree1 <- setdiff(tips_tree3, tips_tree2)
unique_to_tree2 <- setdiff(tips_tree2, tips_tree3)

# Drop unique tips from each tree
pruned_dscam <- drop.tip(dscam, unique_to_tree1)
pruned_control <- drop.tip(control, unique_to_tree2)

# make the cophylo
cophylo <- cophylo(pruned_dscam,pruned_control)
plot(cophylo, link.type="curved",link.lwd=3,link.lty="solid",
     link.col=make.transparent("blue",0.25),fsize=0.5)

# Ensure both trees have tip labels
if (is.null(t1$tip.label) || is.null(t2$tip.label)) {
  stop("Both trees must have tip labels.")
}

# Find shared tip labels
shared_tips <- intersect(pruned_dscam$tip.label, pruned_control$tip.label)

# Match shared tips to their indices in both trees
dscam_indices <- match(shared_tips, pruned_dscam$tip.label)
control_indices <- match(shared_tips, pruned_control$tip.label)

#### gathering branch lengths for DSCAM and species trees ####
dscamdist<-cophenetic(pruned_dscam)
controldist<-cophenetic(pruned_control)
# dividing distances by 2
dscam_dist2<-dscamdist/2
control_dist2<-controldist/2

# the order for tip labels need to be placed in for the heatmap where matrix_order is
matrix_order<-c("Schistocerca_gregaria","Zootermopsis_nevadensis","Bacillus_rossius_redtenbacheri","Thrips_palmi",
                "Planococcus_citri", "Sipha_flava","Rhopalosiphum_padi", "Melanaphis_sacchari", "Myzus_persicae",
                "Metopolophium_dirhodum","Diuraphis_noxia","Daktulosphaira_vitifoliae", "Adelges_cooleyi", 
                "Bemisia_tabaci", "Onthophagus_taurus", "Tribolium_madens","Tenebrio_molitor","Sitophilus_oryzae",
                "Euwallacea_fornicatus","Cylas_formicarius","Nicrophorus_vespilloides","Chrysoperla_carnea",
                "Vanessa_cardui","Vanessa_atalanta","Nymphalis_io","Melitaea_cinxia","Pararge_aegeria","Pieris_napi",
                "Pieris_brassicae","Zerene_cesonia","Colias_croceus","Leptidea_sinapis","Battus_philenor",
                "Trichoplusia_ni","Spodoptera_litura","Helicoverpa_zea","Bombyx_mandarina","Hyposmocoma_kahamanoa",
                "Cydia_pomonella","Phlebotomus_papatasi","Hermetia_illucens","Episyrphus_balteatus","Glossina_fuscipes",
                "Scaptodrosophila_lebanonensis","Drosophila_tropicalis","Drosophila_serrata","Drosophila_subpulchrella",
                "Drosophila_teissieri","Drosophila_erecta","Drosophila_mauritiana","Drosophila_gunungcola",
                "Drosophila_pseudoobscura","Drosophila_subobscura","Drosophila_guanche","Drosophila_sulfurigaster_albostrigata",
                "Drosophila_nasuta","Drosophila_innubila","Drosophila_virilis","Drosophila_novamexicana",
                "Drosophila_montana","Drosophila_arizonae","Teleopsis_dalmanni","Rhagoletis_pomonella","Ceratitis_capitata",
                "Bactrocera_tryoni","Bactrocera_neohumeralis","Bactrocera_latifrons","Bradysia_coprophila",
                "Uranotaenia_lowii","Culex_quinquefasciatus","Wyeomyia_smithii","Sabethes_cyaneus","Malaya_genurostris",
                "Armigeres_subalbatus","Aedes_aegypti","Anopheles_marshallii","Anopheles_gambiae","Anopheles_arabiensis",
                "Anopheles_darlingi","Culicoides_brevitarsis","Ctenocephalides_felis","Ooceraea_biroi","Pseudomyrmex_gracilis",
                "Pogonomyrmex_barbatus","Vollenhovia_emeryi","Temnothorax_curvispinosus","Wasmannia_auropunctata",
                "Trachymyrmex_zeteki","Trachymyrmex_cornetzi","Trachymyrmex_septentrionalis","Acromyrmex_echinatior",
                "Atta_cephalotes","Cyphomyrmex_costatus","Nylanderia_fulva","Formica_exsecta","Cataglyphis_hispanica",
                "Odontomachus_brunneus","Dinoponera_quadriceps","Megalopta_genalis","Dufourea_novaeangliae",
                "Colletes_gigas","Osmia_lignaria","Megachile_rotundata","Habropoda_laboriosa","Frieseomelitta_varia",
                "Bombus_pascuorum","Bombus_pyrosoma","Bombus_bifarius","Bombus_affinis","Eufriesea_mexicana",
                "Apis_florea","Apis_dorsata","Leptopilina_boulardi","Ceratosolen_solmsi_marchali","Cephus_cinctus",
                "Diprion_similis","Ischnura_elegans")

# Initialize difference_matrix using matrix_order for row and column names
difference_matrix <- matrix(
  NA,
  nrow = length(matrix_order),
  ncol = length(matrix_order),
  dimnames = list(matrix_order, matrix_order)
)

# Iterate through matrix_order to ensure values are aligned correctly
for (row_name in matrix_order) {
  for (col_name in matrix_order) {
    # Check if the names exist in both controldist and dscam_dist2
    if (row_name %in% rownames(control_dist2) && col_name %in% colnames(control_dist2) &&
        row_name %in% rownames(dscam_dist2) && col_name %in% colnames(dscam_dist2)) {
      
      # Perform the subtraction and store it in the appropriate position
      difference_matrix[row_name, col_name] <- control_dist2[row_name, col_name] - dscam_dist2[row_name, col_name]
    }
  }
}


# View the result
difference_matrix

difference_matrix[upper.tri(difference_matrix)] <- NA

####Heat map with a continuous legend  
#The heat map I used
#Color palettes were from https://waldyrious.net/viridis-palette-generator/
viridis <- colorRampPalette(c("#fde725","#5ec962","#31688e","#3b528b","#440154"))(100)
inferno <- colorRampPalette(c("#fbb61a","#ed6925","#bc3754","#781c6d","#320a5e","#000004"))(100)
magma <- colorRampPalette(c("#000004","#51127c","#b73779","#fc8961","#fcfdbf"))(100)
plasma <- colorRampPalette(c("#f0f921", "#f89540", "#cc4778", "#7e03a8","#0d0887"))(100)
bluered <- colorRampPalette(c("#0000CC","#400099","#800066","#bf0033","#FF0000"))(100)
bloodymary <- colorRampPalette(c("#000004","#440154","#68157B","#DD2476","#EE3B53","#FF512F"))(100)
custom<-colorRampPalette(c("#69000C","#7D0025","#DA3500","#F39300"))(100)

heatmap.2(difference_matrix,
          col = viridis,
          symm = TRUE,
          main = "Cophenetic Matrix Heatmap", #Title of graph
          cexRow = 0.6,
          cexCol = 0.6,
          dendrogram = "none",
          Colv = NA,
          Rowv = NA,
          key = TRUE,  # Display the legend
          key.title = "T(control)-T(dscam1)",  # Title of the legend
          key.xlab = "Branch length difference",  # X-label
          key.ylab = "Distribution",# Y-label
          trace = "none",
          mar = c(9, 9)
)


#### testing code for a new function ####
plot.chronogram <- function(t1, t2, rotate = TRUE, assoc = NULL, ...) {
  # Colors for tree visualization
  if (hasArg(colors)) colors <- list(...)$colors
  else colors <- sapply(c("blue", "red"), make.transparent, alpha = 0.5)
  
  if (hasArg(arr.colors)) arr.colors <- list(...)$arr.colors
  else arr.colors <- sapply(c("blue", "red"), make.transparent, alpha = 0.9)
  
  # Check if input trees are of class "phylo"
  if (!inherits(t1, "phylo") || !inherits(t2, "phylo")) 
    stop("t1 & t2 should be objects of class \"phylo\".")
  
  # Create assoc if not provided
  if (is.null(assoc)) {
    assoc <- intersect(t1$tip.label, t2$tip.label)
    assoc <- if (length(assoc) > 0) 
      cbind(assoc, assoc)
    else NULL
    if (is.null(assoc)) {
      cat("No associations provided or found.\n")
      rotate <- FALSE
    }
  }
  
  # Filter assoc to ensure species are in both trees
  ii <- sapply(assoc[, 1], "%in%", t1$tip.label)
  if (any(!ii)) {
    assoc <- assoc[ii, ]
    cat("Some species in assoc[,1] not in t1. Removing species & links.\n")
  }
  ii <- sapply(assoc[, 2], "%in%", t2$tip.label)
  if (any(!ii)) {
    assoc <- assoc[ii, ]
    cat("Some species in assoc[,2] not in t2. Removing species & links.\n")
  }
  
  # Optimize rotation if requested
  if (rotate) {
    cat("Rotating nodes to optimize matching...\n")
    flush.console()
    
    # Rotate tree 1 based on association with tree 2
    x <- setNames(sapply(assoc[, 2], match, table = t2$tip.label), assoc[, 1])
    t1 <- tipRotate(t1, x * Ntip(t1)/Ntip(t2), right = t2, assoc = assoc, ...)
    
    # Rotate tree 2 based on association with tree 1
    x <- setNames(sapply(assoc[, 1], match, table = t1$tip.label), assoc[, 2])
    t2 <- tipRotate(t2, x * Ntip(t2)/Ntip(t1), left = t1, assoc = assoc, ...)
    
    # Iterative rotation until further improvement is minimal
    best.t1 <- Inf
    best.t2 <- Inf
    while ((best.t2 - attr(t2, "minRotate")) > 0 || (best.t1 - attr(t1, "minRotate")) > 0) {
      best.t1 <- attr(t1, "minRotate")
      best.t2 <- attr(t2, "minRotate")
      
      x <- setNames(sapply(assoc[, 2], match, table = t2$tip.label), assoc[, 1])
      t1 <- tipRotate(t1, x * Ntip(t1)/Ntip(t2), right = t2, assoc = assoc, ...)
      
      x <- setNames(sapply(assoc[, 1], match, table = t1$tip.label), assoc[, 2])
      t2 <- tipRotate(t2, x * Ntip(t2)/Ntip(t1), left = t1, assoc = assoc, ...)
    }
    cat("Rotation optimization done.\n")
  }
  
  # Get node heights for alignment
  h1 <- sapply(1:Ntip(t1), nodeheight, tree = t1)
  h2 <- sapply(1:Ntip(t2), nodeheight, tree = t2)
  
  # Set plotting area for overlapping trees
  plotTree(if (max(h1) > max(h2)) t1 else t2, plot = FALSE,
           mar = c(4.1, 1.1, 1.1, 0.1), direction = "leftwards")
  
  xlim <- get("last_plot.phylo", envir = .PlotPhyloEnv)$x.lim[2:1]
  
  # Plot first tree
  par(fg = "transparent")
  plotTree(t1, color = colors[1], mar = c(4.1, 1.1, 1.1, 0.1),
           xlim = xlim, direction = "leftwards", lwd = 3)
  T1 <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  
  # Overlay second tree
  par(fg = "transparent")
  plotTree(t2, color = colors[2], mar = c(4.1, 1.1, 1.1, 0.1),
           xlim = xlim, add = TRUE, direction = "leftwards", ftype = "off", lwd = 3)
  T2 <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  par(fg = "black")
  
  # Draw lines and species labels between matching nodes
  for (i in 1:Ntip(t1)) {
    arrows(T1$xx[i], T1$yy[i], T2$xx[i], T2$yy[i], lwd = 2,
           col = if (T1$xx[i] > T2$xx[i]) arr.colors[2] else arr.colors[1],
           length = 0.1)
    
    # Add species names near each arrow in italics
    text(mean(c(T1$xx[i], T2$xx[i])), mean(c(T1$yy[i], T2$yy[i])),
         labels = bquote(italic(.(t1$tip.label[i]))), pos = 4, cex = 0.8, col = "black")
    
  }
  
  # Add tip labels and connecting lines
  h <- mapply(function(x, y) if (x < y) x else y, x = T1$xx[1:Ntip(t1)], y = T2$xx[1:Ntip(t2)])
  text(rep(min(h), Ntip(t1)), T1$yy[1:Ntip(t1)], labels = t1$tip.label,
       font = 3, pos = 4, offset = 0.1 * max(c(h1, h2)))
  
  for (i in 1:Ntip(t1)) {
    lines(c(h[i] + if (h[i] > min(h)) 0.005 * diff(xlim) else 0, min(h)),
          rep(T1$yy[i], 2), lty = "dotted")
  }
}

plot.chronogram(pruned_control,pruned_dscam)
plot(pruned_control)


