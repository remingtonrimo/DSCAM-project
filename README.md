# DSCAM-project
## Density trees & dating phylogenies
- Collect sequences & construct the phylogeny
- Sample 100 trees from the generated trees file and plot those as a density tree using the DSCAM-project/FigS1.R script
- Date the mitochondrial and DSCAM1 phylogenies using the DSCAM-project/dating_phylogenies.R script
# Comparing mitochondrial and DSCAM1 phylogenies
- Create a cophylo using DSCAM-project/chronogram/FigS2.R script to detemerine discordant species, remove species from chronogram.
- Overlay phylogenies on top of each other for visualization. Subtract mitochondrial branch lengths from DSCAM1 branch lengths. Use DSCAM-project/chronogram/Fig1.R to make a heatmap for these values.
# phyloP analysis
- Conduct a base by base likelihood ratio test using phyloP using bash. Check DSCAM-project/phyloP/workflow for code used.
- Use DSCAM-project/phyloP/Fig2 & DSCAM-project/phyloP/Fig3 to analyze the phyloP scores.
# Ancestral Character State reconstruction
- Use DSCAM-project/phyloP/ancestral_character_states/Fig4 to reconstruct the ancestral character states for sociality.
