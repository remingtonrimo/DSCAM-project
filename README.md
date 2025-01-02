# DSCAM-project
## Collect sequences
- First we collect all the data available on NCBI for DSCAM1
- Grab the entire gene region, and then cut the introns out
- See DSCAM-project/parameters/data_collection.R for the script
## Reconstruct/Assess the quality of Phylogeny
- Reconstruct the DSCAM1 gene tree using these parameters
    - GTR+G+I nucleotide substitution model
    - 8 gamma catagories
    - Optimized relaxed clock model
    - Birth-Death model
    - 50 million generations
- Assess the quality of tree via Tracer.v.2.7.7
- Apply same parameters to exon trees
## Dating phylogeny
- Date the root, Ischnura elegans, on the DSCAM1 phylogeny at 401 mya
- DSCAM-project/*dating_phylogenies.R* script
## Comparing DSCAM1 phylogeny to TimeTree phylogeny
- Create a cophylo using DSCAM-project/chronogram/*FigS2.R* script to detemerine discordant species, remove species from chronogram.
- Overlay phylogenies on top of each other for visualization. Subtract the TimeTree branch lengths from DSCAM1 branch lengths. Use DSCAM-project/chronogram/*Fig1.R* to make a heatmap for these values.
## Compare exon trees to overall gene tree
- 
## phyloP analysis
- Conduct a base by base likelihood ratio test using phyloP using bash. Check DSCAM-project/phyloP/*workflow.txt* for code used.
- Use DSCAM-project/phyloP/*Fig2* & DSCAM-project/phyloP/*Fig3* to analyze the phyloP scores.
