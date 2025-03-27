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
- DSCAM-project/R-files/*dating_phylogenies.R* script
## Comparing DSCAM1 phylogeny to TimeTree phylogeny
- Use DSCAM-project/R-files/*chronogram.R* for chronogram analysis between the species and gene trees
## Compare exon trees to overall gene tree
- Use DSCAM-project/R-files/*exon_comparison.R* for chronogram analysis between the exon and gene trees
## phyloP analysis
- See DSCAM-project/R-files/*phyloP.R* to conduct the likelihood base-by-base test.
-  Note: *tree_doctor* function isn't available in the rPHAST package to label the interior nodes of the tree. Use dscam_named.mod for analysis.
