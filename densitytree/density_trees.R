library(phytools)

#Functions available for densityTree
?densityTree
setwd("c:/Users/remin/Desktop/School/Programs/BEAST.v.2.7.6/DSCAM/dscamrun2")
setwd("c:/Users/remin/Desktop/School/Programs/BEAST.v.2.7.6/DSCAM/mitochondria")
#These are density trees
d1trees<-read.nexus("rem2_4-29-24-10-25_mafft_miss-rem.trees")
mitotrees<-read.nexus("mito-5-2.trees")
#try density tree function. 
#densityTree(mitotrees)
#root the tree on octopus 
# newmito<-root(mitotrees, "Octopus_sinensis") This is for PCDH files
#make a subset of trees, because it's a large file

random.trees<-sample(d1trees,size=100)
random.trees01<-sample(mitotrees,size=100)

# writing out the random trees from the original trees files in order to make a density tree
setwd("c:/Users/remin/Desktop/School/Programs/BEAST.v.2.7.6/DSCAM/densitytree")
#write.tree(random.trees, file = "dscam.tree")
#write.tree(random.trees01, file = "mito.tree")
dscam<-read.tree("dscam.tree")
mito<-read.tree("mito.tree")

# We use the written out tree files because they ladderize better 
# than the other sampled trees
layout(matrix(1:2, 1, 2), widths = c(3, 3))
densityTree(dscam,use.gradient=TRUE, alpha=0.01, fix.depth=TRUE, 
            compute.consesus=FALSE, show.axis=FALSE,fsize=0.7)
title(main = "DSCAM")
#dscamplot <- recordPlot()

densityTree(mito, use.gradient=TRUE, alpha=0.01, fix.depth=TRUE, 
            compute.consesus=FALSE, show.axis=FALSE,fsize=0.7)
title(main = "Mitochondria")
#mitoplot <- recordPlot()

