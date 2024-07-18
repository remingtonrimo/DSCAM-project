# Making a chronogram manually
library(phytools)

# First we have to load in the time trees
setwd("C:/Users/remin/Desktop/School/Programs/BEAST.v.2.7.6/DSCAM/dscamrun2")
dscam<-read.newick("4-29dscamrun2.new")
setwd("C:/Users/remin/Desktop/School/Programs/BEAST.v.2.7.6/DSCAM/mitochondria")
mito<-read.newick("mito5-2.new")

# We will make a chronogram first
d1.tax<-dscam$tip.label
mito.tax<-mito$tip.label
taxa<-cbind(d1.tax, mito.tax)

# Making a csv file with the taxa from both phylogenies
setwd("C:/Users/remin/Desktop/School/Programs/BEAST.v.2.7.6/DSCAM/chronogram")
write.csv(taxa, "tax.csv")

# Before reading in the csv, manually check the names so they match
allignedtax<-read.csv("tax.csv",header=FALSE)
assoc<-as.matrix(allignedtax)

obj<-cophylo(dscam,mito,rotate=TRUE,assoc=assoc)
plot(obj)  

plot(obj,link.type="curved",link.lwd=3,link.lty="solid",
     link.col=make.transparent("#750647",0.50),fsize=0.5)
mtext("DSCAM (left) vs. Mito (right)", side = 3, line = 4, font = 2, cex =0.8)


# Second we drop the tips of the tree with conflicting topology
  print(dscam$tip.label) #This give us a list of our tips if we need it

plot(dscam)
todrop_D<-c("Bicyclus_anynana","Drosophila_elegans","Drosophila_ficusphila",
            "Leptinotarsa_decemlineata","Nomia_melanderi","Onthophagus_taurus",
            "Thrips_palmi","Helicoverpa_zea","Chrysoperla_carnea","Cotesia_glomerata",
            "Vespa_mandarinia","Cephus_cinctus","Neodiprion_lecontei",
            "Sitophilus_oryzae","Photinus_pyralis","Schistocera_cancellata",
            "Tanypteryx_hageni*","Drosophila_sechellia") # This is a list of species we want to drop
dscamdrop<-drop.tip(dscam, todrop_D) # This drops the tips we made list for previously
plot(dscamdrop, cex=0.6) # This plots the resulting phylogeny with a character expansion of 0.6

# We need to repeat these steps for mitochondrion
print(mito$tip.label)  
todrop_M<-c("Drosophila_ananassae","Drosophila_biarmipes","Drosophila_miranda",
            "Drosophila_santomea","Homalodisca_vitripennis","Macrosteles_quadrilineatus",
            "Periplaneta_americana","Platycnemis_pennipes","Reticulitermes_speratus",
            "Helicoverpa_zea","Chrysoperla_carnea","Cotesia_vestalis",
            "Vespa_mandarinia","Cephus_cinctus","Neodiprion_lecontei",
            "Sitophilus_oryzae","Photinus_pyralis","Schistocerca_cancellata",
            "Tanypteryx_hageni","Drosophila_sechellia")
mitodrop<-drop.tip(mito, todrop_M)
plot(mitodrop, cex=0.6)
# Finally we write those into new trees in order to compare node height
write.tree(mitodrop, file="mitodrop02.new")
write.tree(dscamdrop, file="dscamdrop02.new")


# Then we make a heat map
library(ape)
library(gplots)
#setting the working directory for the correct time treees

####gathering branch lengths for DSCAM and mitochondria
dscamdist<-cophenetic(dscamdrop)
mitodist<-cophenetic(mitodrop)
#dividing distances by 2
dscam_dist2<-dscamdist/2
mito_dist2<-mitodist/2
####writing out as cv
write.csv(dscam_dist2, "dscam_coph.csv")
write.csv(mito_dist2, "mito_coph.csv")
####reading in csv file
mito.dscam_dist<-read.csv("mito-dscam.csv")
mito.dscam_dist[upper.tri(mito.dscam_dist)] <- NA
####makinng a heatmap with mito-dscam
m<-as.matrix(mito.dscam_dist[,-1])
rownames(m) <- mito.dscam_dist$ï..Species
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
viridis <- colorRampPalette(c("#440154", "#3b528b", "#21918c", "#5ec962"))(100)
inferno <- colorRampPalette(c("#fcffa4", "#f98e09", "#bc3754", "#57106e", "#000004"))(100)
magma <- colorRampPalette(c("#fcfdbf", "#fc8961", "#b73779", "#51127c", "#000004"))(100)
plasma <- colorRampPalette(c("#f0f921", "#f89540", "#cc4778", "#7e03a8", "#0d0887"))(100)

heatmap.2(m,
          col = inferno,
          symm = TRUE,
          main = "Cophenetic Matrix Heatmap", #Title of graph
          cexRow = 0.6,
          cexCol = 0.6,
          Colv = NA,
          Rowv = NA,
          key = TRUE,  # Display the legend
          key.title = "T(m)-T(d)",  # Title of the legend
          key.xlab = "Branch length difference",  # X-label
          key.ylab = "Distribution"  # Y-label
)


