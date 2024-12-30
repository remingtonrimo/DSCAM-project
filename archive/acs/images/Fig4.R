library(phytools)

##Reading in tree and data 
##make sure to switch your directory tree is in BEAST v1.10.4
setwd("C:/Users/remin/Desktop/School/Programs/BEAST.v.2.7.6/DSCAM/ancestralcharacterstate")
DSCAM<-read.newick("4-29dscamrun2.new")
plot(DSCAM, cex=0.6)
##csv is in Excel Folders
setwd("C:/Users/remin/Desktop/School/Programs/BEAST.v1.10.4/BEAST v1.10.4/dscam-may-10-23/A-final_files/sociality//")
dat<-read.csv("dscam-soc02.csv")
str(dat)
##Getting Sociality and binding it to species name
dat$sociality<-as.character(dat$Sociality)
Sociality<-dat$Sociality
names(Sociality)<-dat$Species

print(DSCAM$tip.label)
print(dat$ï..Species)

##To plot the simmap
trees<-make.simmap(DSCAM, x=Sociality, model="ER", nsim=100, message=FALSE)
#Changing the color
##testcol line is important!
col<-c("Eusocial"="black","Solitary"="black","Presocial"="black")
testcol<-c("Eusocial"="#fde725", "Presocial"="#21918c","Solitary"="#440154")

plotSimmap(trees[[1]], colors=col,fsize=0.7, ftype="i")
testsim<-describe.simmap(trees, plot=FALSE)
nodelabels(pie=testsim$ace, piecol=testcol, cex=0.4)
#edgelabels(text="fuck",tip=1:66)
add.simmap.legend(colors=testcol, prompt=FALSE, x = -1.0*par()$usr[1], y=-10.0*par()$usr[3], fsize=0.9)
title(main = "A")

final_plot <- recordPlot()
final_plot

## testing cool things to do with the figures
test<-ggarrange(line1, line2, seqlogo,ncol=1, nrow=3,heights = c(2.5, 1.3, 1))

install.packages("gridGraphics")
ggarrange(final_plot, test, ncol=2, nrow=1)

# repeat this for mitochondria tree
setwd("C:/Users/remin/Desktop/School/Programs/BEAST.v.2.7.6/DSCAM/ancestralcharacterstate")
mito<-read.newick("mito5-2.new")

# unnecessary code
# testing adding clade labels
nodelabels()
cladelabels(tree=NULL, text="Hymenoptera", node=96, offset=8, wing.length=NULL, cex=1,orientation="horizontal")
cladelabels(tree=NULL, text="Diptera", node=75, offset=8, wing.length=NULL, cex=1,orientation="horizontal")
cladelabels(tree=NULL, text="Coleoptera", node=85, offset=8, wing.length=NULL, cex=1,orientation="horizontal")
cladelabels(tree=NULL, text="Orthoptera", node=127, offset=8, wing.length=NULL, cex=1,orientation="horizontal")
cladelabels(tree=NULL, text="Blattodea", node=125, offset=8, wing.length=NULL, cex=1,orientation="horizontal")
cladelabels(tree=NULL, text="Odonata", node=130, offset=8, wing.length=NULL, cex=1,orientation="horizontal")
cladelabels(tree=NULL, text="Lepidoptera", node=74, offset=8, wing.length=NULL, cex=1,orientation="horizontal")
