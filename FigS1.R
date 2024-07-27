library(phytools)
### Making density trees for DSCAM and Mitochondrial trees 
# Functions available for densityTree
  ?densityTree
  setwd("c:/Users/remin/Desktop/School/Programs/BEAST.v.2.7.6/DSCAM/densitytree")
# These are density trees
  d1trees<-read.nexus("rem2_4-29-24-10-25_mafft_miss-rem.trees")
  mitotrees<-read.nexus("mito-5-2.trees")
  
  random.trees<-sample(d1trees,size=100)
  random.trees01<-sample(mitotrees,size=100)
  
# writing out the random trees from the original trees files in order to make a density tree
# write.tree(random.trees, file = "dscam.tree")
# write.tree(random.trees01, file = "mito.tree")
dscam<-read.tree("dscam.tree")
mito<-read.tree("mito.tree")
  
  layout(matrix(1:2, 1, 2), widths = c(3, 3))
  densityTree(dscam,use.gradient=TRUE, alpha=0.01, fix.depth=TRUE, 
              compute.consesus=FALSE, show.axis=FALSE,fsize=0.7)
  title(main = "DSCAM")
# dscamplot <- recordPlot()
  
  densityTree(mito, use.gradient=TRUE, alpha=0.01, fix.depth=TRUE, 
              compute.consesus=FALSE, show.axis=FALSE,fsize=0.7)
  title(main = "Mitochondria")
# mitoplot <- recordPlot() 
  
# making a line graphs and distribution plots for hymenopterans
  library(ggplot2)
# Set the working directory
  setwd("C:/Users/remin/Desktop/School/Programs/BEAST.v.2.7.6/DSCAM/phyloP/eus/analysis/csv")
# read in the wig files  
  ants<-read.table("ants.wig")
  eus_bees<-read.table("bees.wig")
  sol_bees<-read.table("solbees.wig")
  parasitoids<-read.table("parasites.wig")
  wasps<-read.table("wasps.wig")
  termites<-read.table("termites.wig")
# make the data frames  
  df1 <- data.frame(x = 1:4755, y = ants$V1, group = "Ants")
  df2 <- data.frame(x = 1:4755, y = eus_bees$V1, group = "Eusocial Bees")
  df3 <- data.frame(x = 1:4755, y = sol_bees$V1, group = "Solitary Bees")
  df4 <- data.frame(x = 1:4755, y = parasitoids$V1, group = "Parasitoid Wasps")
  df5 <- data.frame(x = 1:4755, y = wasps$V1, group = "Vespid Waps")
  df6 <- data.frame(x = 1:4755, y = termites$V1, group = "Termites")
# combine the data frames  
  combined_df <- rbind(df1, df2, df3, df4, df5, df6)
# plot the lines and make them an object  
  line1<-ggplot(combined_df, aes(x = x, y = y, color = group)) +
    geom_line(aes(group = group), linewidth = 0.5, color = "black", show.legend = FALSE) +  # Outline
    geom_line(aes(group = group), linewidth = 0.4, show.legend = FALSE) +  # Main line
    geom_rect(aes(xmin = 250, xmax = 450, ymin = -2, ymax = 2),
              fill = NA, color = "red", linetype = "solid", size = 0.4, alpha=0.5) +
    geom_rect(aes(xmin = 1575, xmax = 1700, ymin = -2, ymax = 2),
              fill = NA, color = "red", linetype = "solid", size = 0.4, alpha=0.5) +
    labs(x = "Length of DSCAM", y = "phyloP score", title = "") +
    scale_color_manual(values = c("#009532","#257A4D","#4A5F67","#704582","#952A9C","#BA0FB7")) +  # Specify line colors
    theme_minimal() +
    facet_wrap(~ group, ncol = 1)
  
# make distribution plots for the values for each line
  setwd("C:/Users/remin/Desktop/School/Programs/BEAST.v.2.7.6/DSCAM/phyloP/eus/analysis")
  hymen1<-read.csv("hymen300-400.csv")
  hymen2<-read.csv("hymen1600-1700.csv")

  density_plot1 <- ggplot(hymen1, aes(x = phyloP.score, color = Group)) +
    geom_density(show.legend = FALSE) +
    facet_wrap(~ Group, ncol=1,scales = "free_y") +
    labs(title = "300-400 nt",
               x = "phyloP Score",
               y = "Density") +
    scale_color_manual(values = c("#009532","#257A4D","#4A5F67","#704582","#952A9C","#BA0FB7")) +  # Specify line colors
    theme_minimal()
  
  density_plot2 <- ggplot(hymen2, aes(x = phyloP.score, color = Group)) +
    geom_density(show.legend = FALSE) +
    facet_wrap(~ Group, ncol=1,scales = "free_y") +
    labs(title = "1600-1700 nt",
         x = "phyloP Score",
         y = "") +
    scale_color_manual(values = c("#009532","#257A4D","#4A5F67","#704582","#952A9C","#BA0FB7")) +  # Specify line colors
    theme_minimal()

# Combining the figures together
  library(ggpubr)
  
  ggarrange(line1, density_plot1, density_plot2,ncol=3, nrow=1,widths = c(2.5, 1, 1))
  
  
    
  