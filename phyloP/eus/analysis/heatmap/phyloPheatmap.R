# setting the working directory
  setwd("C:/Users/remin/Desktop/School/Programs/BEAST.v.2.7.6/DSCAM/phyloP/eus/analysis/heatmap")
# loading in the packages
  library(ape)
  library(gplots)
  library(ggplot2)
# read in the csv file
csv<-read.csv("hymen-phyloP.csv")
data<-read.csv("hpP300-400.csv")
m<-as.matrix(csv[,-1])

rownames(m) <- c("Ants", "Eusocial Wasps", "Eusocial Bees", "Solitary Bees", "Parasitoid Wasps", "Termites")
colnames(m) <- 1:4779
colnames(m) <- 300:400

viridis <- colorRampPalette(c("#440154", "#3b528b", "#21918c", "#5ec962","#fde725"))(100)
viridis2 <- colorRampPalette(c("#fde725", "#44bf70","#414487","#440154","#301934","black"))(100)
bluered <- colorRampPalette(c("#0000CC","#400099","#800066","#bf0033","#FF0000"))(100)


heatmap.2(m,
          col = bluered,
          symm = TRUE,
          cexRow = 1,
          cexCol = 0.3,
          Colv = NA,
          Rowv = NA,
          key = FALSE,  # Display the legend
          margins = c(5, 8),# Adjust margins
          lhei = c(0.1,0.5),
          lwid = c(0.2,25),
          trace = "none"  # Remove cell labeling
)


# making a distribution plot for the two different areas
# Density plots with semi-transparent fill
setwd("C:/Users/remin/Desktop/School/Programs/BEAST.v.2.7.6/DSCAM/phyloP/eus/analysis")

dat34<-read.csv("hymen300-400.csv")
dat1617<-read.csv("hymen1600-1700.csv")

viridis <- c("Ants" = "#7ad151", 
                    "Wasps" = "#fde725", 
                    "Parasitic Wasps" = "#440154",
                    "Solitary Bees" = "#414487",
                    "Termites" = "#2a788e",
                    "Bees" = "#22a884")

ggplot(dat34, aes(x=phyloP.score, colour=Group)) + 
  geom_density() +
  scale_colour_manual(values = viridis) +
  theme_minimal()

ggplot(dat1617, aes(x=phyloP.score, colour=Group)) + 
  geom_density() +
  scale_colour_manual(values = viridis) +
  theme_minimal()
