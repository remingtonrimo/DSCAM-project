# making a line graphs and distribution plots for hymenopterans
library(ggplot2)
library(scales)
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
  labs(x = "Length of DSCAM", y = "phyloP score", title = "A") +
  scale_color_manual(values = c("#009532","#257A4D","#4A5F67","#704582","#952A9C","#BA0FB7")) +  # Specify line colors
  theme_minimal() +
  facet_wrap(~ group, ncol = 1)

# make distribution plots for the values for each line
setwd("C:/Users/remin/Desktop/School/Programs/BEAST.v.2.7.6/DSCAM/phyloP/eus/analysis")
hymen1<-read.csv("hymen300-400.csv")
hymen2<-read.csv("hymen1600-1700.csv")

density_plot1 <- ggplot(hymen1, aes(x = phyloP.score, color = Group)) +
  geom_density(show.legend = FALSE) +
  facet_wrap(~ Group, ncol=1, scales="free_y") +
  labs(title = "300-400 nt",
       x = "phyloP Score",
       y = "Density") +
  scale_color_manual(values = c("#009532","#257A4D","#4A5F67","#704582","#952A9C","#BA0FB7")) +  # Specify line colors
  scale_y_continuous(breaks = pretty_breaks(n = 5), labels = function(x) floor(x)) +
  theme_minimal()

histogram1 <- ggplot(hymen1, aes(x = phyloP.score, fill = Group)) +
  geom_histogram(position = "dodge", bins = 50, show.legend = FALSE) +
  facet_wrap(~ Group, ncol = 1, scales = "free_y") +
  labs(title = "B",
       x = "phyloP Score",
       y = "Count") +
  scale_fill_manual(values = c("#009532", "#257A4D", "#4A5F67", "#704582", "#952A9C", "#BA0FB7")) +  # Specify fill colors
  scale_y_continuous(breaks = pretty_breaks(n = 5), labels = function(x) floor(x)) +
  theme_minimal()

density_plot2 <- ggplot(hymen2, aes(x = phyloP.score, color = Group)) +
  geom_density(show.legend = FALSE) +
  facet_wrap(~ Group, ncol=1,scales = "free_y") +
  labs(title = "1600-1700 nt",
       x = "phyloP Score",
       y = "") +
  scale_color_manual(values = c("#009532","#257A4D","#4A5F67","#704582","#952A9C","#BA0FB7")) +  # Specify line colors
  scale_y_continuous(breaks = pretty) +
  theme_minimal()

histogram2 <- ggplot(hymen2, aes(x = phyloP.score, fill = Group)) +
  geom_histogram(position = "dodge", bins = 50, show.legend = FALSE) +
  facet_wrap(~ Group, ncol = 1, scales = "free_y") +
  labs(title = "C",
       x = "phyloP Score",
       y = "") +
  scale_fill_manual(values = c("#009532", "#257A4D", "#4A5F67", "#704582", "#952A9C", "#BA0FB7")) +  # Specify fill colors
  scale_y_continuous(breaks = pretty_breaks(n = 5), labels = function(x) floor(x)) +
  theme_minimal()

# Combining the figures together
library(ggpubr)

obj <- ggarrange(histogram1, histogram2, ncol = 2, nrow = 1)

ggarrange(line1, obj, ncol=2, nrow=1,widths = c(2.5, 2))

