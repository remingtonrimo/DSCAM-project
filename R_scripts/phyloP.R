# We are testing for acceleration/conservation in dscam1 using phyloP
#devtools::install_github("CshlSiepelLab/RPHAST")
library(rphast)
library(ape)
library(ggplot2)
library(ggpubr)
library(forcats)
library(scales)
library(dplyr)
# setting the working directory
setwd("~/projects/dscam/phyloP")

# reading the tree & alignment file in
  dscam1 <- readLines("dscam1.new")
  align <- read.msa("dscam1_align.fa", format = "FASTA")
# we fit a phylogenetic model to our alignment file
  # dscam.mod <- phyloFit(msa = align, tree = dscam1, subst.mod = "REV", nrates = 4)
# I ran this through Tree Annotator on my laptop, this likely isn't needed
  # ---- You can skip to this line if you have the model file ----
  dscam_named <- read.tm("dscam_named.mod")

# we read in our csv file so we can categorize our species by sociality
  csv <- read.csv("~/projects/dscam/dscam.csv")
# create vectors for each category
  eusocial <- csv$Species[csv$Sociality == "Eusocial"]
  presocial <- csv$Species[csv$Sociality == "Presocial"]
  solitary <- csv$Species[csv$Sociality == "Solitary"]
# replace the spaces in the species names with "_"
  eusocial<- gsub(" ", "_", eusocial)
  presocial<- gsub(" ", "_", presocial)
# remove species that were dropped from the tree either for missing data or poor alignment
    presocial <- presocial[presocial != "Neodiprion_pinetum"]
  solitary<- gsub(" ", "_", solitary)
# remove species that were dropped from the tree either for missing data or poor alignment
    solitary <- solitary[!solitary %in% c("Calliphora_vicina", "Neocloeon_triangulifer", "Lucilia_sericata", 
                                      "Cloeon_dipterum")]
    
# running likelihood ratio test on all of the subgroups (i.e. sociality)
  eus_phyloP <- phyloP(dscam_named, align, method = "LRT", mode = "CONACC", branches = eusocial)
  pre_phyloP <- phyloP(dscam_named, align, method = "LRT", mode = "CONACC", branches = presocial)
  sol_phyloP <- phyloP(dscam_named, align, method = "LRT", mode = "CONACC", branches = solitary)

# making data frames for phyloP scores 
  eus_df <- data.frame(x=1:8838, y = eus_phyloP$score, group="Eusocial")
  pre_df <- data.frame(x=1:8838, y = pre_phyloP$score, group="Presocial")
  sol_df <- data.frame(x=1:8838, y = sol_phyloP$score, group="Solitary")

# Combining all of the phyloP scores into one data frame
  combined_df <- rbind(eus_df, pre_df, sol_df)

# making a line graph for the phyloP scores
line1<-ggplot(combined_df, aes(x = x, y = y, color = group)) +
  geom_line(aes(group = group), linewidth = 0.2, color = "black", show.legend = FALSE) + 
  geom_line(aes(group = group), linewidth = 0.1) +
  geom_point(aes(group = group), size = 0.01, shape = 22, fill = "white", stroke = 0, show.legend = TRUE) +
  geom_rect(aes(xmin = 1300, xmax = 2200, ymin = -10, ymax = 5),
            fill = NA, color = "red", linetype = "solid", size = 0.4, alpha=0.5) +
  labs(x = "Nucleotide position within DSCAM1", y = "phyloP score", title = "B", color="Sociality") +
  scale_color_manual(values = c("#90d743","#21918c", "#440154")) +
  theme_minimal() +
  theme(legend.position = "right",
        legend.box.background = element_rect(color = "black", size = 1),
  ) +
  guides(color = guide_legend(
    override.aes = list(
      linetype = 0,  
      shape = 22,    
      size = 5,      
      fill = c("#90d743","#21918c", "#440154"), 
      color = "black" 
    ) 
    )) +
  facet_wrap(~ group, ncol = 1)
plot(line1)

# This graph seems to indicate acceleration for exons 4-6, i.g. Ig domains 2 & 3

# Code for Figure 2C
  eus_4 <- eus_df[1300:1761,]
  pre_4 <- pre_df[1300:1761,]
  sol_4 <- sol_df[1300:1761,]
exon4 <- rbind(eus_4, pre_4, sol_4)

plot1 <- ggplot(exon4, aes(x = x, y = y, color = group)) +
  geom_line(aes(group = group), linewidth = 0.6, color = "black", show.legend = FALSE) +
  geom_line(aes(group = group), linewidth = 0.4) +  
  geom_hline(yintercept = 2, color = "red", linetype = "solid", linewidth = 0.3) +
  geom_hline(yintercept = -2, color = "red", linetype = "solid", linewidth = 0.3) +
  labs(x = "exon 4", y = "phyloP score", title = "C") +
  scale_color_manual(values = c("#90d743", "#21918c", "#440154")) +  
  scale_y_continuous(limits = c(-10, 5), breaks = seq(-10, 5, 5)) +
  guides(color = "none") +  
  theme_minimal() +
  theme(
  ) +
  facet_wrap(~ group, ncol = 1)


  eus_5 <- eus_df[1762:1901,]
  pre_5 <- pre_df[1762:1901,]
  sol_5 <- sol_df[1762:1901,]
exon5 <- rbind(eus_5, pre_5, sol_5)

plot2 <- ggplot(exon5, aes(x = x, y = y, color = group)) +
  geom_line(aes(group = group), linewidth = 0.6, color = "black", show.legend = FALSE) +
  geom_line(aes(group = group), linewidth = 0.4) +
  geom_hline(yintercept = 2, color = "red", linetype = "solid", linewidth = 0.3) +
  geom_hline(yintercept = -2, color = "red", linetype = "solid", linewidth = 0.3) +
  labs(x = "exon 5", y="",title = "") +
  scale_color_manual(values = c("#90d743","#21918c", "#440154")) + 
  scale_y_continuous(limits = c(-10, 5), breaks = seq(-10, 5, 5)) +
  guides(color = "none") +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),  
    axis.ticks.y = element_blank()  
  ) +
  facet_wrap(~ group, ncol = 1)

eus_6 <- eus_df[1902:2200,]
pre_6 <- pre_df[1902:2200,]
sol_6 <- sol_df[1902:2200,]
exon6 <- rbind(eus_6, pre_6, sol_6)

plot3 <- ggplot(exon6, aes(x = x, y = y, color = group)) +
  geom_line(aes(group = group), linewidth = 0.6, color = "black", show.legend = FALSE) +
  geom_line(aes(group = group), linewidth = 0.4) +  
  geom_hline(yintercept = 2, color = "red", linetype = "solid", linewidth = 0.3) +
  geom_hline(yintercept = -2, color = "red", linetype = "solid", linewidth = 0.3) +
  labs(x = "exon 6", y = "", title = "") +
  scale_color_manual(values = c("#90d743","#21918c", "#440154")) + 
  scale_y_continuous(limits = c(-10, 5), breaks = seq(-10, 5, 5)) +
  guides(color = "none") +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
    ) +
  facet_wrap(~ group, ncol = 1)

# group the exons line graphs together
exons <- ggarrange(plot1, plot2, plot3,ncol=3, nrow=1,widths = c(0.44, 0.18, 0.38))
# put them all together
ggarrange(line1, exons,ncol=1, nrow=2)

# Order likelihood ratio test
coleoptera <- csv$Species[csv$Order == "Coleoptera"]
diptera <- csv$Species[csv$Order == "Diptera"]
lepidoptera <- csv$Species[csv$Order == "Lepidoptera"]
orthoptera <- csv$Species[csv$Order == "Orthoptera"]
hymenoptera <- csv$Species[csv$Order == "Hymenoptera"]
hemiptera <- csv$Species[csv$Order == "Hemiptera"]

# replace the spaces in the species names with "_"
coleoptera<- gsub(" ", "_", coleoptera)
diptera<- gsub(" ", "_", diptera)
  diptera <- diptera[!diptera %in% c("Calliphora_vicina", "Lucilia_sericata")]
lepidoptera<- gsub(" ", "_", lepidoptera)
orthoptera<- gsub(" ", "_", orthoptera)
hymenoptera<- gsub(" ", "_", hymenoptera)
  hymenoptera <- hymenoptera[hymenoptera != "Neodiprion_pinetum"]
hemiptera<- gsub(" ", "_", hemiptera)

# Define the groups and their names
groups <- list(
  coleoptera = coleoptera,
  diptera = diptera,
  lepidoptera = lepidoptera,
  orthoptera = orthoptera,
  hymenoptera = hymenoptera,
  hemiptera = hemiptera
)

# Initialize an empty list to store results
phyloP_order_results <- list()

# Loop through each group and run phyloP
for (group in names(groups)) {
  phyloP_order_results[[group]] <- phyloP(dscam_named, align, method = "LRT", mode = "CONACC", branches = groups[[group]])
}

# Initialize a data frame to store scores
order_scores <- data.frame()

# Loop through each group's results
for (group in names(phyloP_order_results)) {
  # Extract the score for the current group
  group_scores <- phyloP_order_results[[group]]$score
  
  # Create a data frame with the group name and scores
  group_df <- data.frame(x=1:8838,
    Group = group,
    Score = group_scores
  )
  
  # Combine with the main data frame
  order_scores <- rbind(order_scores, group_df)
}

# Replace "coleoptera" with "Coleoptera" in the Group column
# Not automated, you have to replace each group name with the capitalized version
order_scores$Group <- ifelse(order_scores$Group == "lepidoptera", "Lepidoptera", order_scores$Group)

line_supp<-ggplot(order_scores, aes(x = x, y = Score, color = Group)) +
  geom_line(aes(group = Group), linewidth = 0.2, color = "black", show.legend = FALSE) +  
  geom_line(aes(group = Group), linewidth = 0.1) +  
  geom_rect(aes(xmin = 1400, xmax = 2200, ymin = -15, ymax = 5),
            fill = NA, color = "red", linetype = "solid", size = 0.4, alpha=0.5) +
  labs(x = "Length of DSCAM", y = "phyloP score", title = "phyloP score by Order") +
  scale_color_manual(values = c("#0b7d39", "#3a51c2", "#e86427", "#f9c000", "#50276e", "#9e103a")) +
  theme_minimal() +
  facet_wrap(~ Group, ncol = 1) 

# Indicates that acceleration is specific to Hymenoptera in exons 4-6, i.g. Ig domains 2 & 3
# grabbing the data from group_results for exons 4,5,6
df1 <- order_scores$Coleoptera[1400:2200,]
df1$Group <- "Coleoptera"
df2 <- order_scores$Diptera[1400:2200,]
df2$Group <- "Diptera"
df3 <- order_scores$Hemiptera[1400:2200,]
df3$Group <- "Hemiptera"
df4 <- order_scores$Hymenoptera[1400:2200,]
df4$Group <- "Hymenoptera"
df5 <- order_scores$Lepidoptera[1400:2200,]
df5$Group <- "Lepidoptera"
df6 <- order_scores$Orthoptera[1400:2200,]
df6$Group <- "Orthoptera"
hist_o <- rbind(df1,df2,df3,df4,df5,df6)
# Add a new column categorizing scores
hist_o <- hist_o %>%
  mutate(score_category = case_when(
    score > 2 ~ "Conserved",
    score < -2 ~ "Accelerated",
    TRUE ~ "Neutral"
  ))
# Calculate proportions per Group
prop_data <- hist_o %>%
  group_by(Group, score_category) %>%
  summarise(count = n()) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup()
prop_data <- prop_data %>%
  mutate(Group = fct_rev(fct_inorder(Group)))
# adjust the size settings as need be
bar_stacked <- ggplot(prop_data, aes(x = Group, y = proportion, fill = score_category)) +
  geom_bar(stat = "identity", position = position_stack(reverse = TRUE), alpha = 0.9) +
  coord_flip() +
  labs(
    title = "B",
    x = "",
    y = "Proportion",
    fill = "Category"
  ) +
  scale_fill_manual(values = c(
    "Conserved" = "black",
    "Neutral" = "gray70", 
    "Accelerated" = "#b80707"   
  )) +
  theme_minimal() +
  theme(
    legend.position = "right",
    axis.text.y = element_blank(),
    # Customize legend box:
    legend.background = element_rect(
      color = "black",   
      fill = "white",    
      size = 0.5,        
      linetype = "solid" 
    ),
    # Adjust legend key appearance (optional):
    legend.key = element_rect(
      color = "gray50",  
      size = 0.3
    )
  )

ggarrange(line_supp, bar_stacked, ncol=2, nrow=1 ,widths = c(1, 1))

# Group analysis of Hymenoptera
# Define the groups in a list
groups <- list(
  Parasitoids = csv$Species[csv$Group == "Parasitoid"],
  Apidae_Megachilidae = csv$Species[csv$Group == "Long-tongue"],
  Colletidae_Halictidae = csv$Species[csv$Group == "Short-tongue"],
  Ants = csv$Species[csv$Group == "Ant"],
  Sawfly = csv$Species[csv$Group == "Sawfly"],
  Vespid = csv$Species[csv$Group == "Vespid"]
)

# Replace spaces with underscores in each group using a for loop
groups <- lapply(groups, function(species) gsub(" ", "_", species))

groups$Sawfly <- groups$Sawfly[groups$Sawfly != "Neodiprion_pinetum"]

# Initialize an empty list to store results
group_results <- list()

# Loop through each group and run phyloP
for (group in names(groups)) {
  group_results[[group]] <- phyloP(dscam_named, align, method = "LRT", mode = "CONACC", branches = groups[[group]])
}

# Initialize a data frame to store scores
group_scores <- data.frame()

# Loop through each group's results
for (group in names(group_results)) {
  # Extract the score for the current group
  scores <- group_results[[group]]$score
  
  # Create a data frame with the group name and scores
  group_df <- data.frame(
    x = 1:8838,
    Group = group,
    Score = scores
  )
  # Combine with the main data frame
  group_scores <- rbind(group_scores, group_df)
}

# the entire length of the gene for Hymenoptera
hymen <- ggplot(group_scores, aes(x = x, y = Score, color = Group)) +
  geom_line(aes(group = Group), linewidth = 0.2, color = "black", show.legend = FALSE) +
  geom_line(aes(group = Group), linewidth = 0.1) +
  geom_rect(aes(xmin = 1400, xmax = 2200, ymin = -8, ymax = 5),
            inherit.aes = FALSE, fill = NA, color = "red", linetype = "solid", size = 0.4, alpha = 0.5) +
  labs(x = "Length of DSCAM", y = "phyloP score", title = "DSCAM1 for Hymenoptera") +
  scale_color_manual(values = c("#febd2a", "#f48849", "#db5c68", "#b83289", "#8b0aa5", "#5302a3")) +
  theme_minimal() +
  theme(legend.position = "none") +
  facet_wrap(~ Group, ncol = 1)


# This seems to indicate that eusocial wasps and eusocial bees experienced acceleration
# grabbing the data from group_results for exons 4,5,6
df1 <- group_results$Ants[1400:2200,]
df1$Group <- "Ants"
df2 <- group_results$Apidae_Megachilidae[1400:2200,]
df2$Group <- "Apidae and Megachilidae"
df3 <- group_results$Colletidae_Halictidae[1400:2200,]
df3$Group <- "Colletidae and Halictidae"
df4 <- group_results$Parasitoids[1400:2200,]
df4$Group <- "Parasitoids"
df5 <- group_results$Sawfly[1400:2200,]
df5$Group <- "Sawflies"
df6 <- group_results$Vespid[1400:2200,]
df6$Group <- "Vespids"
hist1 <- rbind(df1,df2,df3,df4,df5,df6)



hymen_exons<- ggplot(hist1, aes(x = coord, y = score, color = Group)) +
  # Add reference lines first (will appear behind the data)
  geom_hline(yintercept = 2, color = "red", linetype = "solid", linewidth = 0.3) +
  geom_hline(yintercept = -2, color = "red", linetype = "solid", linewidth = 0.3) +
  geom_line(aes(group = Group), linewidth = 0.4, color = "black", show.legend = FALSE) +
  geom_line(aes(group = Group), linewidth = 0.3) +
  labs(x = "", y = "phyloP score", title = "A") +
  scale_color_manual(values = c("#febd2a", "#f48849", "#db5c68", "#b83289", "#8b0aa5", "#5302a3")) + 
  scale_y_continuous(
    limits = c(-7.5, 5),
    breaks = seq(-7.5, 5, 2.5),
    labels = seq(-7.5, 5, 2.5)
  ) +
  guides(color = "none") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),  
    axis.ticks.y = element_line(),        
    strip.text = element_text(size = 10)
  ) +
  facet_wrap(~ Group, ncol = 1)

# Here we will show the proportion of scores that are nuetral, accelerated and conserved 
# since phyloP is the -log10(p value) we used the threshold of pval = 0.01 or phyloP = 2 & -2
# Add a new column categorizing scores
hist1 <- hist1 %>%
  mutate(score_category = case_when(
    score > 2 ~ "Conserved",
    score < -2 ~ "Accelerated",
    TRUE ~ "Neutral"
  ))
# Calculate proportions per Group
prop_data <- hist1 %>%
  group_by(Group, score_category) %>%
  summarise(count = n()) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup()
prop_data <- prop_data %>%
  mutate(Group = fct_rev(fct_inorder(Group)))
bar_stacked <- ggplot(prop_data, aes(x = Group, y = proportion, fill = score_category)) +
  geom_bar(stat = "identity", position = position_stack(reverse = TRUE), alpha = 0.9) +
  coord_flip() +
  labs(
    title = "B",
    x = "",
    y = "Proportion",
    fill = "Category"
  ) +
  scale_fill_manual(values = c(
    "Conserved" = "black",  # Purple for conserved
    "Neutral" = "gray70",  # Gray for neutral
    "Accelerated" = "#b80707"   # Gold for accelerated
  )) +
  theme_minimal() +
  theme(
    legend.position = "right",
    axis.text.y = element_blank(),
    # Customize legend box:
    legend.background = element_rect(
      color = "black",   # Border color
      fill = "white",    # Background fill (optional)
      size = 0.5,        # Border thickness
      linetype = "solid" # Solid line
    ),
    # Adjust legend key appearance (optional):
    legend.key = element_rect(
      color = "gray50",  # Border around legend keys
      size = 0.3
    )
  )

ggarrange(hymen_exons, bar_stacked, ncol=2, nrow=1 ,widths = c(1, 1))


# old visualization code below
violin<-ggplot(hist1, aes(x = Group, y = score, fill = Group)) +
  geom_violin(alpha = 0.9, trim = FALSE, width = 0.8) +  # Adjust width for spacing
  # 2. Box plot (overlaid)
  geom_boxplot(
    width = 0.15,           # Narrow width to fit inside violins
    fill = "white",         # White fill for contrast
    alpha = 0.8,            # Slight transparency
    outlier.shape = NA,     # Hide outliers (already shown in violins)
    show.legend = FALSE) +
  # 3. Optional: Add median markers
  stat_summary(
    fun = median, 
    geom = "point", 
    shape = 18,             # Diamond shape
    size = 3, 
    color = "black", 
    show.legend = FALSE) +
  labs(title = "B",
       x = "",
       y = "Score") +
  scale_fill_manual(values = c("#febd2a", "#f48849", "#db5c68", "#b83289", "#8b0aa5", "#5302a3")) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-labels if needed


# data and plot for exon 4
df1 <- group_results$Ants[1300:1761,]
  df1$Group <- "Ants"
df2 <- group_results$EusocialBee[1300:1761,]
  df2$Group <- "Eusocial Bees"
df3 <- group_results$Parasitoids[1300:1761,]
  df3$Group <- "Parasitoids"
df4 <- group_results$Sawfly[1300:1761,]
  df4$Group <- "Sawflies"
df5 <- group_results$SolitaryBee[1300:1761,]
  df5$Group <- "Solitary Bees"
df6 <- group_results$Wasp[1300:1761,]
  df6$Group <- "Wasps"
hist1 <- rbind(df1,df2,df3,df4,df5,df6)

hist_filtered1 <- hist1 %>% filter(score < 0)

dense1 <-ggplot(hist_filtered1, aes(x = score, fill = Group)) +
  geom_density(alpha = 0.9, show.legend = TRUE) +
  facet_wrap(~ Group, ncol = 1, scales = "fixed") +
  labs(title = "B exon 4",
       x = "",
       y = "Density") +
  scale_fill_manual(values = c("#febd2a", "#f48849", "#db5c68", "#b83289", "#8b0aa5", "#5302a3")) +  # Specify fill colors
  theme_minimal() +
  theme(legend.position = "none")

# data and plot for hymenoptera exon 5

df21 <- group_results$Ants[1762:1901,]
  df21$Group <- "Ants"
df22 <- group_results$EusocialBee[1762:1901,]
  df22$Group <- "Eusocial Bees"
df23 <- group_results$Parasitoids[1762:1901,]
  df23$Group <- "Parasitoids"
df24 <- group_results$Sawfly[1762:1901,]
  df24$Group <- "Sawflies"
df25 <- group_results$SolitaryBee[1762:1901,]
  df25$Group <- "Solitary Bees"
df26 <- group_results$Wasp[1762:1901,]
  df26$Group <- "Wasps"
hist2 <- rbind(df21,df22,df23,df24,df25,df26)

hist_filtered2 <- hist2 %>% filter(score < 0)

dense2 <-ggplot(hist_filtered2, aes(x = score, fill = Group)) +
  geom_density(alpha = 0.9, show.legend = TRUE) +
  facet_wrap(~ Group, ncol = 1, scales = "fixed") +
  labs(title = "exon 5",
       x = "phyloP Scores below 0",
       y = "") +
  scale_fill_manual(values = c("#febd2a", "#f48849", "#db5c68", "#b83289", "#8b0aa5", "#5302a3")) +  # Specify fill colors
  theme_minimal() +
  theme(legend.position = "none")

# data and plot for exon 6

df31 <- group_results$Ants[1902:2200,]
  df31$Group <- "Ants"
df32 <- group_results$EusocialBee[1902:2200,]
  df32$Group <- "Eusocial Bees"
df33 <- group_results$Parasitoids[1902:2200,]
  df33$Group <- "Parasitoids"
df34 <- group_results$Sawfly[1902:2200,]
  df34$Group <- "Sawflies"
df35 <- group_results$SolitaryBee[1902:2200,]
  df35$Group <- "Solitary Bees"
df36 <- group_results$Wasp[1902:2200,]
  df36$Group <- "Wasps"
hist3 <- rbind(df31,df32,df33,df34,df35,df36)

hist_filtered3 <- hist3 %>% filter(score < 0)

dense3 <-ggplot(hist_filtered3, aes(x = score, fill = Group)) +
  geom_density(alpha = 0.9, show.legend = TRUE) +
  facet_wrap(~ Group, ncol = 1, scales = "fixed") +
  labs(title = "exon 6",
       x = "",
       y = "") +
  scale_fill_manual(values = c("#febd2a", "#f48849", "#db5c68", "#b83289", "#8b0aa5", "#5302a3")) +  # Specify fill colors
  theme_minimal()

dense <- ggarrange(dense1, dense2, dense3, ncol=3, nrow=1,widths = c(0.75, 0.75, 1))

ggarrange(hymen, dense, ncol=1, nrow=2 ,heights = c(1, 1))


