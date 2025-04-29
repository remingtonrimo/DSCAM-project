library(deeptime)
library(ggplot2)
library(dplyr)
library(ggtree)
library(phytools)
library(rphylopic)
# reading in tree
setwd("~/dscam")
dscam1 <- read.newick("iqtree_run/dscam1.new")
plot.phylo(dscam1,type = "fan",cex=0.5)
nodelabels()

get_uuid(name = "Aphididae")

clade_data <- data.frame(
  node = c(378, 305, 201, 264, 289, 366),  
  annote = c("Orthoptera", "Hymenoptera", "Diptera","Lepidoptera", "Coleoptera", "Hemiptera"),
  image=c("0de35750-d9ba-472a-b4f2-3c630777fcc3","10daab45-21e7-4d2a-97e0-9809e0a3eb5b",
          "0cd6cc9f-683c-470e-a4a6-3b68beb826fa","4207925a-1a59-47cd-acc5-e75da4e5e128",
          "81db0738-e72b-472a-9091-91a2082244e2","cd57935f-334e-4067-8078-3c0dc1f44bf0"),
  color = c("#9e103a", "#f9c000","#3a51c2", "#50276e", "#0b7d39", "#e86427")
)

p <- revts(ggtree(dscam1)) +
  coord_geo_radial(
    dat = list("stages", "periods"), alpha = .5, lty = "dashed",
    prop = list(0.66, .34), start = 2 * pi, end = 1.75 * pi, direction = 1,
  ) +
  scale_x_continuous(breaks = seq(-401, 0, 50), labels = abs(seq(-401, 0, 50)),
                     expand = expansion(mult = c(0.05, 0))) +
  scale_y_continuous(guide = NULL, expand = expansion(mult = c(0.01, 0.01))) +
  theme_classic()

# Add clade labels with proper formatting
p + geom_cladelab(
  data = clade_data,
  mapping = aes(node = node, label = annote, image=image, color=annote),
  color = clade_data$color, 
  angle = "auto",
  horizontal = FALSE,
  offset = 5,                
  offset.text = 1,          
  barsize = 2.5,             
  fontsize = 0,            
  extend = 0.3,              
  align = TRUE,              
  show.legend = TRUE        
) +
  scale_color_manual(values = setNames(clade_data$color, clade_data$annote),name = "Order") +
  theme(legend.position = "right",
        legend.background = element_rect(
          color = "black",  
          fill = "white",     
          linewidth = 0.5,    
          linetype = "solid"
        ))

imgs <- sapply(clade_data$image, function(uuid) get_phylopic(uuid = uuid))

# Add the images via GIMP
plot(x=1,y=1, type = "n", ann= FALSE)
add_phylopic_base(img = imgs$`0de35750-d9ba-472a-b4f2-3c630777fcc3`, x = 0.65, y = 0.65, height = 0.1, fill = "#9e103a") 
add_phylopic_base(img = imgs$`10daab45-21e7-4d2a-97e0-9809e0a3eb5b`, x = 0.7, y = 0.7, height = 0.1, fill ="#f9c000") 
add_phylopic_base(img = imgs$`0cd6cc9f-683c-470e-a4a6-3b68beb826fa`, x = 0.8, y = 0.8, height = 0.1, fill ="#3a51c2") 
add_phylopic_base(img = imgs$`4207925a-1a59-47cd-acc5-e75da4e5e128`, x = 0.9, y = 0.9, height = 0.1, fill ="#50276e") 
add_phylopic_base(img = imgs$`81db0738-e72b-472a-9091-91a2082244e2`, x = 1, y = 1, height = 0.1, fill ="#0b7d39") 
add_phylopic_base(img = imgs$`cd57935f-334e-4067-8078-3c0dc1f44bf0`, x = 1.1, y = 1.1, height = 0.1, fill ="#e86427") 


