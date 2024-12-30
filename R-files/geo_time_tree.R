library(deeptime)
library(ggplot2)
library(dplyr)
BiocManager::install("ggtree")
library(ggtree)

# reading in tree
setwd("~/dscam")
dscam1 <- read.newick("iqtree_run/dscam1.new")

revts(ggtree(dscam1)) +
  coord_geo_radial(
    dat = list("stages", "periods"), alpha = .5, lty = "dashed",
    prop = list(0.66, .34), start = 2 * pi, end = 1.75 * pi, direction = 1,
  ) +
  scale_x_continuous(breaks = seq(-401, 0, 50), labels = abs(seq(-401, 0, 50)),
                     expand = expansion(mult = c(0.05, 0))) +
  scale_y_continuous(guide = NULL, expand = expansion(mult = c(0.01, 0.01))) +
  theme_classic()
