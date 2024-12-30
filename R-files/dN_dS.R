library(ape)

# doing a dnds analysis
setwd("~/dscam")
  csv <- read.csv("dscam.csv")
setwd("~/dscam/dnds_analysis/")
  align <- read.dna("dscam_mafft.")
# split up 'align' into sociality catagories and compare between them
  eusocial <- csv$Species[csv$Sociality == "Eusocial"]
  presocial <- csv$Species[csv$Sociality == "Presocial"]
  solitary <- csv$Species[csv$Sociality == "Solitary"]
  
  
  
# doing the dnds analysis
dnds(eusocial, code = 1, codonstart = 1, quiet = FALSE,
     details = FALSE, return.categories = FALSE)
dnds(presocial, code = 1, codonstart = 1, quiet = FALSE,
     details = FALSE, return.categories = FALSE)
dnds(solitary, code = 1, codonstart = 1, quiet = FALSE,
     details = FALSE, return.categories = FALSE)

# compare them via bar graphs