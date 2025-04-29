library(ape)

# doing a dnds analysis
setwd("~/dscam")
  csv <- read.csv("dscam.csv")
setwd("~/dscam/dnds_analysis/")
  align <- read.FASTA("align_dscam.fasta")
  
# split up 'align' into sociality catagories and compare between them
  csv$Species <- gsub(" ", "_", csv$Species)
  names(align) <- gsub("[\t ]", "", names(align))
  eusocial <- csv$Species[csv$Sociality == "Eusocial"]
  presocial <- csv$Species[csv$Sociality == "Presocial"]
  solitary <- csv$Species[csv$Sociality == "Solitary"]
  
  eusocial_align <- align[names(align) %in% eusocial]
  presocial_align <- align[names(align) %in% presocial]
  solitary_align <- align[names(align) %in% solitary]
  
  eus_align_matrix <- as.character(as.matrix(eusocial_align))
# doing the dN/dS analysis
eus_dnds<-dnds(eus_align_matrix, code = 5, codonstart = 3, quiet = FALSE,
     details = FALSE, return.categories = FALSE)
pre_dnds<-dnds(presocial_align, code = 1, codonstart = 3, quiet = FALSE,
     details = FALSE, return.categories = FALSE)
sol_dnds<-dnds(solitary_align, code = 1, codonstart = 3, quiet = FALSE,
     details = FALSE, return.categories = FALSE)
# It looks like the alignment has too many gaps for a dN/dS analysis to be conducted

# compare them via bar graphs









#testing how this works

data(woodmouse)
test<-dnds(woodmouse, code = 2, quiet = TRUE)







