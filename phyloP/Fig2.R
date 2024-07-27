install.packages("ggseqlogo")
library(ggseqlogo)
library(ggplot2)
?geom_logo
# Making a sequence logo 
# Example list of DNA sequences
seq1 <- c("tgcgcgcagtggtgtcccagcactacatcacggaggccgaaaacgagtacgttatccggggcaactcggcggtcatgaagtgcaagatccccagcttcgt",
               "ttcgcgcagttgtttcgcagttctacataaccgaggccgagaacgagtatgtaatacgcggaaactcggccgtcctcaagtgcaagataccgagtttcgt",
               "ttagagcagttgttacacagtactacagtacagaggcagaaaacgaatatgttataaaaggtaattccgtcgtaatgaaatgcaaaataccgagcttcgt",
               "ttagagccgttgtgtcgcagtactatgtgaccgaggccgaaaatgaatacgtcatccgcggaaattcggccgtcatgaaatgcaagatacctagtttcgt",
               "tgagagcagtcgtctcgcaattctatgtgacggaagccgaaaacgaatacgtcattcggggaaactcagccgtaatgaagtgcaagataccgagcttcgt",
               "gtaaaatagttgtgtcgcagtactatgtaaccgaggccgaaaacgaatacgtaatactaggaaatgctgctattatgaaatgcaagattccatcatttgt",
               "tgagggcagttgtttctcaatcttacaacgcgaatgtgatggacgaaagtgtccttaaaggaaacacggccatcttcaagtgccacattcccagttttgt",
               "ttagggccgttgtgtcccaattctacgatatagacgtgaacaagcagtacgtaatcaggggcaacgcggccattttaaaatgcgaaattccctcgtttgt",
               "tcagagcagtcgtgaatcaattctacgagaccagagttactgacgaatttgtcctgaaaggcaacaccgggatcttaaagtgcgtagtaccaagtttcgt",
               "tcagagcagtcgtaaaccagtactacgacacggacgtgaacaaggagtacgcgatcagaggcaacgccgccatactgaagtgtcaggtgccgtcgttcgt",
               "ttcgagcagttgtaaaccagtactacgaagccgaagtcgtttcagagtacgtaattcgaggcaataccgctgttatcaaatgcaatattccgtcgtttgt",
               "tgcgtgcagtggtaaatcagttctataaagcggaaatcctgaccgagtacgttatccgagggaataccgcgattctgaaatgtagcataccgtcttttgt",
               "taagagccgttgtgatgcaaaactatctaccagaaatattaactgagtacgcgattcgtggtaatagtgcaatacttaaatgtagtatacccagttacat",
               "ttcgagcagttgttgcacaaacttatcagccagaaataatgacggaatatgtcattagaggaaacagtgctattctgaagtgtagtataccgagttatat",
               "ttcgcgccgttgtcgcacaaacctatcaacccgaaataatgactgagtacgtcattagaggaaacagcgcaatcttaaaatgtagcataccgagttacat",
               "ttcgagctgttgtgtcacaatattttgaagtacaagtgtacgatgtctttgcgattcgtggtaatgctgcgatattcaagtgtcaggtaccgtcatttgt",
               "tcagggccgtggtccatcaattttacgaaacgagaataacggacgaattcgttctgcgagggaatacggcgacccttaaatgcattgtgcccagctttgt",
               "tcagggccgtcgttcatcagtactatcaatccgaggtaaataacgagtacgttattcgaggaaacgcggcgattctaaaatgtagtatcccaagttttgt",
               "tcagggccgtcgttcatcagtactatcaatccgaggtaaataacgagtacgttattcgaggaaacgcggctattctaaaatgtagtatcccaagttttgt",
               "tcagggccgtcgttcatcagtactatcaatccgaggtaaataacgagtacgttattcgaggaaacgcggctattctaaaatgtagtatcccaagttttgt",
               "tcagggccgtcgtcacccagtactacgaggctgaggtggtttcggagtatgtgattcgcggtaacgctgcgattgtcaagtgcaccatacccagtttcgt",
               "taagagccgtggtacatcaatactatcagtcagaagtcaacaacgaatacgttattcgtggcaatgcggcgattctaaagtgcagcataccgagttttgt",
               "taagagccgtcgtgactcagtattacgaggcagaggtggtatcggaatacgtaataaggggaaatgcagccatcttgaagtgcaccataccgagctttgt",
               "taagagccgttgtgactcagtattacgaggcagaggtggtatcggaatacgtaataaggggaaatgcagccatcttgaagtgcaccataccgagctttgt",
               "taagagccgttgtgactcagtattacgaggcagaggtggtatcggaatatgtgataaggggaaatgcagccatcttgaaatgcaccataccgagcttcgt",
               "taagagccgtggtcgggcagtactttgaagttcaggtttacgaccagttcgcgatacgcggcaacgcggccattttcaaatgtcaagttccttccttcgt",
               "taagagccgttgtgactcagtattacgaggcagaggtagtatcggaatatgtgataaggggaaatgcagccatcttgaagtgcaccataccgagcttcgt",
               "taagagccgtcgtcgcgcagccttaccaaccggaaatcctgacagagtatgtaataaggggtaacagcgcgatcttgaaatgcagcatcccgagttacat",
               "tgagagccgtcgtgactcagtactacgaagcggaggtcgtatccgaatacgtaatacgcggaaacgcagccatcgtgaagtgcaccataccgagcttcgt",
               "tgagagccgtcgtgactcagtactacgaagcggaggtcgtatccgaatacgtaatacgcggaaacgcagccatcgtgaagtgcaccataccgagcttcgt",
               "tgagagccgtcgtgactcagtactacgaagcggaggtcgtttccgaatacgtgatacgcggaaatgcagccatcttgaaatgcaccatacccagcttcgt",
               "tgagagccgtggtggcacaatactacgacacagacgtcaacaaggagtacgcgattcgtggaaacagtgccattttgaagtgcgtcgttccgtcctttgt",
               "tgagagccgtggtggcacagtactacgacacagacgtcaacaaggagtacgcgattcgtggaaacagtgccattttgaagtgcgtcgttccatcctttgt",
               "tgagagccgtggtgcaccagtattatcagtccgaagtcaacaacgagtacgtaattcgcggcaacgcggcgattctcaagtgcagcatacccagctttgt",
               "tgagagctgtcgtcactcagccttacaatccagagatattgactgaatacgttataagaggaaacagcgcgatcttaaaatgcagcattcccagttacat",
               "tgcgagctgtcgtgactcaatattacgaggcggaagttgtatcggagtatgtaatacgtggaaatgcggctattttgaaatgcaccataccgagcttcgt",
               "tgagagctgtcgtcactcagccttacaatccggagatcttgaccgaatacgttataaggggaaacagcgcgatattaaaatgcagcatccccagttacat",
               "taagagccgtcgtgactcaatattatgaggcggaagttgtatcggaatacgtgatacgtggaaatgcggctattttaaaatgcaccataccgagctttgt",
               "tgagagccgtcgtcactcagccatacaatccggagatcttaactgagtacgttataaggggaaacagcgcgatcttaaaatgcagcatccccagttacat",
               "tcagggccgtggtggcgcaatactacgacaccgacgtcaacaaggaatacgccataagggggaacagcgccattctcaaatgtgtcgttccttcgttcgt",
               "tccgagcagttgttcaacaattttacgatactgacgtaaataaggaatatgcaattcgtggtaattcagctgttatgaaatgtgttgtgccatcttttgt",
               "tcagagcagttgtttcacagtactacgttacagaggctgaaaatgagtatgtgatccgtggtaacagtgctgttatgaaatgcaagattcctagtttcgt",
               "tcagagcagttgttcagcagttttaccagactgaagttaacaacgaatatgtgatacgtggaaattcagcagttttgaaatgcagcattccatcatttgt",
               "tcagagcagttgtttcacagtactacgttacagaggctgaaaatgagtatgtgatccgtggtaacagtgctgttatgaaatgcaagattcctagtttcgt",
               "tcagagcagttgttcagcagttttaccagactgaagttaacaacgaatatgtgatacgtggaaattcagcagttttgaaatgcagtattccatcatttgt",
               "ttcgagctgttgtatctcagtactacatcactgaagctgagaacgagtacgttattcgtgggaacagtgctgtaatgaagtgcaagattcccagtttcgt",
               "ttcgagctgttgtttctcaatactacatcaccgaggctgagaacgagtatgtgattcgtgggaacagtgccattatgaagtgcaagattcccagtttcgt",
               "ttagagctgttgtttctcagtactacatcaccgaggctgaaaacgagtatgtgatacgtggaaacagtgctgtaatgaagtgcaagattccaagctttgt",
               "tgagagctgtggtagctcagttctatgttacggaagccgagaatgaatatgttatccgtggcaacagtgcagtgatgaaatgcaagattccaagttttgt",
               "tagaagctgtggtaactcagttctatgtcactgaagcagaaaacgagtatgtgattcgtggcaacagcgcagtgatgaaatgcaagattccgagttttgt",
               "ttagagccgttgtaccccagtcctatacagtcaatgttatggatgaggctgttttgcgtggcaatagtgccattttaaaatgtcatattccttcatttgt",
               "tgcgagccgtcgtttcccaacactatgaagaggatatacacaaggcttttgtcatccgcggcaattccgccatcttgaaatgtgatataccgtcgttcgt",
               "tgcgagccgtcgtttcccaacactacgaagaagatatacacaaggcatttgtcatccgcggcaattccgccatattgaaatgcgatataccgtcgtttgt",
               "tgcgagccgtcgtttcccaacactacgaagaagatatacacaaggcatttgtcatccgcggcaattccgccatattgaaatgcgatataccgtcgtttgt",
               "tgcgagccgtcgtttcccaacactacgaagaagatatacacaaggcatttgtcatccgcggcaattccgccatattgaaatgcgatataccgtcgtttgt",
               "tgcgagccgtcgtttcccaacactacgaagaagatatacacaaggcatttgtcatccgcggcaattccgccatattgaaatgcgatataccgtcgtttgt",
               "ttcgagccgtcgtttcccaacactacgaagaagatatacacaaggcatttgtcatacgcggcaattccgccatattaaaatgtgatataccgtcgtttgt",
               "tgcgagccgtcgtttcccaacactacgaagaagatatacacaaggcatttgtcatccgcggcaattccgccatattgaaatgtgatattccgtcgtttgt",
               "tgcgagccgtcgtttcgcaacactacgaagaggatatacacaaggcttttgtcatccgcggcaattccgccatcttgaaatgtgatataccatcgttcgt",
               "tgcgcgccgtcgtttcccaacactacgaagaggatatacacaaggcttttgtcattcgcggcaattccgccatcttgaaatgtgatataccgtcatttgt",
               "tgcgagccgtcgtttcgcaacactacgaagaggatatacacaaggcgtttgtcatacgcggcaattccgctatcttgaaatgtgatataccatcgtttgt",
               "tacgagccgtggttgcgcagcattacgacaccgacgtcaataaggagtatgtcattactggcaacagtatcgtgctcaagtgtcaggtcccatctttcgt",
               "tgcgagctgtggttgcgcagaattacgacactgatgttaacaaggagtatgttataatgggaaacagcatcattctgaagtgtcaggtcccatcctttgt")

seq2<-c("gttccacccaagatgatgccgttctctttcggagagatgt",
        "gtggctcctcaagtacttcctttcgaatttggcgaaacgt",
        "gtccaaccccaaattttgcctttcgattttggagaaacgt",
        "gttcctccccaaatcgttccttttgattttggagaaatgt",
        "gttcccccacagattactccttttgacttcggtaacatgt",
        "gtaccacctcaaattgttccatttgattttggcgaaacgt",
        "gtaccaccccaaattgttccttttgaattcggcgaaacgt",
        "gtactgccgcatataactccttttgagttcgaaaacaggt",
        "gttctacccaaaattcaccccttcagtttcggtgataggt",
        "gtattgcccagaattactccattctactttgaagacatgt",
        "gtacttccacgaatatctccctttactttcgaagaaacgt",
        "gtgaaaccggacgttgtcccgtttagttttggaaagatgt",
        "gttcctccacagattttgcattttgatttcggagatatgt",
        "gtcgcaccacagataggtccattcttaatgagtgatatgt",
        "gtcacgccacaaatagctccattttcaatgagtgatatgt",
        "gtttcgccacaaataggtacgatagcgtttaccgatatgt",
        "gtcgctccgcagatattacccttgagttttggcgatatgt",
        "gtcgctccacagattcttccattcgagtttggtgaaacgt",
        "gtcgctccacagattcttccattcacttttggtgaaatgt",
        "gtcgctccacagattcttccattcacttttggtgaaatgt",
        "gtcgctccacagattcttccattcacttttggtgaaatgt",
        "gtcgctccacagattcttccattcacttttggtgaaatgt",
        "gttgcaccacagatctttcccttttcgtttggcgacacgt",
        "gttgcaccgcaaatcatacccttcgattttggtgatatgt",
        "gtggcaccacagatcgggcccttctcgttcggtgacatgt",
        "gttgcaccgcaaatcatacccttcgattttggtgacatgt",
        "gtctcaccagagatcgctcctttcgatatcgcggaaatgt",
        "gtctcaccagagatcgctcccttcgatatcgcggaaatgt",
        "gtggtaccccatataggtcctttttcgatcagtgacacgt",
        "gtcgcaccgcagatagctccgttcgtaatcaccgagatgt",
        "gtcgcaccgcagatagctccgttcgtaatcaccgagatgt",
        "gtgtctccgcagatcgcaccaatctccttcggcgacacgt",
        "gtgtctccgcagatcgcaccaatctccttcggcgacacgt",
        "gtcgccccgcaaatactgcccttcatgttcggcgacacgt",
        "gttgcgccgcatataggaccgttctcgatcagtgacacgt",
        "gttgcgccgcacatagggcctttcacgattagcgatacgt",
        "gttgccccgcaaatactgcctttcgctttcggggatacgt",
        "gtcgcaccgcaaatcgcaccgtttgaaatcgccgagacgt",
        "gttgcgccgcacatagggccttttgcgatcagtgacatgt",
        "gtcgctccgcaaatcgctccattttctatcagtgaaatgt",
        "gttgctccgcagatcttcccgttcacgttcggcgacatgt",
        "gtcgctccacagataggcccgttctcattcggcgagacgt",
        "gttccaccacaattaatacccttcgatttcggtgatatgt",
        "gttcctccgcagatcctacctttctcatttggagacatgt",
        "gttcctccgcagatcctacctttctcatttggagacatgt",
        "gttcctccaaaaacgatgcctttttcatttggggatatgt",
        "gttcctccgcagatcctacctttctcatttggagacatgt",
        "gttccacctcaaatttcaccatttacctttggagatacgt",
        "ggtcctctttaagttttattacttgcttgtaaaaatacgt",
        "gttccaccagaaattcttcctttcacttttggagagatat",
        "gtcctaccacagattctagctttctcatttggcaat--gt",
        "gttataccacagatcctgccgttctcgtttggcgctacgt",
        "gtgttaccgcagatccaagcattctcgtttggcgataagt",
        "gtgccaccccaagttttaccttttaatatcggcgatatgt",
        "gttctgccgcgcatcattcccttcgccttcgaggagacgt",
        "gttcttccgcgcatcatacccttcgccttcgaggagacgt",
        "gttcttccgcgcatcatacccttcgccttcgaggagacgt",
        "gttcttccgcgcatcatacccttcgccttcgaggagacgt",
        "gttcttccgcgcatcatacccttcgccttcgaggagacgt",
        "gttcttccgcgcatcatacccttcgccttcgaggagacgt",
        "gttctgccgcgcatcataccctttgccttcgaggagacgt",
        "gttcttccgcgcatcatacccttcgattttgaggaaacgt",
        "gttctgccgcgcatcatacccttcgccttcgaggagacgt",
        "gttctgccgcgcatcatacccttcgccttcgaggagatgt",
        "gttgtaccgcaggttacgcccttcgatttcggggaaatgt",
        "gttcctccgcaaatattgccattcgagtttggtgacacgt")

createPFM <- function(seq1) {
  library(Biostrings)
  aln <- DNAStringSet(seq1)
  pfm <- consensusMatrix(aln, as.prob=TRUE, baseOnly = TRUE)
  pfm <- pfm[1:4,]  # Select only A, C, G, T rows
  return(pfm)
}
pfm <- createPFM(seq1)
pfm2<- createPFM(seq2)

#making the sequence logo an object
seqlogo<-ggplot() + 
      geom_logo(pfm2) + 
      theme_logo() + 
      labs(title ="D") +
      scale_fill_manual(values = c("A" = "orange", 
                                   "C" = "hotpink", 
                                   "G" = "purple", 
                                   "T" = "blue"))

# making a line plot for the phyloP scores for the entire gene
setwd("C:/Users/remin/Desktop/School/Programs/BEAST.v.2.7.6/DSCAM/phyloP/nrate=4/ma")

solitary<-read.table("sol-ma.wig")
presocial<-read.table("sub-ma.wig")
eusocial<-read.table("eus-ma.wig")

df1 <- data.frame(x = 1:4779, y = solitary$V1, group = "Solitary")
df2 <- data.frame(x = 1:4779, y = presocial$V1, group = "Pre-Social")
df3 <- data.frame(x = 1:4779, y = eusocial$V1, group = "Eusocial")

combined_df <- rbind(df1, df2, df3)

line1<-ggplot(combined_df, aes(x = x, y = y, color = group)) +
  geom_line(aes(group = group), size = 0.8, color = "black", show.legend = FALSE) +  # Outline
  geom_line(aes(group = group), size = 0.5) +  # Main line
  geom_rect(aes(xmin = 1600, xmax = 1700, ymin = -2, ymax = 2),
            fill = NA, color = "red", linetype = "solid", size = 0.4, alpha=0.5) +
  labs(x = "Length of DSCAM", y = "phyloP score", title = "B") +
  scale_color_manual(values = c("#fde725","#21918c", "#440154")) +  # Specify line colors
  theme_minimal() +
  facet_wrap(~ group, ncol = 1)  # Stack graphs on top of each other

# Making a line graph for the zoomed in area 1600-1700

setwd("C:/Users/remin/Desktop/School/Programs/BEAST.v.2.7.6/DSCAM/phyloP/nrate=4/ma/1600-1700")

solzoom<-read.table("sol1620-1660.wig")
prezoom<-read.table("sub1620-1660.wig")
euszoom<-read.table("eus1620-1660.wig")

df4 <- data.frame(x = 1:40, y = solzoom$V1, group = "Solitary")
df5 <- data.frame(x = 1:40, y = prezoom$V1, group = "Pre-Social")
df6 <- data.frame(x = 1:40, y = euszoom$V1, group = "Eusocial")

combined_df2 <- rbind(df4, df5, df6)

line2<-ggplot(combined_df2, aes(x = x, y = y, color = group)) +
  geom_line(aes(group = group), size = 1, color = "black", show.legend = FALSE) +  # Outline
  geom_line(aes(group = group), size = 0.7) +  # Main line
  #geom_hline(yintercept = -2, linetype = "dashed", color = "red", size = 0.3) +  # Red dotted line at y = 0
  labs(x="",y = "phyloP score", title ="C") +
  scale_color_manual(values = c("#fde725","#21918c", "#440154")) +  # Specify line colors
  theme_minimal()+
  theme(legend.position = "none",axis.text.x = element_blank())


# making a figure for the paper
library(ggpubr)

ggarrange(line1, line2, seqlogo,ncol=1, nrow=3,heights = c(2.5, 1.3, 1))
