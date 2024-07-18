if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("seqLogo")

library(seqLogo)
# Function to create a PFM from a list of sequences
createPFM <- function(sequences) {
  library(Biostrings)
  aln <- DNAStringSet(sequences)
  pfm <- consensusMatrix(aln, as.prob=TRUE, baseOnly = TRUE)
  pfm <- pfm[1:4,]  # Select only A, C, G, T rows
  return(pfm)
}

# Example list of DNA sequences
sequences <- c("tgcgcgcagtggtgtcccagcactacatcacggaggccgaaaacgagtacgttatccggggcaactcggcggtcatgaagtgcaagatccccagcttcgt",
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
# Create PFM from sequences
pfm <- createPFM(sequences)
pfm2<- createPFM(seq2)
# Convert the PFM to a PWM
pwm <- makePWM(pfm)
pwm2<- makePWM(pfm2)
#?makePWM
# Plot the sequence logo
seqLogo(pwm)
seqLogo(pwm2, fill=c(A='slateblue', C='orchid', G='green', T='blue'))
?seqLogo


par(mfrow = c(2, 2))
plot(iris$Sepal.Length~iris$Petal.Length)
plot(mtcars$mpg~mtcars$disp)
