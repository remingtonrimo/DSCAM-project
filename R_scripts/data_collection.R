# load the neccessary libraries
library(rentrez)
library(xml2)
library(dplyr)
library(Biostrings)
library(stringr)

# start here for searching NCBi for DSCAM1 gene
search_results <- entrez_search(
  db = "gene",
  term = "DSCAM1[Gene Name] AND Insecta[Organism]",
  retmax = 1000  # Adjust retmax based on the expected number of results
)

# Fetch summaries for all gene IDs
gene_summaries <- entrez_summary(db = "gene", id = search_results$ids)

# Initialize a list to store information
gene_info_list <- list()

# Loop over summaries and extract relevant information
for (i in 1:length(gene_summaries)) {
  gene_info <- gene_summaries[[i]]
  gene_id <- gene_info$uid
  assembly_accession <- gene_info$assembly
  species_name <- gene_info$organism
  location <- gene_info$location
  
  # Store the information in a list
  gene_info_list[[i]] <- list(
    GeneID = gene_id,
    AssemblyAccession = assembly_accession,
    Species = species_name,
    Location = location
  )
}


# making the data more managable
results_matrix <- matrix(nrow = length(gene_info_list), ncol = 7)
colnames(results_matrix) <- c("Species", "CommonName", "TaxID","GeneID", "Source","Start","End")

for (i in 1:length(gene_info_list)) {
  species <- gene_info_list[[i]]$Species$scientificname
  common_name <- gene_info_list[[i]]$Species$commonname
  tax_id <- gene_info_list[[i]]$Species$taxid
  gene_id <- gene_info_list[[i]]$GeneID
  source <- gene_info_list[[i]]$Location$chraccver
  start <- gene_info_list[[i]]$Location$chrstart
  end <- gene_info_list[[i]]$Location$chrstop
# Create a matrix for the current entry
  # Check if all items are character vectors
  if (length(species) == 1 && length(common_name) == 1 && length(tax_id) == 1 && length(gene_id) == 1 && length(start) == 1 && length(end) == 1) {
    # Store data in the matrix
    results_matrix[i, ] <- c(species, common_name, tax_id, gene_id, source, start, end)
  } else {
    warning(paste("Data for index", i, "is not of length 1:", species, common_name, tax_id, gene_id, source, start, end))
  }
}

# Remove rows where GeneID (3rd column) is NA
cleaned_matrix <- results_matrix[!is.na(results_matrix[, "GeneID"]), ]

# writing a csv file for later use                                          
# write.csv(cleaned_matrix, file="dscam.csv", row.names = FALSE)

# This is for fetching the data
# Define the output FASTA file name
output_file <- "dscam1.fasta"

# Open a connection to the file
fasta_conn <- file(output_file, open = "wt")

# Loop through the cleaned_matrix to fetch sequences
for (i in 1:nrow(cleaned_matrix)) {
  message("Fetching sequence ", i, " of ", nrow(cleaned_matrix))
  tryCatch({
    # Fetch the sequence
    fasta_sequence <- entrez_fetch(
      db = "nucleotide", 
      id = cleaned_matrix[,"Source"][i], 
      rettype = "fasta",
      seq_start = cleaned_matrix[,"Start"][i], 
      seq_stop = cleaned_matrix[,"End"][i]
    )
    
    # Write the fetched sequence to the FASTA file
    writeLines(fasta_sequence, fasta_conn)
    
  }, error = function(e) {
    message("Error fetching ID ", cleaned_matrix[,"Source"][i], ": ", e$message)
  })
}

# Close the file connection
close(fasta_conn)


#### Gathering all exon positions ####
# Initialize an empty list to store gene records
gene_records <- list()

# Loop over each row of cleaned_matrix
for (i in 1:nrow(cleaned_matrix)) {
  # Extract values for the current row
  accession_id <- cleaned_matrix[i, "Source"]
  start_pos <- as.integer(cleaned_matrix[i, "Start"])
  end_pos <- as.integer(cleaned_matrix[i, "End"])
  
  # Fetch the gene record for the current accession and positions
  tryCatch({
    gene_record <- entrez_fetch(
      db = "nucleotide",
      id = accession_id,  # Current accession ID
      rettype = "gb",
      seq_start = start_pos,
      seq_stop = end_pos
    )
    
    # Store the gene record in the list, indexed by the accession ID
    gene_records[[accession_id]] <- gene_record
    
  }, error = function(e) {
    message("Error fetching gene record for accession ", accession_id, ": ", e$message)
  })
}

# Now, gene_records contains all fetched records, indexed by accession ID
# Initialize a list to store unique exon ranges for each gene record
all_unique_exons <- list()

# Iterate over each gene record
for (accession_id in names(gene_records)) {
  # Get the gene record for the current accession ID
  gene_record <- gene_records[[accession_id]]
  
  # Match all mRNA join lines in the record, accounting for complement cases
  mRNA_exons <- str_extract_all(gene_record, "(mRNA\\s+)?(complement\\s+)?join\\([^)]*\\)")
  
  # Initialize a list to store unique exon ranges for this gene record
  unique_exons_list <- list()
  
  # Check if any mRNA join lines were found
  if (length(mRNA_exons) > 0 && length(mRNA_exons[[1]]) > 0) {
    # Iterate over matched mRNA lines to capture exon ranges
    for (mRNA in mRNA_exons[[1]]) {
      # Extract only the numeric ranges within the join statement
      exons <- str_extract_all(mRNA, "\\d+\\.\\.\\d+")[[1]]  # Only captures ranges
      
      # Remove duplicates by converting to unique values
      unique_exons <- unique(exons)
      unique_exons_list <- c(unique_exons_list, list(unique_exons))
    }
    
    # Concatenate all exon ranges from this gene record into a single list
    all_exons <- unlist(unique_exons_list)
    
    # Remove duplicates across all variants
    unique_exons <- unique(all_exons)
    
    # Store the unique exons for this accession ID
    all_unique_exons[[accession_id]] <- unique_exons
  } else {
    message("No mRNA join lines found for accession ", accession_id)
  }
}

# Display or process the unique exons for all gene records
print(all_unique_exons)

# Initialize a list to store the sorted unique exons for each accession ID
sorted_unique_exons <- list()
# Loop over each accession ID in all_unique_exons
for (accession_id in names(all_unique_exons)) {
  # Get the current vector of unique exons
  unique_exons_vector <- all_unique_exons[[accession_id]]
  
  # Extract starting numbers and sort based on them
  sorted_exons <- unique_exons_vector[order(as.numeric(sub("\\..*", "", unique_exons_vector)))]
  
  # Store the sorted exons in the new list
  sorted_unique_exons[[accession_id]] <- sorted_exons
}

# Display the sorted unique exons for all gene records
print(sorted_unique_exons)

# Make an empty list
exon_matrices <- list()

# Iterate over each accession ID in sorted_unique_exons
for (accession_id in names(sorted_unique_exons)) {
  # Get the exon ranges for the current accession
  exon_ranges <- sorted_unique_exons[[accession_id]]
  
  # Split the ranges into start and stop positions
  split_ranges <- str_split(exon_ranges, "\\.\\.")
  
  # Create a matrix to hold start and stop values
  exon_matrix <- do.call(rbind, lapply(split_ranges, function(range) {
    as.numeric(range)  # Convert to numeric
  }))
  
  # Set column names
  colnames(exon_matrix) <- c("exon_start", "exon_stop")
  
  # Store the matrix in the list with the accession ID as the name
  exon_matrices[[accession_id]] <- exon_matrix
}

# Display the resulting list of matrices
print(exon_matrices)

#### Filtered exons fasta file ####
fasta_sequences <- readDNAStringSet("dscam1.fasta")

# Use sorted_unique_exons as the exon ranges
exon_ranges <- sorted_unique_exons

# Output file for filtered sequences
output_file <- "filtered_exons.fasta"
fasta_conn <- file(output_file, open = "wt")

# Check if fasta_sequences are loaded correctly
if (length(fasta_sequences) == 0) {
  stop("No sequences found in the FASTA file.")
}

# Loop through each sequence in the FASTA file
for (accession in names(fasta_sequences)) {
  # Get the sequence for the current accession
  sequence <- fasta_sequences[[accession]]
  
  # Extract base accession number (the part before the first ':')
  base_accession <- unlist(strsplit(accession, ":"))[1]
  
  # Find the species name from cleaned_matrix
  species_name <- cleaned_matrix[cleaned_matrix[,"Source"] == base_accession, "Species"]
  
  # Debugging output to check the current accession being processed
  message("Processing accession: ", base_accession)
  
  # Check if there are exon ranges for the current base accession
  if (base_accession %in% names(exon_ranges)) {
    ranges <- exon_ranges[[base_accession]]
    
    # Initialize a character vector to hold the exon sequences
    exon_sequence_parts <- character()
    
    # Loop through each exon range and extract sequences
    for (range in ranges) {
      # Split the range into start and end positions
      positions <- as.numeric(unlist(strsplit(range, "\\.\\.")))
      start <- positions[1]
      end <- positions[2]
      
      # Ensure start and end positions are valid
      if (start > 0 && end <= nchar(as.character(sequence))) {
        # Extract the nucleotides corresponding to the exon range
        exon_sequence <- subseq(sequence, start = start, end = end)
        exon_sequence_parts <- c(exon_sequence_parts, as.character(exon_sequence))
      } else {
        message("Invalid range for accession ", base_accession, ": ", range)
      }
    }
    
    # Combine all parts into a single sequence
    if (length(exon_sequence_parts) > 0) {
      combined_exon_sequence <- paste(exon_sequence_parts, collapse = "")
      
      # Write to the new FASTA file with species name in the header
      writeLines(paste(">", base_accession, species_name, sep = " "), fasta_conn)
      writeLines(combined_exon_sequence, fasta_conn)
    } else {
      message("No valid exon sequences found for accession: ", base_accession)
    }
  } else {
    message("No exon ranges found for accession: ", base_accession)
  }
}

# Close the file connection
close(fasta_conn)

# Final message to indicate completion
message("Filtered exons written to ", output_file)














