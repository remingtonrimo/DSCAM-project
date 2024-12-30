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




gene_info <- entrez_summary(db = "gene", id = 1271065)
exons <- gene_info$exons  # Access the exon information












####Weird code, I'm not exactly sure what this is actually fetching####

# Load your data (example for species names and positions)
species_names <- cleaned_matrix[, "Species"]  # Assuming cleaned_matrix is a data frame
all_accessions <- names(sorted_unique_exons)  # Assuming sorted_unique_exons is a list with positions per accession

# Initialize FASTA output file
output_fasta <- file("all_concatenated_sequences2.fasta", "w")

# Define a function to extract and concatenate specific positions from a sequence
get_subsequences <- function(sequence, positions) {
  concatenated_sequence <- ""
  for (pos in positions) {
    start <- as.numeric(strsplit(pos, "\\.\\.")[[1]][1])
    end <- as.numeric(strsplit(pos, "\\.\\.")[[1]][2])
    concatenated_sequence <- paste0(concatenated_sequence, substr(sequence, start, end))
  }
  return(concatenated_sequence)
}

# Function to fetch sequence data with retry logic
fetch_sequence_with_retry <- function(accession, retries = 3, delay = 5) {
  for (attempt in 1:retries) {
    result <- tryCatch({
      entrez_fetch(db = "nuccore", id = accession, rettype = "fasta", retmode = "text")
    }, error = function(e) {
      message("Attempt ", attempt, " failed for ", accession, ": ", e$message)
      Sys.sleep(delay)
      NULL  # Return NULL if there's an error
    })
    if (!is.null(result)) return(result)  # If successful, return the result
  }
  stop("Failed to fetch sequence after ", retries, " attempts for ", accession)
}

# Open output FASTA file
output_fasta <- file("all_concatenated_sequences2.fasta", "w")

# Iterate over each accession and retrieve concatenated sequences
for (i in seq_along(all_accessions)) {
  accession <- all_accessions[i]
  species_name <- species_names[i]
  positions <- sorted_unique_exons[[accession]]
  
  # Fetch the entire sequence with retry logic
  sequence_data <- fetch_sequence_with_retry(accession)
  sequence_lines <- strsplit(sequence_data, "\n")[[1]]
  sequence <- paste(sequence_lines[-1], collapse = "")
  
  # Extract and concatenate all specified positions
  concatenated_sequence <- get_subsequences(sequence, positions)
  
  # Write to FASTA file with accession and species name in the header
  writeLines(paste0(">", accession, " ", species_name, "\n", concatenated_sequence), output_fasta)
}

# Close the FASTA file after writing all records
close(output_fasta)






####new code####
# first I need to make sorted_unique_exons a list of matrices

class(exon_start)

# Initialize a list to store sequences
all_sequences <- list()

# Loop through each accession in cleaned_matrix
for (i in 1:nrow(cleaned_matrix)) {
  gene_start <- cleaned_matrix[,"Start"][i]
  sequence <- cleaned_matrix[,"Source"][i]
  gene_start<- as.numeric(gene_start)
  # Access the exon matrix for the current sequence
  exon_matrix_sequence <- exon_matrices[[sequence]]
  
  # Get exon start and stop positions
  exon_start <- exon_matrix_sequence[,"exon_start"]
  exon_stop <- exon_matrix_sequence[,"exon_stop"]
  
  # Loop through each exon to fetch sequences
  for (j in 1:length(exon_start)) {
    start <- gene_start + exon_start[j]
    end <- gene_start + exon_stop[j]
    
    # Fetch the sequence from NCBI
    dscam <- entrez_fetch(
      db = "nucleotide", 
      id = sequence, 
      rettype = "fasta",
      seq_start = start, 
      seq_stop = end
    )
    
    # Concatenate sequences based on accession and species name
    seq_name <- paste(sequence, cleaned_matrix[,"Species"][i], sep = "_")
    
    # Add the fetched sequence to the list
    if (!is.null(all_sequences[[seq_name]])) {
      all_sequences[[seq_name]] <- paste0(all_sequences[[seq_name]], "\n", dscam)
    } else {
      all_sequences[[seq_name]] <- dscam
    }
  }
}

# Print or process the concatenated sequences
# Specify the output file path
output_file <- "all_sequences.fasta"

# Open a connection to write to the file
file_conn <- file(output_file, open = "w")

# Iterate over each sequence in all_sequences
for (seq_name in names(all_sequences)) {
  # Write the header with the sequence name
  writeLines(paste0(">", seq_name), file_conn)
  
  # Write the sequence itself
  writeLines(all_sequences[[seq_name]], file_conn)
}

# Close the file connection
close(file_conn)

# Confirmation message
cat("FASTA file created:", output_file, "\n")

