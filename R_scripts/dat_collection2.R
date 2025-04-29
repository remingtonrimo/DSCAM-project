#### Gathering all of the transcript variants for each DSCAM gene, aligning 
#### them with mafft and then cutting out introns with BMGE
# Initialize a list to store the transcript variants for each species
transcript_variants_list <- list()

# Loop over each row in cleaned_matrix
for (i in 1:nrow(cleaned_matrix)) {
  
  # Get the species name and GeneID
  species_name <- cleaned_matrix[i, "Species"]
  gene_id <- cleaned_matrix[i, "GeneID"]
  
  # Fetch linked nucleotide records associated with the GeneID
  tryCatch({
    links <- entrez_link(dbfrom = "gene", id = gene_id, db = "nucleotide")
    
    # Check if nucleotide records are found
    if (!is.null(links$links$gene_nuccore_refseqrna)) {
      # Retrieve summaries for all linked nucleotide records
      summaries <- entrez_summary(db = "nucleotide", id = links$links$gene_nuccore_refseqrna)
      
      # Store each transcript variant summary under the species name
      transcript_variants_list[[species_name]] <- summaries
    } else {
      message("No nucleotide links found for GeneID ", gene_id, " in ", species_name)
    }
    
  }, error = function(e) {
    message("Error fetching links for GeneID ", gene_id, ": ", e$message)
  })
}

# Display the collected transcript variant summaries
print(transcript_variants_list)


# Loop over each species and its corresponding transcript variants
for (species_name in names(transcript_variants_list)) {
  message("Processing species: ", species_name)
  
  # Get the summaries for the transcript variants
  transcript_variants <- transcript_variants_list[[species_name]]
  
  if (length(transcript_variants) == 0) {
    message("No transcript variants found for species: ", species_name)
    next  # Skip to the next species if no variants are found
  }
  
  # Create a connection to a FASTA file for this species
  fasta_file_path <- file.path("temp_fasta_files", paste0(species_name, "_transcripts.fasta"))
  
  # Use tryCatch to handle errors with file connection
  fasta_conn <- tryCatch({
    file(fasta_file_path, open = "wt")
  }, error = function(e) {
    message("Error opening file for species ", species_name, ": ", e$message)
    return(NULL)  # Return NULL to indicate failure to open file
  })
  
  # If the file connection failed, skip this species
  if (is.null(fasta_conn)) {
    next
  }
  
  # Fetch sequences for each transcript variant
  for (variant in transcript_variants) {
    # Check if variant is a list and has 'uid'
    if (is.list(variant) && !is.null(variant$uid)) {
      variant_id <- variant$uid  # Assuming 'uid' holds the sequence ID
    } else {
      message("Variant UID is missing or invalid for species: ", species_name)
      next  # Skip this variant if UID is missing or invalid
    }
    
    # Fetch the transcript sequence
    tryCatch({
      fasta_sequence <- entrez_fetch(
        db = "nucleotide", 
        id = variant_id, 
        rettype = "fasta"
      )
      
      # Write the fetched sequence to the FASTA file
      writeLines(fasta_sequence, fasta_conn)
      
    }, error = function(e) {
      message("Error fetching variant ID ", variant_id, " for species ", species_name, ": ", e$message)
    })
  }
  
  # Close the file connection
  close(fasta_conn)
}

# Notify completion
message("All sequences have been fetched and written to FASTA files.")












# Create a directory to store temporary FASTA files if it doesn't exist
dir.create("temp_fasta_files", showWarnings = FALSE)

# Loop over each species and its corresponding transcript variants
for (species_name in names(transcript_variants_list)) {
  message("Processing species: ", species_name)
  
  # Get the summaries for the transcript variants
  transcript_variants <- transcript_variants_list[[species_name]]
  
  if (length(transcript_variants) == 0) {
    message("No transcript variants found for species: ", species_name)
    next  # Skip to the next species if no variants are found
  }
  
  # Create a connection to a FASTA file for this species
  fasta_file_path <- file.path("temp_fasta_files", paste0(species_name, "_transcripts.fasta"))
  
  # Use tryCatch to handle errors with file connection
  fasta_conn <- tryCatch({
    file(fasta_file_path, open = "wt")
  }, error = function(e) {
    message("Error opening file for species ", species_name, ": ", e$message)
    return(NULL)  # Return NULL to indicate failure to open file
  })
  
  # If the file connection failed, skip this species
  if (is.null(fasta_conn)) {
    next
  }
  
  # Fetch sequences for each transcript variant
  for (variant in transcript_variants) {
    # Check if variant is a list and has 'uid'
    if (is.list(variant) && !is.null(variant$uid)) {
      variant_id <- variant$uid  # Assuming 'uid' holds the sequence ID
    } else {
      message("Variant UID is missing or invalid for species: ", species_name)
      next  # Skip this variant if UID is missing or invalid
    }
    
    # Fetch the transcript sequence
    tryCatch({
      fasta_sequence <- entrez_fetch(
        db = "nucleotide", 
        id = variant_id, 
        rettype = "fasta"
      )
      
      # Write the fetched sequence to the FASTA file
      writeLines(fasta_sequence, fasta_conn)
      
    }, error = function(e) {
      message("Error fetching variant ID ", variant_id, " for species ", species_name, ": ", e$message)
    })
  }
  
  # Fetch the entire gene sequence using Source, Start, and End information
  gene_info <- cleaned_matrix[which(cleaned_matrix[, "Species"] == species_name), ]
  
  # Fetch the entire gene sequence if needed
  tryCatch({
    # Access 'Source', 'Start', and 'End' directly from gene_info
    gene_sequence <- entrez_fetch(
      db = "nucleotide", 
      id = gene_info["Source"], 
      rettype = "fasta",
      seq_start = as.integer(gene_info["Start"]), 
      seq_stop = as.integer(gene_info["End"])
    )
    
    # Write the entire gene sequence to the FASTA file
    writeLines(gene_sequence, fasta_conn)
    
  }, error = function(e) {
    message("Error fetching gene sequence for species ", species_name, ": ", e$message)
  })
  
  # Close the file connection
  close(fasta_conn)
}

# Notify completion
message("All sequences have been fetched and written to FASTA files.")

