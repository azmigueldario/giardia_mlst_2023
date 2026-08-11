#!/usr/bin/env Rscript

# Facilitate running form the CLI
suppressPackageStartupMessages(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop("Usage: parse_annotations.R <interpro_tsv_file> <diamond_tsv_outfile> <script_output_file>")
}

interpro_file <- args[1]
diamond_file  <- args[2]
output_file   <- args[3]

# -------------------------------------------------------------------------------------
# Define column names, create function to merge annotations
# -------------------------------------------------------------------------------------

interpro_cols <- c("locus_id", "md5", "seq_len", 
                   "analysis_source", "signature_accession", "signature_description",
                   "start", "stop", "e_value",
                   "status", "date", "interpro_accession",
                   "interpro_description")

# Preserve the source of every annotation and merge annotation results
#   1. Verify that columns are not empty
#   2. Extract info about annotation source and result
#   3. Merge as '[source1] annotation1, [source2] annotation2...'
parse_interpro_results <-  function(sources, values){

  valid  <- !is.na(values) & values != "" & values != "-" 
  if (!any(valid)) return(NA_character_)

  source <- sources[valid]
  value <- values[valid]

  output  <- 
    paste0("[", source, "] ", value) |>
    unique() |>
    paste(collapse = "; ")
  
  return(output)
}

# -------------------------------------------------------------------------------------
# Load Interpro annotation results and collapse annotations from the same type
# -------------------------------------------------------------------------------------

interpro_data <- 
  read_tsv(interpro_file, 
           col_names = interpro_cols, 
           na = c("", "-", "NA"), show_col_types = FALSE) |>
  filter(!is.na(interpro_accession)) |>
  group_by(locus_id) |>
  summarise(
    interpro_domains        = parse_interpro_results(analysis_source, interpro_description),
    interpro_accession_ids  = parse_interpro_results(analysis_source, interpro_accession),
    .groups = "drop"
  )

# -------------------------------------------------------------------------------------
# Load DIAMOND - GiardiaDB results and parse them
# -------------------------------------------------------------------------------------

diamond_cols <- c("locus_id", "giardia_db_id", "percent_identity", 
                  "align_length", "evalue", "protein_title")


# Function to parse and clean annotations from giardiaDB
parse_giardiadb  <-  function(column, regex, 
                              prefix = sub("=.*", "=", regex)){
   value <- str_extract(column, regex)
   clean_value  <- str_remove(value, prefix)
   return(clean_value)
}

diamond_data <- 
  read_tsv(diamond_file, 
           col_names = diamond_cols, 
           show_col_types = FALSE) |>
  arrange(locus_id, evalue) |>
  distinct(locus_id, .keep_all = TRUE) |>
  mutate(
    # Extract gene_product, protein title, gene_id
    giardiadb_gene_product  = parse_giardiadb(protein_title, "gene_product=([^|]+)"),
    giardiadb_protein_title = if_else(is.na(giardiadb_gene_product), str_trim(protein_title), str_trim(giardiadb_gene_product)),
    giardiadb_gene_id       = parse_giardiadb(protein_title, "gene=([^|]+)"),
    giardiadb_organism      = parse_giardiadb(protein_title, "organism=([^|])+"),
    giardiadb_proteinlength = parse_giardiadb(protein_title, "protein_length=([^|]+)")
  )

# -------------------------------------------------------------------------------------
# Merge parsed results
# -------------------------------------------------------------------------------------

final_annotations <- 
  full_join(diamond_data, interpro_data, by = "locus_id") |>
  select(locus_id, giardia_db_id, interpro_domains, 
         percent_identity, evalue, align_length,
         interpro_accession_ids,  giardiadb_gene_id,
         giardiadb_gene_product, giardiadb_organism, giardiadb_proteinlength)

write_tsv(final_annotations, output_file)
cat("Merge complete! Total unique loci annotated:", nrow(final_annotations), "\n")