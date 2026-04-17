# ====== hyperthermophiles =======

library(G4SNVHunter)
library(Biostrings)
library(tidyverse)

setwd('/home/dimitri/project4/ERC/hyperthermophiles')
root_dir <- "."

is_chromosome <- function(seqnames) {
  is_plasmid <- grepl("plasmid|p[0-9]+", seqnames, ignore.case = FALSE) |
    grepl("plasmid", seqnames, ignore.case = TRUE)
  return(!is_plasmid)
}

run_g4hunt_chr <- function(fasta_path, out_file, window=25, threshold=1.2) {
  all_seqs <- Biostrings::readDNAStringSet(fasta_path)
  seq_names <- names(all_seqs)
  chr_idx <- which(is_chromosome(seq_names))
  if (length(chr_idx) == 0) { warning("No chromosomes in: ", fasta_path); return(NULL) }
  chr_seqs <- all_seqs[chr_idx]
  clean_strings <- gsub("[^ACGTUN]", "N", toupper(as.character(chr_seqs)))
  clean_chr_seqs <- Biostrings::DNAStringSet(clean_strings)
  names(clean_chr_seqs) <- seq_names[chr_idx]
  message("  Analyzing ", length(clean_chr_seqs), " chromosome(s)...")
  tryCatch({
    g4res <- G4HunterDetect(clean_chr_seqs, window_size = window, threshold = threshold)
    exportG4(g4res, filename = out_file, include_metadata = TRUE)
  }, error = function(e) message("  Error: ", e$message))
}

# --- Run G4Hunter ---
dirs <- list.dirs(root_dir, recursive = FALSE)
for (d in dirs) {
  fasta <- list.files(d, pattern="\\.fna$", full.names = TRUE)
  if (length(fasta) == 0) next
  base <- basename(d)
  outfile <- file.path(root_dir, paste0(base, "_chr_G4Hunter.csv"))
  message("Processing ", base)
  run_g4hunt_chr(fasta_path = fasta, out_file = outfile)
}

# --- Build the table ---
name_map <- c(
  "GCF_000223395.1" = "Pyrolobus fumarii",
  "GCF_000007185.1" = "Methanopyrus kandleri",
  "GCF_000008625.1" = "Aquifex aeolicus"
)

g4_files <- list.files(pattern="_chr_G4Hunter.csv")


# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# --------- USE reformatting_csv.sh !!!!!!!!!!!!!!!
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

results_table <- map_df(g4_files, function(file) {
  df <- read_csv(file)
  # df <- read_csv(file, comment = "#")
  genome_name <- str_remove(file, "_chr_G4Hunter.csv")
  fasta_file <- list.files(genome_name, pattern="\\.fna$", full.names = TRUE)
  seqs <- readDNAStringSet(fasta_file)
  chr_seqs <- seqs[which(is_chromosome(names(seqs)))]
  genome_len <- sum(width(chr_seqs))
  gc_pct <- sum(letterFrequency(chr_seqs, "GC")) / genome_len * 100
  
  tibble(
    Name         = name_map[genome_name],
    Genome_size  = genome_len,
    GC_percent   = round(gc_pct, 2),
    G4_count     = nrow(df),
    G4_density   = round(nrow(df) / (genome_len / 1e6), 2)
  )
})

print(results_table)
write_csv(results_table, "hyperthermophiles_G4_summary.csv")


# ====== ASGARD =======

# first create a list of accessions and download all the necessary genomes for the analysis

library(G4SNVHunter)
library(Biostrings)
library(tidyverse)

setwd('/home/dimitri/project4/ERC/asgard')
getwd()
root_dir <- "."

is_chromosome <- function(seqnames) {
  is_plasmid <- grepl("plasmid|p[0-9]+", seqnames, ignore.case = FALSE) |
    grepl("plasmid", seqnames, ignore.case = TRUE)
  return(!is_plasmid)
}

run_g4hunt_chr <- function(fasta_path, out_file, window=25, threshold=1.2) {
  all_seqs <- Biostrings::readDNAStringSet(fasta_path)
  seq_names <- names(all_seqs)
  chr_idx <- which(is_chromosome(seq_names))
  if (length(chr_idx) == 0) { warning("No chromosomes in: ", fasta_path); return(NULL) }
  chr_seqs <- all_seqs[chr_idx]
  clean_strings <- gsub("[^ACGTUN]", "N", toupper(as.character(chr_seqs)))
  clean_chr_seqs <- Biostrings::DNAStringSet(clean_strings)
  names(clean_chr_seqs) <- seq_names[chr_idx]
  message("  Analyzing ", length(clean_chr_seqs), " chromosome(s)...")
  tryCatch({
    g4res <- G4HunterDetect(clean_chr_seqs, window_size = window, threshold = threshold)
    exportG4(g4res, filename = out_file, include_metadata = TRUE)
  }, error = function(e) message("  Error: ", e$message))
}

# --- Run G4Hunter ---
dirs <- list.dirs(root_dir, recursive = FALSE)
for (d in dirs) {
  fasta <- list.files(d, pattern="\\.fna$", full.names = TRUE)
  if (length(fasta) == 0) next
  base <- basename(d)
  outfile <- file.path(root_dir, paste0(base, "_chr_G4Hunter.csv"))
  message("Processing ", base)
  run_g4hunt_chr(fasta_path = fasta, out_file = outfile)
}

# --- Build the table ---
name_map <- c(
  "GCF_000223395.1" = "Pyrolobus fumarii",
  "GCF_000007185.1" = "Methanopyrus kandleri",
  "GCF_000008625.1" = "Aquifex aeolicus"
)

g4_files <- list.files(pattern="_chr_G4Hunter.csv")


# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# --------- USE reformatting_csv.sh !!!!!!!!!!!!!!!
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

results_table <- map_df(g4_files, function(file) {
  df <- read_csv(file)
  # df <- read_csv(file, comment = "#")
  genome_name <- str_remove(file, "_chr_G4Hunter.csv")
  fasta_file <- list.files(genome_name, pattern="\\.fna$", full.names = TRUE)
  seqs <- readDNAStringSet(fasta_file)
  chr_seqs <- seqs[which(is_chromosome(names(seqs)))]
  genome_len <- sum(width(chr_seqs))
  gc_pct <- sum(letterFrequency(chr_seqs, "GC")) / genome_len * 100
  
  tibble(
    Name         = name_map[genome_name],
    Genome_size  = genome_len,
    GC_percent   = round(gc_pct, 2),
    G4_count     = nrow(df),
    G4_density   = round(nrow(df) / (genome_len / 1e6), 2)
  )
})

print(results_table)
write_csv(results_table, "hyperthermophiles_G4_summary.csv")



