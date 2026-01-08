# install.packages(c(
#   "rio",
#   "writexl"
# ))

setwd('/home/dimitri/project4')
getwd()

library(tidyr)
library(dplyr)
library(ggplot2)
library(readxl)
library(rio)
library(openxlsx)
library(writexl)
library(gplots)





thermobase <- read.csv(
  # "D:/project4_X_Lionel/tb.csv", 
  "tb.csv", 
  header = TRUE,
  sep = ",",
  fileEncoding = "latin1" 
)

# cleaning the database
thermobase <- thermobase[, colSums(!is.na(thermobase)) > 0]

# Subsetting the database
library(dplyr)
library(stringr)

names_to_match <- c(
  "Haloferax volcanii",
  "Pyrococcus abyssi",
  "Pyrococcus yayanosii",
  "Saccharolobus solfataricus",
  "Thermococcus barophilus"
)

thermobase_sub <- thermobase %>%
  filter(
    str_detect(Name, paste(names_to_match, collapse = "|"))
  )

thermobase_sub


colnames(thermobase_sub)

# -------------------------------  Running G4hunter 


library(G4SNVHunter)
library(Biostrings)
library(GenomicRanges)

#!/usr/bin/env Rscript

setwd('/home/dimitri/project4/genomes/')
getwd()


# directory with genome subdirs
root_dir <- "."

# list genome directories
dirs <- list.dirs(root_dir, recursive = FALSE)

# function to run detection and save results
run_g4hunt <- function(fasta_path, out_file, window=25, threshold=1.2) {
  seqs <- loadSequence(seq_path = fasta_path)
  # detect G4 sequences
  g4res <- G4HunterDetect(seqs,
                          window_size = window,
                          threshold = threshold)
  # export as simple text with seqnames, start, end, score
  exportG4(g4res, filename = out_file, include_metadata = TRUE)
}

for (d in dirs) {
  fasta <- list.files(d, pattern="\\.fna$", full.names = TRUE)
  if (length(fasta) == 0) next
  
  base <- basename(d)
  outfile <- file.path(root_dir, paste0(base, "_G4Hunter.csv"))
  
  
  
  message("Processing ", base)
  run_g4hunt(fasta_path = fasta, out_file = outfile)
}


################# ONLY CHROMOSOMES ###################################



#!/usr/bin/env Rscript

library(G4SNVHunter)
library(Biostrings)
library(GenomicRanges)

# set working directory to the parent containing genome folders
setwd('/home/dimitri/project4/genomes/')

root_dir <- "."

dirs <- list.dirs(root_dir, recursive = FALSE)

# helper to filter out plasmids
is_chromosome <- function(seqnames) {
  # common patterns for plasmid names: "plasmid", "p", "circle", etc.
  # adjust as needed for your datasets
  !grepl("plasmid|p[0-9]+$|p[0-9]+\\.", seqnames, ignore.case = TRUE)
}

# run G4Hunter on chromosome sequences only
run_g4hunt_chr <- function(fasta_path, out_file, window=25, threshold=1.2) {
  
  # read all sequences
  all_seqs <- readDNAStringSet(fasta_path)
  seq_names <- names(all_seqs)
  
  # pick chromosome sequences only
  chr_idx <- which(is_chromosome(seq_names))
  if (length(chr_idx) == 0) {
    warning("No chromosome sequences detected in ", fasta_path)
    return(NULL)
  }
  
  chr_seqs <- all_seqs[chr_idx]
  
  # detect G4s
  g4res <- G4HunterDetect(chr_seqs,
                          window_size = window,
                          threshold = threshold)
  
  # export
  exportG4(g4res, filename = out_file, include_metadata = TRUE)
}

for (d in dirs) {
  # find fasta
  fasta <- list.files(d, pattern="\\.fna$", full.names = TRUE)
  if (length(fasta) == 0) next
  
  base <- basename(d)
  outfile <- file.path(root_dir, paste0(base, "_chr_G4Hunter.csv"))
  
  message("Processing ", base)
  run_g4hunt_chr(fasta_path = fasta, out_file = outfile)
}


# ---- REFORMAT output .CSV files using the reformatting_csv.sh script: ----------


# -----------------------------------------------------------

setwd('/home/dimitri/project4/genomes/')
getwd()

haloferax <- read.csv("haloferax_volcanii_chr_G4Hunter.csv", header = TRUE, stringsAsFactors = FALSE)
barophilus <- read.csv("thermococcus_barophilus_chr_G4Hunter.csv", header = TRUE, stringsAsFactors = FALSE)
yayanosii <- read.csv("pyrococcus_yayanosii_chr_G4Hunter.csv", header = TRUE, stringsAsFactors = FALSE)
abyssi <- read.csv("pyrococcus_abyssi_chr_G4Hunter.csv", header = TRUE, stringsAsFactors = FALSE)


colnames(haloferax)

# ----------- CLUSTERING ---------------


thermobase_sub$salinity_num <-
  as.numeric(thermobase_sub$Avg..Opt..salinity....)

thermobase_sub <- thermobase_sub[!is.na(thermobase_sub$salinity_num), ]


install.packages("ggdendro")

library(tidyverse)
library(Biostrings)
library(ggdendro)

# set working directory
setwd("/home/dimitri/project4/genomes/")

# list all chr G4Hunter files
g4_files <- list.files(pattern="_chr_G4Hunter.csv")

# compute features per genome
g4_features <- map_df(g4_files, function(file) {
  df <- read_csv(file)
  
  genome_name <- str_remove(file, "_chr_G4Hunter.csv")
  
  # load corresponding genome sequence to get length and GC%
  fasta_file <- list.files(genome_name, pattern="\\.fna$", full.names = TRUE)
  seqs <- readDNAStringSet(fasta_file)
  chrom <- seqs[which.max(width(seqs))]
  genome_len <- sum(width(chrom))
  gc_pct <- sum(letterFrequency(chrom, "GC")) / genome_len * 100
  
  tibble(
    genome      = genome_name,
    g4_count    = nrow(df),
    g4_density  = nrow(df) / (genome_len / 1e6),  # per Mb
    mean_score  = mean(df$score, na.rm = TRUE),
    GC_percent  = gc_pct
  )
})

mean(haloferax$score)

# load your ecological metadata
# thermobase <- read_csv("thermobase_sub.csv")
thermobase_sub

thermobase_sub$Name
g4_features$genome
g4_features$g4_count


name_map <- c(
  "Pyrococcus yayanosii CH1"              = "pyrococcus_yayanosii",
  "Pyrococcus abyssi GE5"                 = "pyrococcus_abyssi",
  "Thermococcus barophilus MP, DSM 11836" = "thermococcus_barophilus",
  "Haloferax volcanii"                    = "haloferax_volcanii",
  "Saccharolobus solfataricus"            = "saccharolobus_solfataricus"
)

thermobase_sub$genome <- name_map[thermobase_sub$Name]


# join by genome name
dat <- g4_features %>%
  left_join(thermobase_sub, by = c("genome"))

print(dat)

# select variables to cluster on (standardize them)
clust_data <- dat %>%
  select(g4_density, mean_score, GC_percent, salinity_num) %>%
  scale()

# hierarchical clustering
hc <- hclust(dist(clust_data), method = "ward.D2")

# plot dendrogram
plot(hc, labels = dat$genome, main = "Cluster Dendrogram (G4 + Salinity)", xlab="", sub="")

# ---- ANOTHER PLOT - G4 Density vs Salinity ----

library(ggplot2)

ggplot(dat, aes(x = salinity_num, y = g4_density, label = genome)) +
  geom_point(size = 3) +
  geom_text(nudge_y = 0.1) +
  labs(x = "Avg Optimal Salinity (%)", y = "G4 Density (per Mb)",
       title = "G4 Density vs Salinity") +
  theme_minimal()



# --- consensus pattern for HALOFERAX ----

install.packages("ggseqlogo")
library(ggseqlogo)

g4seqs <- haloferax$sequence
g4seqs <- g4seqs[nchar(g4seqs) >= 15]
ggseqlogo(g4seqs)

table(nchar(g4seqs))
g4seqs25 <- g4seqs[nchar(g4seqs) == 25]
ggseqlogo(g4seqs25)

# --- consensus pattern for BAROPHILUS ----

g4seqs <- barophilus$sequence
g4seqs <- g4seqs[nchar(g4seqs) >= 15]

table(nchar(g4seqs))
g4seqs25 <- g4seqs[nchar(g4seqs) == 25]
ggseqlogo(g4seqs25)

# --- consensus pattern for YAYANOSII ----

g4seqs <- yayanosii$sequence
g4seqs <- g4seqs[nchar(g4seqs) >= 15]

table(nchar(g4seqs))
g4seqs25 <- g4seqs[nchar(g4seqs) == 25]
ggseqlogo(g4seqs25)

# --- consensus pattern for ABYSSI ----

g4seqs <- abyssi$sequence
g4seqs <- g4seqs[nchar(g4seqs) >= 15]

table(nchar(g4seqs))
g4seqs25 <- g4seqs[nchar(g4seqs) == 25]
ggseqlogo(g4seqs25)



thermobase_sub$Min..pH



library(dplyr)

pca <- prcomp(dat %>% select(g4_density, mean_score, GC_percent), scale=TRUE)

install.packages("ggfortify")

library(ggfortify)
autoplot(pca, data=dat, colour="salinity_num")


dat$genome



# # --- consensus pattern for HALOFERAX 23, 24, 26----
# 
# install.packages("ggseqlogo")
# library(ggseqlogo)
# 
# g4seqs <- haloferax$sequence
# g4seqs <- g4seqs[nchar(g4seqs) >= 23 & nchar(g4seqs) <= 26]
# ggseqlogo(g4seqs)
# 
# table(nchar(g4seqs))
# g4seqs25 <- g4seqs[nchar(g4seqs) == 23]
# ggseqlogo(g4seqs25)
# 
# g4seqs25 <- g4seqs[nchar(g4seqs) == 24]
# ggseqlogo(g4seqs25)
# 
# g4seqs25 <- g4seqs[nchar(g4seqs) == 25]
# ggseqlogo(g4seqs25)



