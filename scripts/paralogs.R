library(tidyverse)
library(argparse)
library(tools)

# Define the command-line arguments
#parser <- ArgumentParser(description = "Paralogs sum")
#parser$add_argument("--cov_pcg", help = "Path to cov_pcg.tsv file")
#parser$add_argument("--blast_matrix", help = "Path similarity matrix")
#parser$add_argument("--output", help = "Path to output folder")

# Parse the command-line arguments
#args <- parser$parse_args()

# Read the input files
#pcg_cov <- read_tsv(args$cov_pcg, col_names = c("gene", "gc"), show_col_types = FALSE)
#matrix <- read_tsv(args$blast_matrix, col_names = c("gene", "similar_to", "pid", "id_bases", "gene_length"), show_col_types = FALSE)

pcg_cov <- read_tsv("", col_names = c("position", "gene", "raw_coverage", "copynumber", "gc-copynumber"), show_col_types = FALSE)
matrix <- read_tsv("", col_names = c("gene", "similar_to", "pid", "id_bases", "gene_length"), show_col_types = FALSE)

# Analyze
with_paralogs <- filter(pcg_cov, gene %in% matrix$gene)

summed_paralogs <- tibble(position = character(),
                          gene = character(),
                          raw_coverage = double(),
                          copynumber = double(),
                          gc_copynumber = double(),
                          expected_copynumber = n())

for (query in matrix$gene){
  paralogs <- filter(matrix, gene==query)
  coverage_paralogs <- filter(pcg_cov, gene==query | gene %in% paralogs$similar_to) %>%
    summarise(position = paste0(position, collapse = "-"),
              gene = paste0(gene, collapse = "-"),
              raw_coverage = sum(raw_coverage),
              copynumber = sum(copynumber),
              gc_copynumber = sum(gc_copynumber),
              expected_copynumber = n())
  summed_paralogs <- bind_rows(summed_paralogs, coverage_paralogs)
}

without_paralogs <- pcg_cov %>% filter(!(gene %in% matrix$gene)) %>% mutate(expected_copynumer=1)

final_file <- bind_rows(summed_paralogs, without_paralogs)

write_tsv(final_file, "")



# Write the output files
#basename <- file_path_sans_ext(basename(args$pcg_cn))

#write_tsv(normalized_pcg, paste0(args$output, "/", basename, "-gcbias.tsv"))
#ggsave(correlation_plot, file=paste0(args$output, "/", basename, "-gcbias.png"), dpi=300)