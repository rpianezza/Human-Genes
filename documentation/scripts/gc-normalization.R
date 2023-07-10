library(tidyverse)
library(argparse)
library(tools)

# Define the command-line arguments
parser <- ArgumentParser(description = "GC-bias normalization")
parser$add_argument("--gc_scg", help = "Path to gc_scg file")
parser$add_argument("--cov_scg", help = "Path to cov_scg file")
parser$add_argument("--pcg_cn", help = "Path to pcg_cn file (normalized)")
parser$add_argument("--gc_pcg", help = "Path to gc_pcg file")
parser$add_argument("--output", help = "Path to output folder")

# Parse the command-line arguments
args <- parser$parse_args()

# Read the input files
gc_scg <- read_tsv(args$gc_scg, col_names = c("gene", "gc"), show_col_types = FALSE)
cov_scg <- read_tsv(args$cov_scg, col_names = c("gene", "cov"), show_col_types = FALSE) %>%
  select(gene, cov) %>% filter(!startsWith(gene, "chrX")) %>% filter(cov > 0)
pcg_cn <- read_tsv(args$pcg_cn, show_col_types = FALSE)
gc_pcg <- read_tsv(args$gc_pcg, show_col_types = FALSE)

# Perform calculations
normalizer <- cov_scg %>% summarise(mean_cov = mean(cov)) %>% pull()
cov_scg_normalized <- cov_scg %>% mutate(cov = cov / normalizer)
gc_table <- inner_join(gc_scg, cov_scg_normalized, by = "gene")

fitting <- function(gc_data) {
  fit <- lm(cov ~ poly(gc, 2, raw = TRUE), data = gc_data)
  coefficients <- coef(fit)
  return(coefficients)
}

coefficients <- fitting(gc_table)
a <- coefficients[1]
b <- coefficients[2]
c <- coefficients[3]

correlation_plot <- ggplot(gc_table, aes(x = gc, y = cov)) +
  geom_point() +
  geom_smooth(method = "lm", color="grey", se=T, formula = y~poly(x,2)) +
  xlab("gc") + ylab("cov")+
  geom_text(data = data.frame(gc = mean(gc_table$gc), cov = mean(gc_table$cov), 
                              label = paste0("cov = ", as.character(round((a),5)), " + (", as.character(round((b),5)), ") x gc + (", as.character(round((c),5)), ") x gc^2")),
            aes(label = label), hjust = 0.7, vjust = -30)

normalized_pcg <- pcg_cn %>%
  inner_join(gc_pcg, by = "gene") %>%
  mutate(expected_1_copy = (a + (b * gc) + c * (gc * gc)),
         diff_from_1 = 1 - expected_1_copy,
         gc_copynumber = round((copynumber + (diff_from_1 * copynumber)), 2)) %>%
  select(-expected_1_copy, -scg_mean, -diff_from_1)

# Write the output files
basename <- file_path_sans_ext(basename(args$pcg_cn))

write_tsv(normalized_pcg, paste0(args$output, "/", basename, "-gcbias.tsv"))
ggsave(correlation_plot, file=paste0(args$output, "/", basename, "-gcbias.png"), dpi=300)

