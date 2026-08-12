
# Cancer Chromosome Count, Display and Cross-Morphology Comparison ----


# 0. LIBRARY LOADING ----

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)

# 1. DATA PREPARATION ----

# The csv file below was downloaded from the Mitelman Database's "Cases Cytogenetics Searcher", found here: "https://mitelmandatabase.isb-cgc.org/case_search" 

# The searcher may be used to retrieve similar files for other cancer types and morphologies.

# Ensure the working directory is set to this script's file location, and the csv file to be studied is both specified and in the relevant directory.

infile <- "csv/Ast_Chrom_Numbers.csv"
df <- read.csv(infile, stringsAsFactors = FALSE)


# Trim any NA entries in our column of concern (KaryShort)
df <- df %>% filter(!is.na(KaryShort), str_trim(KaryShort) != "")



# If samples have clones, split each case's KaryShort into one row per clone,
clones <- df %>%
  mutate(row_id = row_number()) %>%
  separate_rows(KaryShort, sep = "/") %>%
  mutate(KaryShort = str_trim(KaryShort))


# Extract the chromosome count number for the KaryShort column
clones <- clones %>%
  mutate(modal_token = str_extract(KaryShort, "^[^,]+"))


# Input sanitization 1: Converts chromosome ranges into paired numberic entries
average_modal <- function(token) {
  if (str_detect(token, "^\\d+$")) {
    return(as.numeric(token))
  }
  m <- str_match(token, "^(\\d+)-(\\d+)$")
  if (!is.na(m[1, 1])) {
    lo <- as.numeric(m[1, 2])
    hi <- as.numeric(m[1, 3])
    return((lo + hi) / 2)
  }
  return(NA_real_)
}

# Input sanitization 2: Averages numeric ranges into a single entry, then floors it to an integer to avoid introducing artificial half-chromosomes into the data.
chrom_counts <- clones %>%
  rowwise() %>%
  mutate(chrom_n = floor(average_modal(modal_token))) %>%
  ungroup()

#Input Sanitization 3: Flag any entries dropped from the data set.
n_dropped <- sum(is.na(chrom_counts$chrom_n))
if (n_dropped > 0) {
  message(n_dropped, " sample(s) could not be parsed and were dropped: ",
          paste(unique(chrom_counts$modal_token[is.na(chrom_counts$chrom_n)]), collapse = ", "))
}
chrom_counts <- chrom_counts %>% filter(!is.na(chrom_n))


# 2. MORPHOLOGY GROUPING ----

# Bin into whole-chromosome-number counts, then convert each Morph group's counts into a of that Morph group's total population, allowing direct comparison despite different sample sizes.
binned <- chrom_counts %>%
  mutate(chrom_bin = round(chrom_n)) %>%
  count(Morph, chrom_bin, name = "n") %>%
  group_by(Morph) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  ungroup()


# 3. PLOTTING ----

# Phenotype Plot, all groups overlaid on the same axes, with markers at 2c, 3c, and 4c counts.
ggplot(binned, aes(x = chrom_bin, y = pct, fill = Morph)) +
  geom_col(position = "identity", alpha = 0.55, width = 1) +
  geom_vline(xintercept = 46, linetype = "dashed", color = "red") +
  geom_vline(xintercept = 69, linetype = "dashed", color = "violet") +
  geom_vline(xintercept = 92, linetype = "dashed", color = "blue") +
  xlim(20,125) +
  labs(
    title = "Astrocytoma Chromosome Number Distribution by Tumor Phenotype",
    subtitle = "Marker lines placed at chromsome counts 46 (2c), 69 (3c), and 92 (4c)",
    x = "Chromosome number (ranges averaged)",
    y = "Percent of sample population within morphology",
    fill = "Phenotype"
  ) +
  theme_minimal() +
  theme(legend.position = "top")



# Plot ignoring phenotype.
ggplot(chrom_counts, aes(x = chrom_n)) +
  geom_histogram(binwidth = 1, fill = "#4C72B0", color = "white") +
  geom_vline(xintercept = 46, linetype = "dashed", color = "red") +
  geom_vline(xintercept = 69, linetype = "dashed", color = "violet") +
  geom_vline(xintercept = 92, linetype = "dashed", color = "blue") +
  xlim(20,125) +
  labs(
    title = "Astrocytoma Chromosome Number Distribution",
    subtitle = "Marker lines placed at chromsome counts 46 (2c), 69 (3c), and 92 (4c)",
    x = "Chromosome number (ranges averaged)",
    y = "Sample count",
    fill = "Phenotype"
  ) +
  theme_minimal() +
  theme(legend.position = "top")


# 4. Data Interpretation ====

# Tally of cancers by outcome
outcome_count <- count(chrom_counts, Morph)
print(outcome_count)

out55_count <- count(chrom_counts, Morph, chrom_n > 54)
print(out55_count)


ast_ben_pct <- (out55_count[2,3] / outcome_count[1,2])
print((paste("Fraction of total Astrocytoma (grade I-II) samples over 55 chromosomes:",ast_ben_pct)))

ast_mal_pct <- (out55_count[4,3] / outcome_count[2,2])
print(paste("Fraction of total Astrocytoma (grade III-IV/Glioblastoma) samples over 55 chromosomes:", ast_mal_pct))

ast_comp <- (ast_mal_pct / ast_ben_pct)
print(paste("Relative fraction of malign samples with 55+ chromosomes over benign samples with 55+ chromosomes:", ast_comp))




out75_count <- count(chrom_counts, Morph, chrom_n > 74)
print(out75_count)


ast_ben_pct2 <- (out75_count[2,3] / outcome_count[1,2])
print((paste("Fraction of total Astrocytoma (grade I-II) samples over 75 chromosomes:",ast_ben_pct2)))

ast_mal_pct2 <- (out75_count[4,3] / outcome_count[2,2])
print(paste("Fraction of total Astrocytoma (grade III-IV/Glioblastoma) samples over 75 chromosomes:", ast_mal_pct2))

ast_comp2 <- (ast_mal_pct2 / ast_ben_pct2)
print(paste("Relative fraction of malign samples with 75+ chromosomes over benign samples with 75+ chromosomes:", ast_comp2))
