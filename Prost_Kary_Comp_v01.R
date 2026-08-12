
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

infile <- "csv/Prostrate_adeno_chrom_num.csv"
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



# Plot ignoring phenotype, all groups overlaid on the same axes, with markers at 2c, 3c, and 4c counts.
ggplot(chrom_counts, aes(x = chrom_n)) +
  geom_histogram(binwidth = 1, fill = "purple4", color = "white") +
  geom_vline(xintercept = 46, linetype = "dashed", color = "red") +
  geom_vline(xintercept = 69, linetype = "dashed", color = "violet") +
  geom_vline(xintercept = 92, linetype = "dashed", color = "blue") +
  xlim(20,125) +
  labs(
    title = "Prostrate Chromosome Number Distribution",
    subtitle = "Marker lines placed at chromsome counts 46 (2c), 69 (3c), and 92 (4c)",
    x = "Chromosome number (ranges averaged)",
    y = "Sample count",
    fill = "Phenotype"
  ) +
  theme_minimal() +
  theme(legend.position = "top")



