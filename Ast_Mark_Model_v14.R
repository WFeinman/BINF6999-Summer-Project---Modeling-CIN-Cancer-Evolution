
# Dual-Population Markov Chain Model for Chromosome Missegregation  Incorporating WGD, Binomial Missegregation, Asymmetric Fitness, & Parallelization ----

# 0. LIBRARY LOADING ----

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(foreach)
library(doParallel)


# 1. DATA PREPARATION ----

# The csv file below was downloaded from the Mitelman Database's "Cases Cytogenetics Searcher", found here: "https://mitelmandatabase.isb-cgc.org/case_search" 

# The searcher may be used to retrieve similar files for other cancer types and morphologies.

# Ensure the working directory is set to this script's file location, and the csv file to be studied is both specified and in the relevant directory.


# Reading input file into R
infile <- "csv/Ast_Chrom_Numbers.csv"
df <- read.csv(infile, stringsAsFactors = FALSE)

# Trim any NA entries in our column of concern (KaryShort)
df <- df %>% filter(!is.na(KaryShort), str_trim(KaryShort) != "")

# Extract the chromosome count number for the KaryShort column
clones <- df %>%
  mutate(row_id = row_number()) %>%
  separate_rows(KaryShort, sep = "/") %>%
  mutate(KaryShort = str_trim(KaryShort)) %>%
  mutate(modal_token = str_extract(KaryShort, "^[^,]+"))

# Input sanitization 1: Converts chromosome ranges into paired numberic entries
average_modal <- function(token) {
  if (str_detect(token, "^\\d+$")) return(as.numeric(token))
  m <- str_match(token, "^(\\d+)-(\\d+)$")
  if (!is.na(m[1, 1])) return((as.numeric(m[1, 2]) + as.numeric(m[1, 3])) / 2)
  return(NA_real_)
}

# Input sanitization 2: Averages numeric ranges into a single entry, then floors it to an integer to avoid introducing artificial half-chromosomes into the data.
chrom_counts <- clones %>%
  rowwise() %>%
  mutate(chrom_n = floor(average_modal(modal_token))) %>%
  ungroup() %>%
  filter(!is.na(chrom_n))

#Input Sanitization 3: Flag any entries dropped from the data set.
n_dropped <- sum(is.na(chrom_counts$chrom_n))
if (n_dropped > 0) {
  message(n_dropped, " sample(s) could not be parsed and were dropped: ",
          paste(unique(chrom_counts$modal_token[is.na(chrom_counts$chrom_n)]), collapse = ", "))
}

# Get the observed distribution
obs_binned <- chrom_counts %>%
  mutate(chrom_bin = round(chrom_n)) %>%
  count(chrom_bin, name = "n") %>%
  mutate(pct = n / sum(n))


# 2. MARKOV CHAIN SETUP WITH ASYMMETRIC INDEPENDENT FITNESS ----

N_MAX <- 150 # Maximum chromosome number to track

# Function to build the transition matrix for missegregation (Binomial)
build_transition_matrix <- function(p_misseg, N_max) {
  M <- matrix(0, nrow = N_max, ncol = N_max)
  
  for (i in 1:N_max) {
    M[i, i] <- dbinom(0, i, p_misseg)
    
    for (k in 1:i) {
      prob <- dbinom(k, i, p_misseg) / 2
      
      if (i + k <= N_max) M[i, i + k] <- M[i, i + k] + prob
      else M[i, N_max] <- M[i, N_max] + prob
      
      if (i - k >= 1) M[i, i - k] <- M[i, i - k] + prob
      else M[i, 1] <- M[i, 1] + prob 
    }
  }
  return(M)
}

# Simulate the dual population with asymmetric fitness peak variables
simulate_MC <- function(pA, pB, w, sigmaA_gain, sigmaA_loss, sigmaB_gain, sigmaB_loss, peakA, peakB, t_max, N_max = N_MAX) {
  M_A <- build_transition_matrix(pA, N_max)
  M_B <- build_transition_matrix(pB, N_max)
  
  # Define independent Gaussian fitness landscapes F(x) with asymmetric gains/losses
  x_vals <- 1:N_max
  
  F_A <- ifelse(x_vals > peakA,
                exp(-0.5 * ((x_vals - peakA) / sigmaA_gain)^2),
                exp(-0.5 * ((x_vals - peakA) / sigmaA_loss)^2))
                
  F_B <- ifelse(x_vals > peakB,
                exp(-0.5 * ((x_vals - peakB) / sigmaB_gain)^2),
                exp(-0.5 * ((x_vals - peakB) / sigmaB_loss)^2))
  
  V_A <- rep(0, N_max)
  V_B <- rep(0, N_max)
  V_A[46] <- 1 # Assuming the initial population still starts at normal diploid (46)
  
  for (step in 1:t_max) {
    # a) Whole Genome Doubling (WGD)
    WGD_mass <- rep(0, N_max)
    for (i in 1:floor(N_max/2)) {
      WGD_mass[2 * i] <- V_A[i] * w
    }
    V_A_rem <- V_A * (1 - w) 
    
    # b) Missegregation (Matrix multiplication)
    V_A_miss <- as.numeric(V_A_rem %*% M_A)
    V_B_miss <- as.numeric(V_B %*% M_B) + WGD_mass
    
    # v) Apply Independent Fitness F(x) & Normalize
    V_A_fit <- V_A_miss * F_A
    V_B_fit <- V_B_miss * F_B
    
    total_mass <- sum(V_A_fit) + sum(V_B_fit)
    
    V_A <- V_A_fit / total_mass
    V_B <- V_B_fit / total_mass
  }
  
  final_dist <- V_A + V_B
  return(final_dist / sum(final_dist))
}


# 3. MAXIMUM LIKELIHOOD ESTIMATION, GRID SETUP ----


calc_log_likelihood <- function(obs_data, pred_dist) {
  ll <- 0
  #1e^-3 used for the epsilon parameter, as it best matches the smallest possibility this dataset can distinguish from zero (1/939).
  for (i in 1:nrow(obs_data)) {
    x <- obs_data$chrom_bin[i]
    count <- obs_data$n[i]
    if (x >= 1 && x <= length(pred_dist)) {
      p <- pred_dist[x] + epsilon
      ll <- ll + count * log(p)
    }
  }
  return(ll)
}


# Note: The specified paramater range below is specific to the Mitelman Astrocytoma data, and represent a compromise between a fully exhaustive grid search and ease of reproducibility without a large computing cluster. This example grid search iterates through 165888 parameter combinations. On a single 4.70 GHz CPU using 23 threads, this took about 55 minutes to perform. Increasing the parameter range or fidelity will result in a better fit to the data, but take exponentially more computational resources.

# pA: Missegregation chance per chromosome of Group A cells.
# pB: Missegregation chance per chromosome of Group B cells.
# w: Chance of Group A cells to undergo WGD, doubling their chromosome count and transitioning to Group B.

# peakA: Chromosome count for optimal fitness of Group A cells.
# peakA: Chromosome count for optimal fitness of Group B cells.

# sigmaA/B/_gain/loss: Selection parameter representing directional tolerance to copy number change. Small number = Sharp fitness dropoff. Large number = Shallow fitness dropoff.

# sigmaA_gain: Selection strength acting on chromosome counts above peakA. 
# sigmaA_loss: Selection strength acting on chromosome counts below peakA.

# sigmaB_gain: Selection strength acting on chromosome counts above peakB.
# sigmaB_gain: Selection strength acting on chromosome counts below peakB.

# t: The number of cycles the model runs for.


# Expand Grid Setup, Set Parameters
grid <- expand.grid(
  pA = seq(0.005, 0.02, length.out = 3),   
  pB = seq(0.01, 0.04, length.out = 3),    
  w  = seq(0.001, 0.002, length.out = 3),   
  peakA  = c(45.5),              
  peakB  = seq(75, 90, 3),
  sigmaA_gain = c(2,3),                  
  sigmaA_loss = c(3.5, 5),                  
  sigmaB_gain = seq(1, 61, 4),                 
  sigmaB_loss = seq(1, 61, 4),                 
                  
  t  = (5000)                      
)

# 4.  PARALLEL GRID SEARCH ----
# Set up cluster to use one less than total available cores, keep the system responsive
num_cores <- parallel::detectCores() - 1 
cl <- parallel::makeCluster(num_cores)
doParallel::registerDoParallel(cl)

print(paste("Starting Parallel Grid Search over", nrow(grid), "parameter combinations using", num_cores, "cores..."))

# Execute the parallel grid search
# .combine = rbind ensures the results are combined into a clean matrix row-by-row
search_results <- foreach(
  i = 1:nrow(grid), 
  .combine = rbind, 
  .export = c("simulate_MC", "build_transition_matrix", "calc_log_likelihood", "N_MAX", "obs_binned")
) %dopar% {
  
  params <- grid[i, ]
  
  pred_dist <- simulate_MC(
    params$pA, params$pB, params$w, 
    params$sigmaA_gain, params$sigmaA_loss,
    params$sigmaB_gain, params$sigmaB_loss, 
    params$peakA, params$peakB, 
    params$t, N_MAX
  )
  
  ll <- calc_log_likelihood(obs_binned, pred_dist)
  
  # Return the index and the calculated log-likelihood 
  c(index = i, log_likelihood = ll)
}

# Close the cluster
parallel::stopCluster(cl)

# 5. RESULT PROCESSING ----
# Locate the row with the maximum log-likelihood
best_idx <- search_results[which.max(search_results[, "log_likelihood"]), "index"]
best_ll <- search_results[which.max(search_results[, "log_likelihood"]), "log_likelihood"]
best_params <- grid[best_idx, ]

print("Best Parameters Found:")
print(best_params)
print(paste("Maximum Log-Likelihood:", best_ll))

# Re-run the model once with the absolute best parameters to get the final line for graphing
best_pred <- simulate_MC(
  best_params$pA, best_params$pB, best_params$w, 
  best_params$sigmaA_gain, best_params$sigmaA_loss,
  best_params$sigmaB_gain, best_params$sigmaB_loss, 
  best_params$peakA, best_params$peakB, 
  best_params$t, N_MAX
)


# 6. PLOTTING ----


pred_df <- data.frame(
  chrom_bin = 1:N_MAX,
  pct = best_pred
)

ggplot() +
  geom_col(data = obs_binned, aes(x = chrom_bin, y = pct, fill = "Observed Data"), 
           alpha = 0.5, width = 1) +
  geom_line(data = pred_df, aes(x = chrom_bin, y = pct, color = "Model Prediction"), 
            linewidth = 1.2) +
  geom_vline(xintercept = c(best_params$peakA, best_params$peakB), linetype = "dashed", color = "black", alpha = 0.5) +
  geom_vline(xintercept = 46, linetype = "dashed", color = "red") +
  geom_vline(xintercept = 69, linetype = "dashed", color = "violet") +
  geom_vline(xintercept = 92, linetype = "dashed", color = "blue") +
  scale_fill_manual(name = "", values = c("Observed Data" = "#4C72B0")) +
  scale_color_manual(name = "", values = c("Model Prediction" = "red")) +
  xlim(20, 125) +
  labs(
    title = "Markov Chain Model Fit vs. Observed Astrocytoma Data",
    subtitle = sprintf("Best parameters: pA=%.3f, pB=%.3f, w=%.3f, sAg=%.1f, sAl=%.1f, sBg=%.1f, sBl=%.1f, pkA=%.1f, pkB=%.1f, t=%.0f", 
                       best_params$pA, best_params$pB, best_params$w, 
                       best_params$sigmaA_gain, best_params$sigmaA_loss,
                       best_params$sigmaB_gain, best_params$sigmaB_loss, 
                       best_params$peakA, best_params$peakB, best_params$t),
    x = "Chromosome Number",
    y = "Proportion of Population"
  ) +
  theme_minimal() +
  theme(legend.position = "top")
