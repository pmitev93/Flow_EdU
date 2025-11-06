# ==============================================================================
# FLOW CYTOMETRY ANALYSIS - MASTER SCRIPT
# ==============================================================================

# Libraries ####
library(flowCore)
library(tidyverse)
library(sp) 

# Working Directory ####
setwd("/Users/petar.mitev/Library/CloudStorage/OneDrive-KarolinskaInstitutet/Experiments/Flow/R_New")

# Paths ####
master_path <- "Experiments/"
#experiment_folders <- list.dirs(master_path, recursive = FALSE, full.names = TRUE)

cat("Found", length(experiment_folders), "experiments\n")
print(experiment_folders)

# Channel Mapping ####
CHANNELS <- list(
  FSC_A = "FSC-A",
  FSC_H = "FSC-H",
  SSC_A = "SSC-A",
  EdU = "FL1-A",      # EdU B525-A
  HA = "FL3-A",       # HA R660-A
  DCM = "FL5-A",      # DCM R780-A
  FxCycle = "FL6-A",  # FxCycle V450-A
  SAMD9 = "FL10-A"    # SAMD9/SAMD9L Y585-A (alternative to HA)
)

# Gate Definitions ----

GATES <- list()

# Debris removal gate (FSC-A vs SSC-A)
GATES$debris <- matrix(c(
  1.4e6, 4.1e6,
  4.0e6, 9.4e6,
  8.7e6, 14.0e6,
  12.0e6, 14.0e6,
  13.0e6, 13.0e6,
  14.0e6, 11.0e6,
  13.0e6, 7.3e6,
  8.7e6, 2.4e6,
  4.8e6, 0.384e6,
  2.7e6, 0.068e6,
  1.8e6, 0.565e6,
  1.1e6, 1.7e6,
  1.4e6, 4.1e6  # Close the polygon
), ncol = 2, byrow = TRUE)
colnames(GATES$debris) <- c("FSC-A", "SSC-A")

# Singlet gate (FSC-A vs FSC-H)
GATES$singlet <- matrix(c(
  5.8e6, 2.5e6,
  9.5e6, 2.5e6,
  6.4e6, 1.3e6,
  4.5e6, 0.740e6,
  3.0e6, 0.461e6,
  1.7e6, 0.489e6,
  1.5e6, 0.525e6,
  1.6e6, 0.658e6,
  1.9e6, 0.987e6,
  3.3e6, 1.6e6,
  5.8e6, 2.5e6  # Close the polygon
), ncol = 2, byrow = TRUE)
colnames(GATES$singlet) <- c("FSC-A", "FSC-H")

# Dead cell removal gate (DCM-A vs SSC-A)
GATES$live_cells <- matrix(c(
  122, 181000,
  18000, 181000,
  18000, 13000000,
  122, 13000000,
  122, 181000  # Close the rectangle
), ncol = 2, byrow = TRUE)
colnames(GATES$live_cells) <- c("DCM-A", "SSC-A")

# S-phase gating - Stage 1: Remove outliers (FxCycle-A vs EdU-A)
GATES$s_phase_outliers <- matrix(c(
  1.5e6, 100,
  9e6, 100,
  9e6, 5e6,
  1.5e6, 5e6,
  1.5e6, 100  # Close the rectangle
), ncol = 2, byrow = TRUE)
colnames(GATES$s_phase_outliers) <- c("FxCycle-A", "EdU-A")

# S-phase gating - Stage 2: FxCycle quantile filter (dynamic, 1%-99%)
# This gate is calculated dynamically for each sample
GATES$fxcycle_quantile <- list(
  type = "quantile_range",
  parameter = CHANNELS$FxCycle,  # Use the channel mapping we defined earlier
  probs = c(0.01, 0.90),
  description = "Remove FxCycle outliers (keep 1st-90th percentile)"
)

# Gating Functions ####

# Apply quantile-based gate
apply_quantile_gate <- function(fcs_data, gate_params) {
  # Extract the channel values
  channel_name <- gate_params$parameter
  values <- exprs(fcs_data)[, channel_name]
  
  # Calculate quantile bounds
  quantile_limits <- quantile(values, probs = gate_params$probs, na.rm = TRUE)
  lower_bound <- quantile_limits[1]
  upper_bound <- quantile_limits[2]
  
  # Create filter
  filter_logical <- values >= lower_bound & values <= upper_bound
  
  # Apply filter
  gated_data <- Subset(fcs_data, filter_logical)
  
  # Return results with metadata (use simple names)
  return(list(
    data = gated_data,
    n_before = nrow(fcs_data),
    n_after = nrow(gated_data),
    n_removed = nrow(fcs_data) - nrow(gated_data),
    bounds = c(lower = unname(lower_bound), upper = unname(upper_bound))
  ))
}

# S-phase gating - Stage 3: Top 50% EdU + FxCycle range (dynamic)
# Select cells with highest EdU incorporation AND appropriate FxCycle range
GATES$edu_fxcycle_sphase <- list(
  type = "dual_quantile",
  edu_parameter = CHANNELS$EdU,
  edu_prob = 0.50,  # Keep cells above 50th percentile
  fxcycle_parameter = CHANNELS$FxCycle,
  fxcycle_probs = c(0.01, 0.90),  # Keep between 1st-90th percentile
  description = "Select top 50% EdU AND FxCycle 1st-90th percentile"
)


# Apply dual quantile gate (EdU threshold + FxCycle range)
apply_dual_quantile_gate <- function(fcs_data, gate_params) {
  # EdU threshold (top 50%)
  edu_values <- exprs(fcs_data)[, gate_params$edu_parameter]
  edu_threshold <- quantile(edu_values, probs = gate_params$edu_prob, na.rm = TRUE)
  edu_filter <- edu_values >= edu_threshold
  
  # FxCycle range (1%-90%)
  fxcycle_values <- exprs(fcs_data)[, gate_params$fxcycle_parameter]
  fxcycle_quantiles <- quantile(fxcycle_values, probs = gate_params$fxcycle_probs, na.rm = TRUE)
  fxcycle_filter <- fxcycle_values >= fxcycle_quantiles[1] & fxcycle_values <= fxcycle_quantiles[2]
  
  # Combine both filters
  combined_filter <- edu_filter & fxcycle_filter
  
  # Apply filter
  gated_data <- Subset(fcs_data, combined_filter)
  
  return(list(
    data = gated_data,
    n_before = nrow(fcs_data),
    n_after = nrow(gated_data),
    n_removed = nrow(fcs_data) - nrow(gated_data),
    edu_threshold = unname(edu_threshold),
    fxcycle_bounds = unname(fxcycle_quantiles)
  ))
}

# Apply quantile threshold gate (keep cells above a percentile)
apply_quantile_threshold <- function(fcs_data, gate_params) {
  # Extract the channel values
  channel_name <- gate_params$parameter
  values <- exprs(fcs_data)[, channel_name]
  
  # Calculate threshold
  threshold <- quantile(values, probs = gate_params$prob, na.rm = TRUE)
  
  # Create filter (keep cells >= threshold)
  filter_logical <- values >= threshold
  
  # Apply filter
  gated_data <- Subset(fcs_data, filter_logical)
  
  # Return results with metadata
  return(list(
    data = gated_data,
    n_before = nrow(fcs_data),
    n_after = nrow(gated_data),
    n_removed = nrow(fcs_data) - nrow(gated_data),
    threshold = unname(threshold)
  ))
}

# HA gating (dynamic, per-experiment)
# Threshold set by 98th percentile of empty vector control
GATES$ha_positive <- list(
  type = "control_based_threshold",
  parameter = CHANNELS$HA,
  control_pattern = "Empty_Vector",  # Pattern to find control sample
  prob = 0.98,
  description = "HA-positive cells (above 98th percentile of empty vector control)"
)

# Find control sample in an experiment's metadata
# Find control sample in an experiment's metadata
find_control_sample <- function(metadata, pattern = "Empty_Vector_Dox-") {
  # Look for pattern in sample names (case-insensitive, flexible matching)
  matches <- grep(pattern, metadata$sample_name, ignore.case = TRUE)
  
  if(length(matches) == 0) {
    warning(sprintf("No control sample found matching pattern '%s'", pattern))
    return(NULL)
  }
  
  if(length(matches) > 1) {
    match_info <- paste(sprintf("%s (%s)", 
                                metadata$well[matches], 
                                metadata$sample_name[matches]), 
                        collapse = ", ")
    warning(sprintf("Multiple control samples found matching '%s': %s. Using first one.", 
                    pattern, match_info))
  }
  
  cat(sprintf("Using control sample: %s (%s)\n", 
              metadata$well[matches[1]], 
              metadata$sample_name[matches[1]]))
  
  return(matches[1])  # Return index of control sample
}

# Calculate threshold from control sample
calculate_control_threshold <- function(control_fcs, gate_params) {
  channel_name <- gate_params$parameter
  values <- exprs(control_fcs)[, channel_name]
  
  threshold <- quantile(values, probs = gate_params$prob, na.rm = TRUE)
  
  return(unname(threshold))
}

# Apply control-based threshold gate
apply_control_threshold <- function(fcs_data, threshold, channel_name) {
  values <- exprs(fcs_data)[, channel_name]
  filter_logical <- values >= threshold
  gated_data <- Subset(fcs_data, filter_logical)
  
  return(list(
    data = gated_data,
    n_before = nrow(fcs_data),
    n_after = nrow(gated_data),
    n_removed = nrow(fcs_data) - nrow(gated_data),
    threshold = threshold
  ))
}

# Data Loading Functions ####

# Load all FCS files from a single experiment folder
load_experiment <- function(experiment_path) {
  experiment_name <- basename(experiment_path)
  cat(sprintf("Loading experiment: %s\n", experiment_name))
  
  # Get all FCS files in the folder
  fcs_files <- list.files(experiment_path, pattern = "\\.fcs$", 
                          full.names = TRUE, ignore.case = TRUE)
  
  if(length(fcs_files) == 0) {
    warning(sprintf("No FCS files found in %s", experiment_path))
    return(NULL)
  }
  
  cat(sprintf("  Found %d FCS files\n", length(fcs_files)))
  
  # Parse filenames to get metadata
  metadata <- data.frame(
    full_filename = basename(fcs_files),
    file_path = fcs_files,
    stringsAsFactors = FALSE
  )
  
  # Extract well and sample info from filenames
  # Format: 01-Well-E9_624_SAMD9_I1553T_Dox+.fcs
  metadata$well <- str_extract(metadata$full_filename, "[A-H][0-9]+")
  
  # Extract sample name (everything after first underscore, before .fcs)
  metadata$sample_name <- str_split_fixed(metadata$full_filename, "_", 2)[, 2]
  metadata$sample_name <- gsub("\\.fcs$", "", metadata$sample_name, ignore.case = TRUE)
  
  metadata$experiment <- experiment_name
  
  # Load FCS files into flowSet
  fs <- read.flowSet(fcs_files, transformation = FALSE, truncate_max_range = FALSE)
  sampleNames(fs) <- metadata$full_filename
  
  # Apply compensation if available
  # Check if compensation matrix exists in first file
  comp_matrix <- keyword(fs[[1]])$SPILL
  if(is.null(comp_matrix)) {
    comp_matrix <- keyword(fs[[1]])$`$SPILLOVER`
  }
  
  if(!is.null(comp_matrix)) {
    cat("  Applying compensation matrix\n")
    fs <- compensate(fs, comp_matrix)
  } else {
    cat("  No compensation matrix found\n")
  }
  
  return(list(
    flowset = fs,
    metadata = metadata,
    experiment_name = experiment_name,
    n_samples = length(fcs_files)
  ))
}

# Load all experiments from master folder
load_all_experiments <- function(experiment_folders) {
  cat(sprintf("\n=== Loading %d experiments ===\n\n", length(experiment_folders)))
  
  all_experiments <- list()
  
  for(i in seq_along(experiment_folders)) {
    exp_data <- load_experiment(experiment_folders[i])
    if(!is.null(exp_data)) {
      all_experiments[[exp_data$experiment_name]] <- exp_data
    }
  }
  
  cat(sprintf("\n=== Successfully loaded %d experiments ===\n", length(all_experiments)))
  
  # Print summary
  cat("\nExperiment Summary:\n")
  for(exp_name in names(all_experiments)) {
    cat(sprintf("  %s: %d samples\n", 
                exp_name, 
                all_experiments[[exp_name]]$n_samples))
  }
  
  return(all_experiments)
}

# Calculate EdU vs HA correlation (log10 scale)
calculate_edu_ha_correlation <- function(fcs_data, edu_channel, ha_channel) {
  # Extract values
  edu_values <- exprs(fcs_data)[, edu_channel]
  ha_values <- exprs(fcs_data)[, ha_channel]
  
  # Apply log10 transformation
  # Add small value to avoid log10(0) = -Inf
  edu_log <- log10(edu_values + 1)
  ha_log <- log10(ha_values + 1)
  
  # Calculate Pearson correlation
  correlation <- cor(ha_log, edu_log, use = "complete.obs")
  
  # Also return the data for plotting
  return(list(
    correlation = correlation,
    n_cells = length(edu_values),
    ha_log = ha_log,
    edu_log = edu_log
  ))
}







# Full Analysis Pipeline ####

# Apply all gates sequentially to a single sample ----
process_single_sample <- function(fcs_data, sample_name, ha_threshold, gates = GATES, channels = CHANNELS) {
  
  cat(sprintf("\nProcessing: %s\n", sample_name))
  
  # Initialize results tracking
  results <- list(
    sample_name = sample_name,
    cell_counts = list(),
    gates_applied = list()
  )
  
  # Starting count
  results$cell_counts$initial <- nrow(fcs_data)
  current_data <- fcs_data
  
  # Gate 1: Debris removal (FSC-A vs SSC-A)
  cat("  Gate 1: Debris removal...")
  debris_filter <- point.in.polygon(
    exprs(current_data)[, channels$FSC_A],
    exprs(current_data)[, channels$SSC_A],
    gates$debris[, 1],  # FSC-A column
    gates$debris[, 2]   # SSC-A column
  ) > 0
  current_data <- Subset(current_data, debris_filter)
  results$cell_counts$after_debris <- nrow(current_data)
  cat(sprintf(" %s cells\n", format(nrow(current_data), big.mark = ",")))
  
  # Gate 2: Singlets (FSC-A vs FSC-H)
  cat("  Gate 2: Singlets...")
  singlet_filter <- point.in.polygon(
    exprs(current_data)[, channels$FSC_A],
    exprs(current_data)[, channels$FSC_H],
    gates$singlet[, 1],  # FSC-A column
    gates$singlet[, 2]   # FSC-H column
  ) > 0
  current_data <- Subset(current_data, singlet_filter)
  results$cell_counts$after_singlets <- nrow(current_data)
  cat(sprintf(" %s cells\n", format(nrow(current_data), big.mark = ",")))
  
  # Gate 3: Live cells (DCM-A vs SSC-A)
  cat("  Gate 3: Live cells...")
  live_filter <- point.in.polygon(
    exprs(current_data)[, channels$DCM],
    exprs(current_data)[, channels$SSC_A],
    gates$live_cells[, 1],  # DCM-A column
    gates$live_cells[, 2]   # SSC-A column
  ) > 0
  current_data <- Subset(current_data, live_filter)
  results$cell_counts$after_live <- nrow(current_data)
  cat(sprintf(" %s cells\n", format(nrow(current_data), big.mark = ",")))
  
  # Gate 4: S-phase outlier removal (FxCycle vs EdU)
  cat("  Gate 4: S-phase outliers...")
  sphase_filter <- point.in.polygon(
    exprs(current_data)[, channels$FxCycle],
    exprs(current_data)[, channels$EdU],
    gates$s_phase_outliers[, 1],  # FxCycle-A column
    gates$s_phase_outliers[, 2]   # EdU-A column
  ) > 0
  current_data <- Subset(current_data, sphase_filter)
  results$cell_counts$after_sphase_outliers <- nrow(current_data)
  cat(sprintf(" %s cells\n", format(nrow(current_data), big.mark = ",")))
  
  # Gate 5: FxCycle quantile (1%-99%)
  cat("  Gate 5: FxCycle quantile...")
  fxcycle_result <- apply_quantile_gate(current_data, gates$fxcycle_quantile)
  current_data <- fxcycle_result$data
  results$cell_counts$after_fxcycle_quantile <- nrow(current_data)
  results$gates_applied$fxcycle_bounds <- fxcycle_result$bounds
  cat(sprintf(" %s cells\n", format(nrow(current_data), big.mark = ",")))
  
  # Gate 6: Top 50% EdU + FxCycle 1-90%
  cat("  Gate 6: Top 50% EdU + FxCycle range...")
  edu_fxcycle_result <- apply_dual_quantile_gate(current_data, gates$edu_fxcycle_sphase)
  current_data <- edu_fxcycle_result$data
  results$cell_counts$after_edu_fxcycle <- nrow(current_data)
  results$gates_applied$edu_threshold <- edu_fxcycle_result$edu_threshold
  results$gates_applied$fxcycle_bounds_stage3 <- edu_fxcycle_result$fxcycle_bounds
  cat(sprintf(" %s cells\n", format(nrow(current_data), big.mark = ",")))
  
  # Gate 7: HA-positive
  cat("  Gate 7: HA-positive...")
  ha_result <- apply_control_threshold(current_data, ha_threshold, channels$HA)
  current_data <- ha_result$data
  results$cell_counts$after_ha_positive <- nrow(current_data)
  results$gates_applied$ha_threshold <- ha_threshold
  cat(sprintf(" %s cells\n", format(nrow(current_data), big.mark = ",")))
  
  # Final: Calculate EdU vs HA correlation
  if(nrow(current_data) > 10) {  # Need at least some cells for correlation
    cat("  Calculating correlation...")
    corr_result <- calculate_edu_ha_correlation(current_data, channels$EdU, channels$HA)
    results$correlation <- corr_result$correlation
    results$correlation_n_cells <- corr_result$n_cells
    results$final_data <- current_data
    results$correlation_data <- list(ha_log = corr_result$ha_log, edu_log = corr_result$edu_log)
    cat(sprintf(" r = %.4f\n", corr_result$correlation))
  } else {
    warning("Too few cells for correlation")
    results$correlation <- NA
    results$correlation_n_cells <- 0
  }
  
  return(results)
}







# Process control sample to get HA threshold
calculate_ha_threshold_from_control <- function(control_fcs, control_name, gates = GATES, channels = CHANNELS) {
  
  cat(sprintf("\nProcessing control sample: %s\n", control_name))
  
  current_data <- control_fcs
  
  # Apply all gates EXCEPT HA gate
  # Gate 1: Debris
  debris_filter <- point.in.polygon(
    exprs(current_data)[, channels$FSC_A],
    exprs(current_data)[, channels$SSC_A],
    gates$debris[, 1], gates$debris[, 2]
  ) > 0
  current_data <- Subset(current_data, debris_filter)
  cat(sprintf("  After debris: %s cells\n", format(nrow(current_data), big.mark = ",")))
  
  # Gate 2: Singlets
  singlet_filter <- point.in.polygon(
    exprs(current_data)[, channels$FSC_A],
    exprs(current_data)[, channels$FSC_H],
    gates$singlet[, 1], gates$singlet[, 2]
  ) > 0
  current_data <- Subset(current_data, singlet_filter)
  cat(sprintf("  After singlets: %s cells\n", format(nrow(current_data), big.mark = ",")))
  
  # Gate 3: Live cells
  live_filter <- point.in.polygon(
    exprs(current_data)[, channels$DCM],
    exprs(current_data)[, channels$SSC_A],
    gates$live_cells[, 1], gates$live_cells[, 2]
  ) > 0
  current_data <- Subset(current_data, live_filter)
  cat(sprintf("  After live: %s cells\n", format(nrow(current_data), big.mark = ",")))
  
  # Gate 4: S-phase outliers
  sphase_filter <- point.in.polygon(
    exprs(current_data)[, channels$FxCycle],
    exprs(current_data)[, channels$EdU],
    gates$s_phase_outliers[, 1], gates$s_phase_outliers[, 2]
  ) > 0
  current_data <- Subset(current_data, sphase_filter)
  cat(sprintf("  After S-phase outliers: %s cells\n", format(nrow(current_data), big.mark = ",")))
  
  # Gate 5: FxCycle quantile
  fxcycle_result <- apply_quantile_gate(current_data, gates$fxcycle_quantile)
  current_data <- fxcycle_result$data
  cat(sprintf("  After FxCycle quantile: %s cells\n", format(nrow(current_data), big.mark = ",")))
  
  # Gate 6: Top 50% EdU + FxCycle range
  edu_values <- exprs(current_data)[, channels$EdU]
  edu_threshold <- quantile(edu_values, probs = 0.50, na.rm = TRUE)
  
  fxcycle_values <- exprs(current_data)[, channels$FxCycle]
  fxcycle_bounds <- quantile(fxcycle_values, probs = c(0.01, 0.90), na.rm = TRUE)
  
  edu_fxcycle_filter <- edu_values >= edu_threshold & 
    fxcycle_values >= fxcycle_bounds[1] & 
    fxcycle_values <= fxcycle_bounds[2]
  current_data <- Subset(current_data, edu_fxcycle_filter)
  cat(sprintf("  After EdU top 50%% + FxCycle range: %s cells\n", format(nrow(current_data), big.mark = ",")))
  
  # Now calculate 98th percentile of HA from gated control
  ha_values <- exprs(current_data)[, channels$HA]
  threshold <- quantile(ha_values, probs = 0.98, na.rm = TRUE)
  
  cat(sprintf("  HA threshold (98th percentile): %s\n", format(round(threshold, 0), big.mark = ",")))
  
  return(list(
    threshold = unname(threshold),
    gated_control_data = current_data,
    n_cells_for_threshold = length(ha_values)
  ))
}



# Helper Functions ----

# Apply sequential gates (used for testing visualizations)
apply_sequential_gates <- function(fcs_data, up_to_gate, gates = GATES, channels = CHANNELS) {
  # Gate 1: Debris
  if(up_to_gate >= 1) {
    debris_filter <- point.in.polygon(
      exprs(fcs_data)[, channels$FSC_A],
      exprs(fcs_data)[, channels$SSC_A],
      gates$debris[, 1], gates$debris[, 2]
    ) > 0
    fcs_data <- Subset(fcs_data, debris_filter)
  }
  
  # Gate 2: Singlets
  if(up_to_gate >= 2) {
    singlet_filter <- point.in.polygon(
      exprs(fcs_data)[, channels$FSC_A],
      exprs(fcs_data)[, channels$FSC_H],
      gates$singlet[, 1], gates$singlet[, 2]
    ) > 0
    fcs_data <- Subset(fcs_data, singlet_filter)
  }
  
  # Gate 3: Live cells
  if(up_to_gate >= 3) {
    live_filter <- point.in.polygon(
      exprs(fcs_data)[, channels$DCM],
      exprs(fcs_data)[, channels$SSC_A],
      gates$live_cells[, 1], gates$live_cells[, 2]
    ) > 0
    fcs_data <- Subset(fcs_data, live_filter)
  }
  
  # Gate 4: S-phase outliers
  if(up_to_gate >= 4) {
    sphase_filter <- point.in.polygon(
      exprs(fcs_data)[, channels$FxCycle],
      exprs(fcs_data)[, channels$EdU],
      gates$s_phase_outliers[, 1], gates$s_phase_outliers[, 2]
    ) > 0
    fcs_data <- Subset(fcs_data, sphase_filter)
  }
  
  return(fcs_data)
}

# Format axis labels (M for millions, K for thousands)
format_axis_labels <- function(x) {
  ifelse(x == 0, "0", 
         ifelse(x >= 1e6, sprintf("%.0fM", x/1e6), sprintf("%.0fK", x/1e3)))
}

# Create density-colored plot data
get_density_colors <- function(x, y, log_x = FALSE, log_y = FALSE) {
  x_plot <- if(log_x) log10(x + 1) else x
  y_plot <- if(log_y) log10(y + 1) else y
  
  densCols(x_plot, y_plot, colramp = colorRampPalette(c("blue", "cyan", "yellow", "red")))
}



# Paths ####
master_path <- "Experiments/"
experiment_folders <- list.dirs(master_path, recursive = FALSE, full.names = TRUE)

# Output folder
OUTPUT_FOLDER <- "Results/"
dir.create(OUTPUT_FOLDER, showWarnings = FALSE, recursive = TRUE)

cat("Found", length(experiment_folders), "experiments\n")
print(experiment_folders)




