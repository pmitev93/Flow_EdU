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

cat("Found", length(experiment_folders), "experiments\n")
print(experiment_folders)

# Channel Mapping ####
CHANNELS <- list(
  FSC_A = "FSC-A",
  FSC_H = "FSC-H",
  SSC_A = "SSC-A",
  EdU = "FL1-A",
  HA = "FL3-A",
  DCM = "FL5-A",
  FxCycle = "FL6-A",
  SAMD9 = "FL10-A"
)

# ==============================================================================
# GATE LOADING
# ==============================================================================

# Specify which gate strategy to use (default)
GATE_STRATEGY_FILE <- "gate_definitions/gates_gdef.r"

# Load the gate definitions
if(!file.exists(GATE_STRATEGY_FILE)) {
  warning(sprintf("Gate definition file not found: %s", GATE_STRATEGY_FILE))
  stop("Please create gate definition file first")
}

cat(sprintf("Loading gate strategy from: %s\n", GATE_STRATEGY_FILE))
source(GATE_STRATEGY_FILE)
if(exists("GATE_STRATEGY")) {
  cat(sprintf("Loaded gate strategy: %s (%s)\n", GATE_STRATEGY$id, GATE_STRATEGY$name))
}

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

# Find control sample in an experiment's metadata
find_control_sample <- function(metadata, pattern = "Control") {
  # First, try to find files ending with "Control.fcs" (case-insensitive)
  control_pattern <- "Control\\.fcs$"
  matches <- grep(control_pattern, metadata$full_filename, ignore.case = TRUE)
  
  if(length(matches) > 0) {
    if(length(matches) > 1) {
      warning(sprintf("Multiple control files found. Using first one: %s", 
                      metadata$full_filename[matches[1]]))
    }
    cat(sprintf("Using control sample: %s (%s)\n", 
                metadata$well[matches[1]], 
                metadata$sample_name[matches[1]]))
    return(matches[1])
  }
  
  # Fallback: try the old pattern matching
  matches <- grep("Empty_Vector_Dox-", metadata$sample_name, ignore.case = TRUE)
  
  if(length(matches) == 0) {
    warning("No control sample found (no file ending in 'Control.fcs' or matching 'Empty_Vector_Dox-')")
    return(NULL)
  }
  
  if(length(matches) > 1) {
    warning(sprintf("Multiple control samples found. Using first one: %s", 
                    metadata$sample_name[matches[1]]))
  }
  
  cat(sprintf("Using control sample: %s (%s)\n", 
              metadata$well[matches[1]], 
              metadata$sample_name[matches[1]]))
  
  return(matches[1])
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
  metadata$well <- str_extract(metadata$full_filename, "[A-H][0-9]+")
  
  # Extract sample name (everything after first underscore, before .fcs)
  metadata$sample_name <- str_split_fixed(metadata$full_filename, "_", 2)[, 2]
  metadata$sample_name <- gsub("\\.fcs$", "", metadata$sample_name, ignore.case = TRUE)
  
  metadata$experiment <- experiment_name
  
  # Load FCS files into flowSet
  fs <- read.flowSet(fcs_files, transformation = FALSE, truncate_max_range = FALSE)
  sampleNames(fs) <- metadata$full_filename
  
  # Apply compensation if available
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

# Apply all gates sequentially to a single sample
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
    gates$debris[, 1],
    gates$debris[, 2]
  ) > 0
  current_data <- Subset(current_data, debris_filter)
  results$cell_counts$after_debris <- nrow(current_data)
  cat(sprintf(" %s cells\n", format(nrow(current_data), big.mark = ",")))
  
  # Gate 2: Singlets (FSC-A vs FSC-H)
  cat("  Gate 2: Singlets...")
  singlet_filter <- point.in.polygon(
    exprs(current_data)[, channels$FSC_A],
    exprs(current_data)[, channels$FSC_H],
    gates$singlet[, 1],
    gates$singlet[, 2]
  ) > 0
  current_data <- Subset(current_data, singlet_filter)
  results$cell_counts$after_singlets <- nrow(current_data)
  cat(sprintf(" %s cells\n", format(nrow(current_data), big.mark = ",")))
  
  # Gate 3: Live cells (DCM-A vs SSC-A)
  cat("  Gate 3: Live cells...")
  live_filter <- point.in.polygon(
    exprs(current_data)[, channels$DCM],
    exprs(current_data)[, channels$SSC_A],
    gates$live_cells[, 1],
    gates$live_cells[, 2]
  ) > 0
  current_data <- Subset(current_data, live_filter)
  results$cell_counts$after_live <- nrow(current_data)
  cat(sprintf(" %s cells\n", format(nrow(current_data), big.mark = ",")))
  
  # Gate 4: S-phase outlier removal (FxCycle vs EdU)
  cat("  Gate 4: S-phase outliers...")
  sphase_filter <- point.in.polygon(
    exprs(current_data)[, channels$FxCycle],
    exprs(current_data)[, channels$EdU],
    gates$s_phase_outliers[, 1],
    gates$s_phase_outliers[, 2]
  ) > 0
  current_data <- Subset(current_data, sphase_filter)
  results$cell_counts$after_sphase_outliers <- nrow(current_data)
  cat(sprintf(" %s cells\n", format(nrow(current_data), big.mark = ",")))
  
  # Gate 5: FxCycle quantile (1%-90%)
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
  if(nrow(current_data) > 10) {
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

  # Calculate HA threshold using percentile from gates parameter
  ha_gate <- gates$ha_positive
  ha_percentile <- ha_gate$prob
  ha_values <- exprs(current_data)[, channels$HA]
  threshold <- quantile(ha_values, probs = ha_percentile, na.rm = TRUE)

  cat(sprintf("  HA threshold (%.0fth percentile): %s\n",
              ha_percentile * 100,
              format(round(threshold, 0), big.mark = ",")))

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

# ==============================================================================
# PAIRED CONTROL FUNCTIONS (for quadrant gating)
# ==============================================================================

# Find paired control for a sample (matches everything before "Dox")
find_paired_control <- function(sample_name, metadata) {
  # Extract the base name (everything before "Dox")
  base_pattern <- gsub("Dox[+-].*", "", sample_name)

  # Find corresponding Dox- control
  control_pattern <- paste0("^", base_pattern, "Dox-")
  matches <- grep(control_pattern, metadata$sample_name, ignore.case = FALSE)

  if(length(matches) == 0) {
    warning(sprintf("No paired Dox- control found for: %s", sample_name))
    return(NULL)
  }

  if(length(matches) > 1) {
    warning(sprintf("Multiple paired controls found for %s. Using first one.", sample_name))
  }

  return(matches[1])
}

# Calculate quadrant thresholds and ratio from paired control
calculate_quadrant_from_paired_control <- function(control_fcs, test_fcs,
                                                     control_name, test_name,
                                                     gates = GATES, channels = CHANNELS) {

  cat(sprintf("\nQuadrant analysis: %s (control) vs %s (test)\n", control_name, test_name))

  # Apply gates 1-6 to control
  control_gated <- control_fcs

  # Gate 1: Debris
  debris_filter <- point.in.polygon(
    exprs(control_gated)[, channels$FSC_A],
    exprs(control_gated)[, channels$SSC_A],
    gates$debris[, 1], gates$debris[, 2]
  ) > 0
  control_gated <- Subset(control_gated, debris_filter)

  # Gate 2: Singlets
  singlet_filter <- point.in.polygon(
    exprs(control_gated)[, channels$FSC_A],
    exprs(control_gated)[, channels$FSC_H],
    gates$singlet[, 1], gates$singlet[, 2]
  ) > 0
  control_gated <- Subset(control_gated, singlet_filter)

  # Gate 3: Live cells
  live_filter <- point.in.polygon(
    exprs(control_gated)[, channels$DCM],
    exprs(control_gated)[, channels$SSC_A],
    gates$live_cells[, 1], gates$live_cells[, 2]
  ) > 0
  control_gated <- Subset(control_gated, live_filter)

  # Gate 4: S-phase outliers (keep cells inside)
  outlier_filter <- point.in.polygon(
    exprs(control_gated)[, channels$FxCycle],
    exprs(control_gated)[, channels$EdU],
    gates$s_phase_outliers[, 1], gates$s_phase_outliers[, 2]
  ) > 0
  control_gated <- Subset(control_gated, outlier_filter)

  # Gate 5: FxCycle quantile
  fxcycle_values <- exprs(control_gated)[, channels$FxCycle]
  fxcycle_limits <- quantile(fxcycle_values, probs = gates$fxcycle_quantile$probs, na.rm = TRUE)
  fxcycle_filter <- fxcycle_values >= fxcycle_limits[1] & fxcycle_values <= fxcycle_limits[2]
  control_gated <- Subset(control_gated, fxcycle_filter)

  # Gate 6: EdU + FxCycle S-phase
  edu_values_g6 <- exprs(control_gated)[, channels$EdU]
  edu_threshold_g6 <- quantile(edu_values_g6, probs = gates$edu_fxcycle_sphase$edu_prob, na.rm = TRUE)
  fxcycle_values_g6 <- exprs(control_gated)[, channels$FxCycle]
  fxcycle_bounds_g6 <- quantile(fxcycle_values_g6, probs = gates$edu_fxcycle_sphase$fxcycle_probs, na.rm = TRUE)
  edu_fxcycle_filter <- edu_values_g6 >= edu_threshold_g6 &
    fxcycle_values_g6 >= fxcycle_bounds_g6[1] &
    fxcycle_values_g6 <= fxcycle_bounds_g6[2]
  control_gated <- Subset(control_gated, edu_fxcycle_filter)

  # Calculate thresholds from control
  ha_threshold <- quantile(exprs(control_gated)[, channels$HA],
                           probs = gates$quadrant$ha_prob, na.rm = TRUE)
  edu_threshold <- quantile(exprs(control_gated)[, channels$EdU],
                            probs = gates$quadrant$edu_prob, na.rm = TRUE)

  cat(sprintf("  Control thresholds: HA = %s, EdU = %s\n",
              format(ha_threshold, big.mark = ","),
              format(edu_threshold, big.mark = ",")))

  # Apply gates 1-6 to test sample
  test_gated <- test_fcs

  # Gate 1: Debris
  debris_filter <- point.in.polygon(
    exprs(test_gated)[, channels$FSC_A],
    exprs(test_gated)[, channels$SSC_A],
    gates$debris[, 1], gates$debris[, 2]
  ) > 0
  test_gated <- Subset(test_gated, debris_filter)

  # Gate 2: Singlets
  singlet_filter <- point.in.polygon(
    exprs(test_gated)[, channels$FSC_A],
    exprs(test_gated)[, channels$FSC_H],
    gates$singlet[, 1], gates$singlet[, 2]
  ) > 0
  test_gated <- Subset(test_gated, singlet_filter)

  # Gate 3: Live cells
  live_filter <- point.in.polygon(
    exprs(test_gated)[, channels$DCM],
    exprs(test_gated)[, channels$SSC_A],
    gates$live_cells[, 1], gates$live_cells[, 2]
  ) > 0
  test_gated <- Subset(test_gated, live_filter)

  # Gate 4: S-phase outliers
  outlier_filter <- point.in.polygon(
    exprs(test_gated)[, channels$FxCycle],
    exprs(test_gated)[, channels$EdU],
    gates$s_phase_outliers[, 1], gates$s_phase_outliers[, 2]
  ) > 0
  test_gated <- Subset(test_gated, outlier_filter)

  # Gate 5: FxCycle quantile
  fxcycle_values <- exprs(test_gated)[, channels$FxCycle]
  fxcycle_limits <- quantile(fxcycle_values, probs = gates$fxcycle_quantile$probs, na.rm = TRUE)
  fxcycle_filter <- fxcycle_values >= fxcycle_limits[1] & fxcycle_values <= fxcycle_limits[2]
  test_gated <- Subset(test_gated, fxcycle_filter)

  # Gate 6: EdU + FxCycle S-phase
  edu_values_g6 <- exprs(test_gated)[, channels$EdU]
  edu_threshold_g6 <- quantile(edu_values_g6, probs = gates$edu_fxcycle_sphase$edu_prob, na.rm = TRUE)
  fxcycle_values_g6 <- exprs(test_gated)[, channels$FxCycle]
  fxcycle_bounds_g6 <- quantile(fxcycle_values_g6, probs = gates$edu_fxcycle_sphase$fxcycle_probs, na.rm = TRUE)
  edu_fxcycle_filter <- edu_values_g6 >= edu_threshold_g6 &
    fxcycle_values_g6 >= fxcycle_bounds_g6[1] &
    fxcycle_values_g6 <= fxcycle_bounds_g6[2]
  test_gated <- Subset(test_gated, edu_fxcycle_filter)

  # Get final HA and EdU values from test sample
  ha_values <- exprs(test_gated)[, channels$HA]
  edu_values <- exprs(test_gated)[, channels$EdU]

  # Calculate quadrant populations
  total_cells <- length(ha_values)

  q1_ha_neg_edu_low <- sum(ha_values < ha_threshold & edu_values < edu_threshold)
  q2_ha_pos_edu_low <- sum(ha_values >= ha_threshold & edu_values < edu_threshold)
  q3_ha_neg_edu_high <- sum(ha_values < ha_threshold & edu_values >= edu_threshold)
  q4_ha_pos_edu_high <- sum(ha_values >= ha_threshold & edu_values >= edu_threshold)

  # Calculate percentages
  q1_pct <- (q1_ha_neg_edu_low / total_cells) * 100
  q2_pct <- (q2_ha_pos_edu_low / total_cells) * 100
  q3_pct <- (q3_ha_neg_edu_high / total_cells) * 100
  q4_pct <- (q4_ha_pos_edu_high / total_cells) * 100

  # Calculate ratio: HA+/EdU-low / (HA+/EdU-low + HA+/EdU-high)
  ratio <- if((q2_pct + q4_pct) > 0) {
    q2_pct / (q2_pct + q4_pct)
  } else {
    NA_real_
  }

  cat(sprintf("  Quadrants: Q1=%.2f%%, Q2=%.2f%%, Q3=%.2f%%, Q4=%.2f%%\n",
              q1_pct, q2_pct, q3_pct, q4_pct))
  cat(sprintf("  Ratio (HA+/EdU-low ratio) = %.3f\n", ratio))

  return(list(
    ha_threshold = ha_threshold,
    edu_threshold = edu_threshold,
    ratio = ratio,
    total_cells = total_cells,
    q1_count = q1_ha_neg_edu_low,
    q2_count = q2_ha_pos_edu_low,
    q3_count = q3_ha_neg_edu_high,
    q4_count = q4_ha_pos_edu_high,
    q1_pct = q1_pct,
    q2_pct = q2_pct,
    q3_pct = q3_pct,
    q4_pct = q4_pct,
    control_name = control_name,
    test_name = test_name
  ))
}

# Wrapper function: Extract correlations with optional quadrant analysis
extract_correlations_with_quadrants <- function(experiment, ha_threshold = NULL,
                                                 gates = GATES, channels = CHANNELS,
                                                 use_quadrant = FALSE) {

  # Load the standard plotting function
  results_df <- extract_correlations(experiment, ha_threshold, gates, channels)

  # Add Ratio column (initialize as NA)
  results_df$Ratio <- NA_real_

  # If not using quadrant strategy, return standard results
  if(!use_quadrant || is.null(gates$quadrant)) {
    return(results_df)
  }

  cat("\n=== QUADRANT ANALYSIS MODE ===\n")

  # For quadrant strategy, calculate ratios for Dox+ samples
  for(i in seq_along(experiment$flowset)) {
    sample_name <- experiment$metadata$sample_name[i]

    # Skip Dox- samples (they're controls)
    if(grepl("Dox-", sample_name, ignore.case = TRUE)) {
      next
    }

    # Find paired control
    control_idx <- find_paired_control(sample_name, experiment$metadata)

    if(is.null(control_idx)) {
      cat(sprintf("WARNING: No paired control for %s - skipping quadrant analysis\n", sample_name))
      next
    }

    # Get control and test FCS data
    control_fcs <- experiment$flowset[[control_idx]]
    test_fcs <- experiment$flowset[[i]]
    control_name <- experiment$metadata$sample_name[control_idx]

    # Calculate quadrant metrics
    tryCatch({
      quadrant_result <- calculate_quadrant_from_paired_control(
        control_fcs, test_fcs,
        control_name, sample_name,
        gates, channels
      )

      # Update the ratio in results_df
      well_match <- which(results_df$Sample == sample_name)
      if(length(well_match) > 0) {
        results_df$Ratio[well_match[1]] <- quadrant_result$ratio
      }

    }, error = function(e) {
      cat(sprintf("ERROR in quadrant analysis for %s: %s\n", sample_name, e$message))
    })
  }

  return(results_df)
}

# Paths ####
master_path <- "Experiments/"
experiment_folders <- list.dirs(master_path, recursive = FALSE, full.names = TRUE)

# Output folder
OUTPUT_FOLDER <- "Results/"
dir.create(OUTPUT_FOLDER, showWarnings = FALSE, recursive = TRUE)

cat("Found", length(experiment_folders), "experiments\n")
print(experiment_folders)