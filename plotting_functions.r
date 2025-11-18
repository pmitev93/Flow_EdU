# ==============================================================================
# FLOW CYTOMETRY GUI - PLOTTING FUNCTIONS
# ==============================================================================

## Quick scan of experiment metadata without loading full FCS data ----
quick_scan_experiment <- function(experiment_path) {
  experiment_name <- basename(experiment_path)
  
  # Get all FCS files
  fcs_files <- list.files(experiment_path, pattern = "\\.fcs$", 
                          full.names = TRUE, ignore.case = TRUE)
  
  if(length(fcs_files) == 0) return(NULL)
  
  # Parse filenames only (no FCS data loading)
  metadata <- data.frame(
    full_filename = basename(fcs_files),
    stringsAsFactors = FALSE
  )
  
  # Extract well and sample info
  metadata$well <- str_extract(metadata$full_filename, "[A-H][0-9]+")
  metadata$sample_name <- str_split_fixed(metadata$full_filename, "_", 2)[, 2]
  metadata$sample_name <- gsub("\\.fcs$", "", metadata$sample_name, ignore.case = TRUE)
  
  # Parse sample components
  metadata$cell_line <- str_extract(metadata$sample_name, "^\\d{3}")
  metadata$gene <- str_extract(metadata$sample_name, "SAMD9L?")
  metadata$gene[is.na(metadata$gene)] <- ""
  
  # Extract mutation
  after_gene <- str_replace(metadata$sample_name, "^.*?SAMD9L?_", "")
  metadata$mutation <- str_replace(after_gene, "_Dox[+-]$", "")
  metadata$mutation[is.na(metadata$mutation) | metadata$mutation == ""] <- ""
  
  # Add placeholders
  metadata$experiment <- experiment_name
  metadata$correlation <- "Not analyzed"
  metadata$ha_pos_pct <- NA_real_
  metadata$n_cells <- "Not analyzed"
  metadata$strength_ratio <- NA_real_
  metadata$notes <- ""

  # Reorder columns
  metadata <- metadata[, c("experiment", "well", "sample_name", "cell_line",
                           "gene", "mutation", "correlation", "ha_pos_pct", "n_cells", "strength_ratio", "notes")]

  # Capitalize column names
  colnames(metadata) <- c("Experiment", "Well", "Sample", "Cell_line",
                          "Gene", "Mutation", "Correlation", "HA_Pos_Pct", "N_cells", "Strength_Ratio", "Notes")
  
  return(metadata)
}

## Gate 1 ----
# Single sample visualization 
plot_debris_gate_single <- function(fcs_data, sample_name, gates = GATES, channels = CHANNELS, show_sample_name = TRUE) {
  # Plot settings
  x <- exprs(fcs_data)[, channels$FSC_A]
  y <- exprs(fcs_data)[, channels$SSC_A]

  # Calculate 2D density
  dens <- densCols(x, y, colramp = colorRampPalette(c("blue", "cyan", "yellow", "red")))

  # Create title
  plot_title <- if(show_sample_name) {
    sprintf("Gate 1: Debris Removal\n%s", sample_name)
  } else {
    "Gate 1: Debris Removal"
  }

  plot(x, y,
       pch = ".",
       col = dens,
       xlab = "FSC-A",
       ylab = "SSC-A",
       main = plot_title,
       xlim = c(0, 20e6),
       ylim = c(0, 20e6),
       mgp = c(3, 0.5, 0),
       xaxs = "i",
       yaxs = "i",
       xaxt = "n",
       yaxt = "n",)
  
  # Add these lines right here:
  format_axis <- function(x) {
    ifelse(x == 0, "0", 
           ifelse(x >= 1e6, sprintf("%.0fM", x/1e6), sprintf("%.0fK", x/1e3)))
  }
  
  axis(1, at = seq(0, 15e6, 5e6), labels = format_axis(seq(0, 15e6, 5e6)), mgp = c(3, 0.5, 0))
  axis(2, at = seq(0, 15e6, 5e6), labels = format_axis(seq(0, 15e6, 5e6)), mgp = c(3, 0.5, 0))
  
  # Add gate polygon
  polygon(gates$debris[, 1], gates$debris[, 2], 
          border = "black", lwd = 1)
  
  # Calculate stats
  total_cells <- nrow(fcs_data)
  inside_gate <- sum(point.in.polygon(
    exprs(fcs_data)[, channels$FSC_A],
    exprs(fcs_data)[, channels$SSC_A],
    gates$debris[, 1],
    gates$debris[, 2]
  ) > 0)
  
  # Add legend
  legend("bottomright", 
         legend = c(sprintf("Total: %s", format(total_cells, big.mark = ",")),
                    sprintf("Inside gate: %s (%.1f%%)", 
                            format(inside_gate, big.mark = ","),
                            100 * inside_gate / total_cells)),
         bty = "n",
         cex = 1.1)
}

# Overview: All samples on one plot 
plot_debris_gate_overview <- function(experiment, gates = GATES, channels = CHANNELS) {
  n_samples <- length(experiment$flowset)
  
  # Calculate grid dimensions
  n_cols <- ceiling(sqrt(n_samples))
  n_rows <- ceiling(n_samples / n_cols)
  
  par(mfrow = c(n_rows, n_cols), mar = c(2, 2, 1.5, 0.5), oma = c(0, 0, 0, 0), mgp = c(3, 0.3, 0))
  
  # Custom axis formatting function
  format_axis <- function(x) {
    ifelse(x == 0, "0", 
           ifelse(x >= 1e6, sprintf("%.0fM", x/1e6), sprintf("%.0fK", x/1e3)))
  }
  
  for(i in seq_along(experiment$flowset)) {
    fcs_data <- experiment$flowset[[i]]
    sample_name <- experiment$metadata$sample_name[i]
    
    # Simplified plot for overview
    x <- exprs(fcs_data)[, channels$FSC_A]
    y <- exprs(fcs_data)[, channels$SSC_A]
    
    # Calculate 2D density
    dens <- densCols(x, y, colramp = colorRampPalette(c("blue", "cyan", "yellow", "red")))
    
    plot(exprs(fcs_data)[, channels$FSC_A],
         exprs(fcs_data)[, channels$SSC_A],
         pch = ".",
         col = dens,
         xlab = "",
         ylab = "",
         main = sample_name,
         xlim = c(0, 20e6),
         ylim = c(0, 20e6),
         cex.main = 1.2,
         cex.axis = 0.7,
         xaxt = "n",
         yaxt = "n",
         xaxs = "i",
         yaxs = "i")
    
    # Add custom axes
    format_axis <- function(x) {
      ifelse(x == 0, "0", 
             ifelse(x >= 1e6, sprintf("%.0fM", x/1e6), sprintf("%.0fK", x/1e3)))
    }
    
    axis(1, at = seq(0, 15e6, 5e6), labels = format_axis(seq(0, 15e6, 5e6)), mgp = c(3, 0.5, 0))
    axis(2, at = seq(0, 15e6, 5e6), labels = format_axis(seq(0, 15e6, 5e6)), mgp = c(3, 0.5, 0))
    
    polygon(gates$debris[, 1], gates$debris[, 2], 
            border = "black", lwd = 1)
    
    # Add % inside gate
    inside_gate <- sum(point.in.polygon(
      exprs(fcs_data)[, channels$FSC_A],
      exprs(fcs_data)[, channels$SSC_A],
      gates$debris[, 1],
      gates$debris[, 2]
    ) > 0)
    
    pct <- 100 * inside_gate / nrow(fcs_data)
    # Color code based on threshold
    text_col <- ifelse(pct < 75, "red", "black")
    text(17e6, 1e6, sprintf("%.1f%%", pct), col = text_col, cex = 0.8, font = 2, pos = 2)
    
    # Add BIG warning if below threshold
    if(pct < 75) {
      # Red border around plot
      box(col = "red", lwd = 4)
      # Big warning text in center
      text(10e6, 10e6, "LOW\nCELL COUNT", col = "red", cex = 1.5, font = 2)
    }
  }
  
  par(mfrow = c(1, 1), mar = c(5, 4, 4, 2), oma = c(0, 0, 0, 0))  # Reset
}

## Gate 2: Singlets (FSC-A vs FSC-H) ####

plot_singlet_gate_single <- function(fcs_data, sample_name, gates = GATES, channels = CHANNELS, show_sample_name = TRUE) {
  # Apply Gate 1 first (debris removal)
  debris_filter <- point.in.polygon(
    exprs(fcs_data)[, channels$FSC_A],
    exprs(fcs_data)[, channels$SSC_A],
    gates$debris[, 1],
    gates$debris[, 2]
  ) > 0
  fcs_data <- Subset(fcs_data, debris_filter)

  # Plot settings
  x <- exprs(fcs_data)[, channels$FSC_A]
  y <- exprs(fcs_data)[, channels$FSC_H]
  dens <- densCols(x, y, colramp = colorRampPalette(c("blue", "cyan", "yellow", "red")))

  # Create title
  plot_title <- if(show_sample_name) {
    sprintf("Gate 2: Singlets\n%s", sample_name)
  } else {
    "Gate 2: Singlets"
  }

  plot(x, y,
       pch = 16,
       cex = 0.3,
       col = dens,
       xlab = "FSC-A",
       ylab = "FSC-H",
       main = plot_title,
       xlim = c(0, 15e6),
       ylim = c(0, 4e6),
       xaxt = "n",
       yaxt = "n",
       mgp = c(3, 0.5, 0),
       xaxs = "i",
       yaxs = "i")
  
  # Custom axis labels
  format_axis <- function(x) {
    ifelse(x == 0, "0", 
           ifelse(x >= 1e6, sprintf("%.0fM", x/1e6), sprintf("%.0fK", x/1e3)))
  }
  
  # Add these two lines here:
  axis(1, at = seq(0, 15e6, 5e6), labels = c("0M", "5M", "10M", "15M"), mgp = c(3, 0.5, 0))
  axis(2, at = seq(0, 4e6, 1e6), labels = c("0M", "1M", "2M", "3M", "4M"), mgp = c(3, 0.5, 0))
  
  # Add gate polygon
  polygon(gates$singlet[, 1], gates$singlet[, 2], 
          border = "black", lwd = 1.5)
  
  # Calculate stats
  total_cells <- nrow(fcs_data)
  inside_gate <- sum(point.in.polygon(
    exprs(fcs_data)[, channels$FSC_A],
    exprs(fcs_data)[, channels$FSC_H],
    gates$singlet[, 1],
    gates$singlet[, 2]
  ) > 0)
  
  # Add legend
  legend("bottomright", 
         legend = c(sprintf("Total: %s", format(total_cells, big.mark = ",")),
                    sprintf("Inside gate: %s (%.1f%%)", 
                            format(inside_gate, big.mark = ","),
                            100 * inside_gate / total_cells)),
         bty = "n",
         cex = 1.1)
}

plot_singlet_gate_overview <- function(experiment, gates = GATES, channels = CHANNELS) {
  # Apply Gate 1 first (debris removal)
  debris_filter <- point.in.polygon(
    exprs(fcs_data)[, channels$FSC_A],
    exprs(fcs_data)[, channels$SSC_A],
    gates$debris[, 1],
    gates$debris[, 2]
  ) > 0
  fcs_data <- Subset(fcs_data, debris_filter)
  
  n_samples <- length(experiment$flowset)
  n_cols <- ceiling(sqrt(n_samples))
  n_rows <- ceiling(n_samples / n_cols)
  
  par(mfrow = c(n_rows, n_cols), mar = c(2, 2, 1.5, 0.5), oma = c(0, 0, 0, 0), mgp = c(3, 0.3, 0))
  
  format_axis <- function(x) {
    ifelse(x == 0, "0", 
           ifelse(x >= 1e6, sprintf("%.0fM", x/1e6), sprintf("%.0fK", x/1e3)))
  }
  
  for(i in seq_along(experiment$flowset)) {
    fcs_data <- experiment$flowset[[i]]
    sample_name <- experiment$metadata$sample_name[i]
    
    # Add these lines here to apply Gate 1 first:
    debris_filter <- point.in.polygon(
      exprs(fcs_data)[, channels$FSC_A],
      exprs(fcs_data)[, channels$SSC_A],
      gates$debris[, 1],
      gates$debris[, 2]
    ) > 0
    fcs_data <- Subset(fcs_data, debris_filter)
    
    x <- exprs(fcs_data)[, channels$FSC_A]
    y <- exprs(fcs_data)[, channels$FSC_H]
    
    dens <- densCols(x, y, colramp = colorRampPalette(c("blue", "cyan", "yellow", "red")))
    
    plot(x, y,
         pch = ".",
         col = dens,
         xlab = "",
         ylab = "",
         main = sample_name,
         xlim = c(0, 15e6),
         ylim = c(0, 4e6),
         cex.main = 1.2,
         cex.axis = 0.7,
         xaxt = "n",
         yaxt = "n",
         mgp = c(3, 0.5, 0),
         xaxs = "i",
         yaxs = "i")
    
    axis(1, at = seq(0, 20e6, 5e6), labels = format_axis(seq(0, 20e6, 5e6)), cex.axis = 0.5)
    axis(2, at = seq(0, 4e6, 1e6), labels = format_axis(seq(0, 4e6, 1e6)), cex.axis = 0.5)
    
    polygon(gates$singlet[, 1], gates$singlet[, 2], 
            border = "black", lwd = 1)
    
    inside_gate <- sum(point.in.polygon(
      exprs(fcs_data)[, channels$FSC_A],
      exprs(fcs_data)[, channels$FSC_H],
      gates$singlet[, 1],
      gates$singlet[, 2]
    ) > 0)
    
    pct <- 100 * inside_gate / nrow(fcs_data)
    text_col <- ifelse(pct < 50, "red", "black")
    text(13e6, 0.5e6, sprintf("%.1f%%", pct), col = text_col, cex = 0.8, font = 2, pos = 2)
    
    if(pct < 70) {
      box(col = "red", lwd = 4)
      text(10e6, 2e6, "LOW\nCELL COUNT", col = "red", cex = 1.5, font = 2)
    }
  }
  
  par(mfrow = c(1, 1), mar = c(5, 4, 4, 2), oma = c(0, 0, 0, 0))
}


## Gate 3: Live Cells (DCM-A vs SSC-A) ####

plot_live_gate_single <- function(fcs_data, sample_name, gates = GATES, channels = CHANNELS, show_sample_name = TRUE) {
  # Apply Gate 1: Debris removal
  debris_filter <- point.in.polygon(
    exprs(fcs_data)[, channels$FSC_A],
    exprs(fcs_data)[, channels$SSC_A],
    gates$debris[, 1],
    gates$debris[, 2]
  ) > 0
  fcs_data <- Subset(fcs_data, debris_filter)
  
  # Apply Gate 2: Singlets
  singlet_filter <- point.in.polygon(
    exprs(fcs_data)[, channels$FSC_A],
    exprs(fcs_data)[, channels$FSC_H],
    gates$singlet[, 1],
    gates$singlet[, 2]
  ) > 0
  fcs_data <- Subset(fcs_data, singlet_filter)
  
  # Plot settings
  # Plot settings with log-scaled DCM
  x <- exprs(fcs_data)[, channels$DCM]
  y <- exprs(fcs_data)[, channels$SSC_A]

  # Log transform DCM for density calculation and plotting
  x_log <- log10(x + 1)  # +1 to avoid log(0)

  dens <- densCols(x_log, y, colramp = colorRampPalette(c("blue", "cyan", "yellow", "red")))

  # Create title
  plot_title <- if(show_sample_name) {
    sprintf("Gate 3: Live Cells\n%s", sample_name)
  } else {
    "Gate 3: Live Cells"
  }

  plot(x, y,
       pch = 16,
       cex = 0.3,
       col = dens,
       xlab = "DCM-A",
       ylab = "SSC-A",
       main = plot_title,
       xlim = c(100, 1000000),
       ylim = c(0, 15e6),
       xaxt = "n",
       yaxt = "n",
       mgp = c(3, 0.5, 0),
       xaxs = "i",
       yaxs = "i",
       log = "x")  # Log scale x-axis
  
  # Custom axis labels for log scale
  format_axis <- function(x) {
    ifelse(x == 0, "0", 
           ifelse(x >= 1e6, sprintf("%.0fM", x/1e6), sprintf("%.0fK", x/1e3)))
  }
  
  axis(1, at = c(100, 1000, 10000, 100000, 1000000), labels = c("100", "1K", "10K", "100K", "1M"), mgp = c(3, 0.5, 0))
  axis(2, at = seq(0, 15e6, 5e6), labels = format_axis(seq(0, 15e6, 5e6)), mgp = c(3, 0.5, 0))
  
  # Add gate polygon
  polygon(gates$live_cells[, 1], gates$live_cells[, 2], 
          border = "black", lwd = 1.5)
  
  # Calculate stats
  total_cells <- nrow(fcs_data)
  inside_gate <- sum(point.in.polygon(
    exprs(fcs_data)[, channels$DCM],
    exprs(fcs_data)[, channels$SSC_A],
    gates$live_cells[, 1],
    gates$live_cells[, 2]
  ) > 0)
  
  # Add legend
  legend("topright", 
         legend = c(sprintf("Total: %s", format(total_cells, big.mark = ",")),
                    sprintf("Inside gate: %s (%.1f%%)", 
                            format(inside_gate, big.mark = ","),
                            100 * inside_gate / total_cells)),
         bty = "n",
         cex = 1.1)
}

plot_live_gate_overview <- function(experiment, gates = GATES, channels = CHANNELS) {
  n_samples <- length(experiment$flowset)
  n_cols <- ceiling(sqrt(n_samples))
  n_rows <- ceiling(n_samples / n_cols)
  
  par(mfrow = c(n_rows, n_cols), mar = c(2, 2, 1.5, 0.5), oma = c(0, 0, 0, 0), mgp = c(3, 0.3, 0))
  
  format_axis <- function(x) {
    ifelse(x == 0, "0", 
           ifelse(x >= 1e6, sprintf("%.0fM", x/1e6), sprintf("%.0fK", x/1e3)))
  }
  
  for(i in seq_along(experiment$flowset)) {
    fcs_data <- experiment$flowset[[i]]
    sample_name <- experiment$metadata$sample_name[i]
    
    # Apply Gate 1: Debris
    debris_filter <- point.in.polygon(
      exprs(fcs_data)[, channels$FSC_A],
      exprs(fcs_data)[, channels$SSC_A],
      gates$debris[, 1],
      gates$debris[, 2]
    ) > 0
    fcs_data <- Subset(fcs_data, debris_filter)
    
    # Apply Gate 2: Singlets
    singlet_filter <- point.in.polygon(
      exprs(fcs_data)[, channels$FSC_A],
      exprs(fcs_data)[, channels$FSC_H],
      gates$singlet[, 1],
      gates$singlet[, 2]
    ) > 0
    fcs_data <- Subset(fcs_data, singlet_filter)
    
    x <- exprs(fcs_data)[, channels$DCM]
    y <- exprs(fcs_data)[, channels$SSC_A]
    
    # Log transform for density
    x_log <- log10(x + 1)
    
    dens <- densCols(x_log, y, colramp = colorRampPalette(c("blue", "cyan", "yellow", "red")))
    
    plot(x, y,
         pch = ".",
         col = dens,
         xlab = "",
         ylab = "",
         main = sample_name,
         xlim = c(100, 1000000),
         ylim = c(0, 15e6),
         cex.main = 1.2,
         cex.axis = 0.7,
         xaxt = "n",
         yaxt = "n",
         log = "x")
    
    axis(1, at = c(100, 1000, 10000, 100000, 1000000), labels = c("100", "1K", "10K", "100K", "1M"), mgp = c(3, 0.5, 0))
    axis(2, at = seq(0, 15e6, 5e6), labels = format_axis(seq(0, 15e6, 5e6)), cex.axis = 0.5)
    
    polygon(gates$live_cells[, 1], gates$live_cells[, 2], 
            border = "black", lwd = 1)
    
    inside_gate <- sum(point.in.polygon(
      exprs(fcs_data)[, channels$DCM],
      exprs(fcs_data)[, channels$SSC_A],
      gates$live_cells[, 1],
      gates$live_cells[, 2]
    ) > 0)
    
    pct <- 100 * inside_gate / nrow(fcs_data)
    text_col <- ifelse(pct < 50, "red", "black")
    text(850000, 14e6, sprintf("%.1f%%", pct), col = text_col, cex = 0.8, font = 2, pos = 2)
    
    if(pct < 50) {
      box(col = "red", lwd = 4)
      text(1000, 7.5e6, "LOW\nCELL COUNT", col = "red", cex = 1.5, font = 2)
    }
  }
  
  par(mfrow = c(1, 1), mar = c(5, 4, 4, 2), oma = c(0, 0, 0, 0))
}

## Gate 4: S-phase Outlier Removal (FxCycle-A vs EdU-A) ----

plot_sphase_outlier_gate_single <- function(fcs_data, sample_name, gates = GATES, channels = CHANNELS, show_sample_name = TRUE) {
  # Apply Gates 1-3
  debris_filter <- point.in.polygon(
    exprs(fcs_data)[, channels$FSC_A],
    exprs(fcs_data)[, channels$SSC_A],
    gates$debris[, 1],
    gates$debris[, 2]
  ) > 0
  fcs_data <- Subset(fcs_data, debris_filter)
  
  singlet_filter <- point.in.polygon(
    exprs(fcs_data)[, channels$FSC_A],
    exprs(fcs_data)[, channels$FSC_H],
    gates$singlet[, 1],
    gates$singlet[, 2]
  ) > 0
  fcs_data <- Subset(fcs_data, singlet_filter)
  
  live_filter <- point.in.polygon(
    exprs(fcs_data)[, channels$DCM],
    exprs(fcs_data)[, channels$SSC_A],
    gates$live_cells[, 1],
    gates$live_cells[, 2]
  ) > 0
  fcs_data <- Subset(fcs_data, live_filter)
  
  # Plot settings
  # Log transform EdU for density
  
  x <- exprs(fcs_data)[, channels$FxCycle]
  y <- exprs(fcs_data)[, channels$EdU]
  y_log <- log10(y + 1)

  dens <- densCols(x, y_log, colramp = colorRampPalette(c("blue", "cyan", "yellow", "red")))

  # Create title
  plot_title <- if(show_sample_name) {
    sprintf("Gate 4: S-phase Outlier Removal\n%s", sample_name)
  } else {
    "Gate 4: S-phase Outlier Removal"
  }

  plot(x, y,
       pch = 16,
       cex = 0.3,
       col = dens,
       xlab = "FxCycle-A",
       ylab = "EdU-A",
       main = plot_title,
       xlim = c(0, 12e6),
       ylim = c(100, 6e6),
       xaxt = "n",
       yaxt = "n",
       mgp = c(3, 0.5, 0),
       xaxs = "i",
       yaxs = "i",
       log= "y")
  
  # Custom axis labels
  format_axis <- function(x) {
    ifelse(x == 0, "0", 
           ifelse(x >= 1e6, sprintf("%.0fM", x/1e6), sprintf("%.0fK", x/1e3)))
  }
  
  axis(1, at = seq(0, 12e6, 3e6), labels = format_axis(seq(0, 12e6, 3e6)), mgp = c(3, 0.5, 0))
  axis(2, at = c(100, 1000, 10000, 100000, 1000000), labels = c("100", "1K", "10K", "100K", "1M"), mgp = c(3, 0.5, 0))
  
  # Add gate polygon
  polygon(gates$s_phase_outliers[, 1], gates$s_phase_outliers[, 2], 
          border = "black", lwd = 1.5)
  
  # Calculate stats
  total_cells <- nrow(fcs_data)
  inside_gate <- sum(point.in.polygon(
    exprs(fcs_data)[, channels$FxCycle],
    exprs(fcs_data)[, channels$EdU],
    gates$s_phase_outliers[, 1],
    gates$s_phase_outliers[, 2]
  ) > 0)
  
  # Add legend
  legend("bottomright", 
         legend = c(sprintf("Total: %s", format(total_cells, big.mark = ",")),
                    sprintf("Inside gate: %s (%.1f%%)", 
                            format(inside_gate, big.mark = ","),
                            100 * inside_gate / total_cells)),
         bty = "n",
         cex = 1.1)
}

plot_sphase_outlier_gate_overview <- function(experiment, gates = GATES, channels = CHANNELS) {
  n_samples <- length(experiment$flowset)
  n_cols <- ceiling(sqrt(n_samples))
  n_rows <- ceiling(n_samples / n_cols)
  
  par(mfrow = c(n_rows, n_cols), mar = c(2, 2, 1.5, 0.5), oma = c(0, 0, 0, 0), mgp = c(3, 0.3, 0))
  
  format_axis <- function(x) {
    ifelse(x == 0, "0", 
           ifelse(x >= 1e6, sprintf("%.0fM", x/1e6), sprintf("%.0fK", x/1e3)))
  }
  
  for(i in seq_along(experiment$flowset)) {
    fcs_data <- experiment$flowset[[i]]
    sample_name <- experiment$metadata$sample_name[i]
    
    # Apply Gates 1-3
    debris_filter <- point.in.polygon(
      exprs(fcs_data)[, channels$FSC_A],
      exprs(fcs_data)[, channels$SSC_A],
      gates$debris[, 1],
      gates$debris[, 2]
    ) > 0
    fcs_data <- Subset(fcs_data, debris_filter)
    
    singlet_filter <- point.in.polygon(
      exprs(fcs_data)[, channels$FSC_A],
      exprs(fcs_data)[, channels$FSC_H],
      gates$singlet[, 1],
      gates$singlet[, 2]
    ) > 0
    fcs_data <- Subset(fcs_data, singlet_filter)
    
    live_filter <- point.in.polygon(
      exprs(fcs_data)[, channels$DCM],
      exprs(fcs_data)[, channels$SSC_A],
      gates$live_cells[, 1],
      gates$live_cells[, 2]
    ) > 0
    fcs_data <- Subset(fcs_data, live_filter)
    
    x <- exprs(fcs_data)[, channels$FxCycle]
    y <- exprs(fcs_data)[, channels$EdU]
    y_log <- log10(y + 1)
    
    dens <- densCols(x, y_log, colramp = colorRampPalette(c("blue", "cyan", "yellow", "red")))
    
    plot(x, y,
         pch = ".",
         col = dens,
         xlab = "",
         ylab = "",
         main = sample_name,
         xlim = c(0, 12e6),
         ylim = c(100, 6e6),
         cex.main = 1.2,
         cex.axis = 0.7,
         xaxt = "n",
         yaxt = "n",
         xaxs = "i",
         yaxs = "i",
         log= "y")
    
    axis(1, at = seq(0, 12e6, 3e6), labels = format_axis(seq(0, 12e6, 3e6)), cex.axis = 0.5)
    axis(2, at = c(100, 1000, 10000, 100000, 1000000), labels = c("100", "1K", "10K", "100K", "1M"), mgp = c(3, 0.5, 0))
    
    polygon(gates$s_phase_outliers[, 1], gates$s_phase_outliers[, 2], 
            border = "black", lwd = 1)
    
    inside_gate <- sum(point.in.polygon(
      exprs(fcs_data)[, channels$FxCycle],
      exprs(fcs_data)[, channels$EdU],
      gates$s_phase_outliers[, 1],
      gates$s_phase_outliers[, 2]
    ) > 0)
    
    pct <- 100 * inside_gate / nrow(fcs_data)
    text_col <- ifelse(pct < 50, "red", "black")
    text(12e6, 3e2, sprintf("%.1f%%", pct), col = text_col, cex = 0.8, font = 2, pos = 2)
    
    if(pct < 70) {
      box(col = "red", lwd = 4)
      text(6e6, 3e5, "LOW\nCELL COUNT", col = "red", cex = 1.5, font = 2)
    }
  }
  
  par(mfrow = c(1, 1), mar = c(5, 4, 4, 2), oma = c(0, 0, 0, 0))
}

## Gate 5: FxCycle Quantile (1%-90%) ----

plot_fxcycle_quantile_gate_single <- function(fcs_data, sample_name, gates = GATES, channels = CHANNELS, show_sample_name = TRUE) {
  # Apply Gates 1-4
  fcs_data <- apply_sequential_gates(fcs_data, up_to_gate = 4, gates = gates, channels = channels)

  # Calculate FxCycle quantile bounds - read from gates
  fxcycle_gate <- gates$fxcycle_quantile
  fxcycle_values <- exprs(fcs_data)[, channels$FxCycle]
  quantile_limits <- quantile(fxcycle_values, probs = fxcycle_gate$probs, na.rm = TRUE)
  lower_bound <- quantile_limits[1]
  upper_bound <- quantile_limits[2]
  
  # Plot settings
  x <- exprs(fcs_data)[, channels$FxCycle]
  y <- exprs(fcs_data)[, channels$EdU]

  dens <- get_density_colors(x, y, log_x = FALSE, log_y = TRUE)

  # Create title
  plot_title <- if(show_sample_name) {
    sprintf("Gate 5: FxCycle Quantile (1%%-90%%)\n%s", sample_name)
  } else {
    "Gate 5: FxCycle Quantile (1%-90%)"
  }

  plot(x, y,
       pch = 16,
       cex = 0.3,
       col = dens,
       xlab = "FxCycle-A",
       ylab = "EdU-A",
       main = plot_title,
       xlim = c(0, 12e6),
       ylim = c(100, 6e6),
       xaxt = "n",
       yaxt = "n",
       mgp = c(3, 0.5, 0),
       xaxs = "i",
       yaxs = "i",
       log = "y")
  
  axis(1, at = seq(0, 12e6, 3e6), labels = format_axis_labels(seq(0, 12e6, 3e6)), mgp = c(3, 0.5, 0))
  axis(2, at = c(100, 1000, 10000, 100000, 1000000), labels = c("100", "1K", "10K", "100K", "1M"), mgp = c(3, 0.5, 0))
  
  # Add vertical lines for quantile bounds
  abline(v = lower_bound, col = "black", lwd = 2, lty = 2)
  abline(v = upper_bound, col = "black", lwd = 2, lty = 2)
  
  # Calculate stats
  total_cells <- nrow(fcs_data)
  inside_gate <- sum(fxcycle_values >= lower_bound & fxcycle_values <= upper_bound)
  
  # Add legend
  legend("bottomright", 
         legend = c(sprintf("Total: %s", format(total_cells, big.mark = ",")),
                    sprintf("Inside gate: %s (%.1f%%)", 
                            format(inside_gate, big.mark = ","),
                            100 * inside_gate / total_cells),
                    sprintf("Bounds: %.1fM - %.1fM", lower_bound/1e6, upper_bound/1e6)),
         bty = "n",
         cex = 0.9)
}

plot_fxcycle_quantile_gate_overview <- function(experiment, gates = GATES, channels = CHANNELS) {
  n_samples <- length(experiment$flowset)
  n_cols <- ceiling(sqrt(n_samples))
  n_rows <- ceiling(n_samples / n_cols)

  par(mfrow = c(n_rows, n_cols), mar = c(2, 2, 1.5, 0.5), oma = c(0, 0, 0, 0), mgp = c(3, 0.3, 0))

  # Read gate parameters
  fxcycle_gate <- gates$fxcycle_quantile

  for(i in seq_along(experiment$flowset)) {
    fcs_data <- experiment$flowset[[i]]
    sample_name <- experiment$metadata$sample_name[i]

    # Apply Gates 1-4
    fcs_data <- apply_sequential_gates(fcs_data, up_to_gate = 4, gates = gates, channels = channels)

    # Calculate quantile bounds - read from gates
    fxcycle_values <- exprs(fcs_data)[, channels$FxCycle]
    quantile_limits <- quantile(fxcycle_values, probs = fxcycle_gate$probs, na.rm = TRUE)
    lower_bound <- quantile_limits[1]
    upper_bound <- quantile_limits[2]
    
    x <- exprs(fcs_data)[, channels$FxCycle]
    y <- exprs(fcs_data)[, channels$EdU]
    
    dens <- get_density_colors(x, y, log_x = FALSE, log_y = TRUE)
    
    plot(x, y,
         pch = ".",
         col = dens,
         xlab = "",
         ylab = "",
         main = sample_name,
         xlim = c(0, 12e6),
         ylim = c(100, 6e6),
         cex.main = 1.2,
         cex.axis = 0.7,
         xaxt = "n",
         yaxt = "n",
         log = "y")
    
    axis(1, at = seq(0, 12e6, 3e6), labels = format_axis_labels(seq(0, 12e6, 3e6)), cex.axis = 0.5)
    axis(2, at = c(100, 1000, 10000, 100000, 1000000), labels = c("100", "1K", "10K", "100K", "1M"), cex.axis = 0.5)
    
    # Add vertical lines
    abline(v = lower_bound, col = "black", lwd = 1, lty = 2)
    abline(v = upper_bound, col = "black", lwd = 1, lty = 2)
    
    inside_gate <- sum(fxcycle_values >= lower_bound & fxcycle_values <= upper_bound)
    pct <- 100 * inside_gate / nrow(fcs_data)
    text_col <- ifelse(pct < 50, "red", "black")
    text(10e6, 200, sprintf("%.1f%%", pct), col = text_col, cex = 0.8, font = 2, pos = 2)
    
    if(pct < 50) {
      box(col = "red", lwd = 4)
      text(6e6, 10000, "LOW\nCELL COUNT", col = "red", cex = 1.5, font = 2)
    }
  }
  
  par(mfrow = c(1, 1), mar = c(5, 4, 4, 2), oma = c(0, 0, 0, 0))
}

## Gate 6: EdU + FxCycle Range (reads percentile from gates) ----

plot_edu_fxcycle_gate_single <- function(fcs_data, sample_name, gates = GATES, channels = CHANNELS, show_sample_name = TRUE) {
  # Apply Gates 1-5
  fcs_data <- apply_sequential_gates(fcs_data, up_to_gate = 4, gates = gates, channels = channels)

  # Apply Gate 5 (FxCycle quantile) - read from gates
  fxcycle_gate <- gates$fxcycle_quantile
  fxcycle_values <- exprs(fcs_data)[, channels$FxCycle]
  fxcycle_limits <- quantile(fxcycle_values, probs = fxcycle_gate$probs, na.rm = TRUE)
  fxcycle_filter <- fxcycle_values >= fxcycle_limits[1] & fxcycle_values <= fxcycle_limits[2]
  fcs_data <- Subset(fcs_data, fxcycle_filter)

  # Calculate EdU threshold and FxCycle bounds - read from gates
  edu_gate <- gates$edu_fxcycle_sphase
  edu_prob <- edu_gate$edu_prob

  edu_values <- exprs(fcs_data)[, channels$EdU]
  edu_threshold <- quantile(edu_values, probs = edu_prob, na.rm = TRUE)

  fxcycle_values <- exprs(fcs_data)[, channels$FxCycle]
  fxcycle_bounds <- quantile(fxcycle_values, probs = edu_gate$fxcycle_probs, na.rm = TRUE)

  # Plot settings
  x <- exprs(fcs_data)[, channels$FxCycle]
  y <- exprs(fcs_data)[, channels$EdU]

  dens <- get_density_colors(x, y, log_x = FALSE, log_y = TRUE)

  # Create title
  plot_title <- if(show_sample_name) {
    sprintf("Gate 6: Top %.0f%% EdU + FxCycle Range\n%s", edu_prob * 100, sample_name)
  } else {
    sprintf("Gate 6: Top %.0f%% EdU + FxCycle Range", edu_prob * 100)
  }

  plot(x, y,
       pch = 16,
       cex = 0.3,
       col = dens,
       xlab = "FxCycle-A",
       ylab = "EdU-A",
       main = plot_title,
       xlim = c(0, 12e6),
       ylim = c(100, 6e6),
       xaxt = "n",
       yaxt = "n",
       mgp = c(3, 0.5, 0),
       xaxs = "i",
       yaxs = "i",
       log = "y")
  
  axis(1, at = seq(0, 12e6, 3e6), labels = format_axis_labels(seq(0, 12e6, 3e6)), mgp = c(3, 0.5, 0))
  axis(2, at = c(100, 1000, 10000, 100000, 1000000), labels = c("100", "1K", "10K", "100K", "1M"), mgp = c(3, 0.5, 0))
  
  # Add gate lines (only in gated region)
  segments(x0 = fxcycle_bounds[1], y0 = edu_threshold, 
           x1 = fxcycle_bounds[2], y1 = edu_threshold, 
           col = "black", lwd = 2, lty = 2)  # Horizontal line between FxCycle bounds
  segments(x0 = fxcycle_bounds[1], y0 = edu_threshold, 
           x1 = fxcycle_bounds[1], y1 = 6e6, 
           col = "black", lwd = 2, lty = 2)  # Left vertical line
  segments(x0 = fxcycle_bounds[2], y0 = edu_threshold, 
           x1 = fxcycle_bounds[2], y1 = 6e6, 
           col = "black", lwd = 2, lty = 2)  # Right vertical line
  
  # Calculate stats
  total_cells <- nrow(fcs_data)
  inside_gate <- sum(edu_values >= edu_threshold & 
                       fxcycle_values >= fxcycle_bounds[1] & 
                       fxcycle_values <= fxcycle_bounds[2])
  
  # Add legend
  legend("bottomright", 
         legend = c(sprintf("Total: %s", format(total_cells, big.mark = ",")),
                    sprintf("Inside gate: %s (%.1f%%)", 
                            format(inside_gate, big.mark = ","),
                            100 * inside_gate / total_cells),
                    sprintf("EdU ≥ %.0fK", edu_threshold/1e3)),
         bty = "n",
         cex = 0.9)
}

plot_edu_fxcycle_gate_overview <- function(experiment, gates = GATES, channels = CHANNELS) {
  n_samples <- length(experiment$flowset)
  n_cols <- ceiling(sqrt(n_samples))
  n_rows <- ceiling(n_samples / n_cols)

  par(mfrow = c(n_rows, n_cols), mar = c(2, 2, 1.5, 0.5), oma = c(0, 0, 0, 0), mgp = c(3, 0.3, 0))

  # Read gate parameters
  fxcycle_gate <- gates$fxcycle_quantile
  edu_gate <- gates$edu_fxcycle_sphase

  for(i in seq_along(experiment$flowset)) {
    fcs_data <- experiment$flowset[[i]]
    sample_name <- experiment$metadata$sample_name[i]

    # Apply Gates 1-5
    fcs_data <- apply_sequential_gates(fcs_data, up_to_gate = 4, gates = gates, channels = channels)

    # Apply Gate 5 - read from gates
    fxcycle_values_temp <- exprs(fcs_data)[, channels$FxCycle]
    fxcycle_limits_temp <- quantile(fxcycle_values_temp, probs = fxcycle_gate$probs, na.rm = TRUE)
    fxcycle_filter <- fxcycle_values_temp >= fxcycle_limits_temp[1] & fxcycle_values_temp <= fxcycle_limits_temp[2]
    fcs_data <- Subset(fcs_data, fxcycle_filter)

    # Calculate thresholds - read from gates
    edu_values <- exprs(fcs_data)[, channels$EdU]
    edu_threshold <- quantile(edu_values, probs = edu_gate$edu_prob, na.rm = TRUE)

    fxcycle_values <- exprs(fcs_data)[, channels$FxCycle]
    fxcycle_bounds <- quantile(fxcycle_values, probs = edu_gate$fxcycle_probs, na.rm = TRUE)
    
    x <- exprs(fcs_data)[, channels$FxCycle]
    y <- exprs(fcs_data)[, channels$EdU]
    
    dens <- get_density_colors(x, y, log_x = FALSE, log_y = TRUE)
    
    plot(x, y,
         pch = ".",
         col = dens,
         xlab = "",
         ylab = "",
         main = sample_name,
         xlim = c(0, 12e6),
         ylim = c(100, 6e6),
         cex.main = 1.2,
         cex.axis = 0.7,
         xaxt = "n",
         yaxt = "n",
         log = "y")
    
    axis(1, at = seq(0, 12e6, 3e6), labels = format_axis_labels(seq(0, 12e6, 3e6)), cex.axis = 0.5)
    axis(2, at = c(100, 1000, 10000, 100000, 1000000), labels = c("100", "1K", "10K", "100K", "1M"), cex.axis = 0.5)
    
    # Add gate lines (only in gated region)
    segments(x0 = fxcycle_bounds[1], y0 = edu_threshold, 
             x1 = fxcycle_bounds[2], y1 = edu_threshold, 
             col = "black", lwd = 2, lty = 2)  # Horizontal line between FxCycle bounds
    segments(x0 = fxcycle_bounds[1], y0 = edu_threshold, 
             x1 = fxcycle_bounds[1], y1 = 6e6, 
             col = "black", lwd = 2, lty = 2)  # Left vertical line
    segments(x0 = fxcycle_bounds[2], y0 = edu_threshold, 
             x1 = fxcycle_bounds[2], y1 = 6e6, 
             col = "black", lwd = 2, lty = 2)  # Right vertical line
    
    inside_gate <- sum(edu_values >= edu_threshold & 
                         fxcycle_values >= fxcycle_bounds[1] & 
                         fxcycle_values <= fxcycle_bounds[2])
    pct <- 100 * inside_gate / nrow(fcs_data)
    text_col <- ifelse(pct < 30, "red", "black")
    text(10e6, 200, sprintf("%.1f%%", pct), col = text_col, cex = 0.8, font = 2, pos = 2)
    
    if(pct < 30) {
      box(col = "red", lwd = 4)
      text(6e6, 10000, "LOW\nCELL COUNT", col = "red", cex = 1.5, font = 2)
    }
  }
  
  par(mfrow = c(1, 1), mar = c(5, 4, 4, 2), oma = c(0, 0, 0, 0))
}

## Gate 7: HA-Positive ----

plot_ha_gate_single <- function(fcs_data, sample_name, ha_threshold, gates = GATES, channels = CHANNELS, show_sample_name = TRUE) {
  # Apply Gates 1-4
  fcs_data <- apply_sequential_gates(fcs_data, up_to_gate = 4, gates = gates, channels = channels)

  # Apply Gate 5: FxCycle quantile - read from gates
  fxcycle_gate <- gates$fxcycle_quantile
  fxcycle_values <- exprs(fcs_data)[, channels$FxCycle]
  fxcycle_limits <- quantile(fxcycle_values, probs = fxcycle_gate$probs, na.rm = TRUE)
  fxcycle_filter <- fxcycle_values >= fxcycle_limits[1] & fxcycle_values <= fxcycle_limits[2]
  fcs_data <- Subset(fcs_data, fxcycle_filter)

  # Apply Gate 6: EdU + FxCycle - read from gates
  edu_gate <- gates$edu_fxcycle_sphase
  edu_values <- exprs(fcs_data)[, channels$EdU]
  edu_threshold <- quantile(edu_values, probs = edu_gate$edu_prob, na.rm = TRUE)

  fxcycle_values2 <- exprs(fcs_data)[, channels$FxCycle]
  fxcycle_bounds <- quantile(fxcycle_values2, probs = edu_gate$fxcycle_probs, na.rm = TRUE)
  
  edu_fxcycle_filter <- edu_values >= edu_threshold & 
    fxcycle_values2 >= fxcycle_bounds[1] & 
    fxcycle_values2 <= fxcycle_bounds[2]
  fcs_data <- Subset(fcs_data, edu_fxcycle_filter)
  
  # Plot settings
  x <- exprs(fcs_data)[, channels$HA]
  y <- exprs(fcs_data)[, channels$EdU]

  # Check if we have enough cells to plot
  if(length(x) < 10) {
    # Create empty plot with warning message
    plot(1, 1, type = "n", xlim = c(100, 1e6), ylim = c(100, 6e6),
         xlab = "HA-A", ylab = "EdU-A",
         main = if(show_sample_name) sprintf("Gate 7: HA-Positive\n%s", sample_name) else "Gate 7: HA-Positive",
         log = "xy", xaxt = "n", yaxt = "n", xaxs = "i", yaxs = "i", mgp = c(3, 0.5, 0))
    axis(1, at = c(100, 1000, 10000, 100000, 1000000),
         labels = c("100", "1K", "10K", "100K", "1M"), mgp = c(3, 0.5, 0))
    axis(2, at = c(100, 1000, 10000, 100000, 1000000),
         labels = c("100", "1K", "10K", "100K", "1M"), mgp = c(3, 0.5, 0))
    text(1000, 50000, sprintf("Insufficient cells for plotting\n(n = %d)", length(x)),
         cex = 1.2, col = "red", font = 2)
    return(invisible(NULL))
  }

  dens <- get_density_colors(x, y, log_x = TRUE, log_y = TRUE)

  # Create title
  plot_title <- if(show_sample_name) {
    sprintf("Gate 7: HA-Positive\n%s", sample_name)
  } else {
    "Gate 7: HA-Positive"
  }

  plot(x, y,
       pch = 16,
       cex = 0.3,
       col = dens,
       xlab = "HA-A",
       ylab = "EdU-A",
       main = plot_title,
       xlim = c(100, 1e6),
       ylim = c(100, 6e6),
       xaxt = "n",
       yaxt = "n",
       mgp = c(3, 0.5, 0),
       xaxs = "i",
       yaxs = "i",
       log = "xy")
  
  axis(1, at = c(100, 1000, 10000, 100000, 1000000),
       labels = c("100", "1K", "10K", "100K", "1M"), mgp = c(3, 0.5, 0))
  axis(2, at = c(100, 1000, 10000, 100000, 1000000),
       labels = c("100", "1K", "10K", "100K", "1M"), mgp = c(3, 0.5, 0))

  # Tint the HA-positive region (only if valid threshold)
  if(!is.null(ha_threshold) && length(ha_threshold) > 0 && !is.na(ha_threshold)) {
    rect(xleft = ha_threshold, ybottom = 100, xright = 1e7, ytop = 6e6,
         col = rgb(0.5, 0.7, 1, 0.25), border = NA)

    # Add vertical threshold line
    abline(v = ha_threshold, col = "black", lwd = 2, lty = 2)
  }
  
  # Calculate stats
  ha_values <- exprs(fcs_data)[, channels$HA]
  total_cells <- nrow(fcs_data)
  inside_gate <- if(!is.null(ha_threshold) && length(ha_threshold) > 0 && !is.na(ha_threshold)) {
    sum(ha_values >= ha_threshold)
  } else {
    0
  }

  # Add legend
  legend_text <- c(sprintf("Total: %s", format(total_cells, big.mark = ",")),
                   sprintf("HA+: %s (%.1f%%)",
                           format(inside_gate, big.mark = ","),
                           100 * inside_gate / total_cells))

  if(!is.null(ha_threshold) && length(ha_threshold) > 0 && !is.na(ha_threshold)) {
    legend_text <- c(legend_text, sprintf("Threshold: %.0fK", ha_threshold/1e3))
  }

  legend("bottomright",
         legend = legend_text,
         bty = "n",
         cex = 0.9)
}

plot_ha_gate_overview <- function(experiment, ha_threshold, gates = GATES, channels = CHANNELS) {
  n_samples <- length(experiment$flowset)
  n_cols <- ceiling(sqrt(n_samples))
  n_rows <- ceiling(n_samples / n_cols)

  par(mfrow = c(n_rows, n_cols), mar = c(2, 2, 1.5, 0.5), oma = c(0, 0, 0, 0), mgp = c(3, 0.3, 0))

  # Read gate parameters
  fxcycle_gate <- gates$fxcycle_quantile
  edu_gate <- gates$edu_fxcycle_sphase

  for(i in seq_along(experiment$flowset)) {
    fcs_data <- experiment$flowset[[i]]
    sample_name <- experiment$metadata$sample_name[i]

    # Apply Gates 1-4
    fcs_data <- apply_sequential_gates(fcs_data, up_to_gate = 4, gates = gates, channels = channels)

    # Apply Gate 5: FxCycle quantile - read from gates
    fxcycle_values <- exprs(fcs_data)[, channels$FxCycle]
    fxcycle_limits <- quantile(fxcycle_values, probs = fxcycle_gate$probs, na.rm = TRUE)
    fxcycle_filter <- fxcycle_values >= fxcycle_limits[1] & fxcycle_values <= fxcycle_limits[2]
    fcs_data <- Subset(fcs_data, fxcycle_filter)

    # Apply Gate 6: EdU + FxCycle - read from gates
    edu_values <- exprs(fcs_data)[, channels$EdU]
    edu_threshold <- quantile(edu_values, probs = edu_gate$edu_prob, na.rm = TRUE)

    fxcycle_values2 <- exprs(fcs_data)[, channels$FxCycle]
    fxcycle_bounds <- quantile(fxcycle_values2, probs = edu_gate$fxcycle_probs, na.rm = TRUE)
    
    edu_fxcycle_filter <- edu_values >= edu_threshold & 
      fxcycle_values2 >= fxcycle_bounds[1] & 
      fxcycle_values2 <= fxcycle_bounds[2]
    fcs_data <- Subset(fcs_data, edu_fxcycle_filter)
    
    x <- exprs(fcs_data)[, channels$HA]
    y <- exprs(fcs_data)[, channels$EdU]

    # Check if we have enough cells to plot
    if(length(x) < 10) {
      # Create empty plot with warning
      plot(1, 1, type = "n", xlim = c(100, 1e6), ylim = c(100, 6e6),
           xlab = "", ylab = "", main = sample_name, log = "xy",
           cex.main = 1.2, cex.axis = 0.7, xaxt = "n", yaxt = "n")
      axis(1, at = c(100, 1000, 10000, 100000, 1000000),
           labels = c("100", "1K", "10K", "100K", "1M"), cex.axis = 0.5)
      axis(2, at = c(100, 1000, 10000, 100000, 1000000),
           labels = c("100", "1K", "10K", "100K", "1M"), cex.axis = 0.5)
      text(1000, 50000, sprintf("n = %d", length(x)), cex = 0.8, col = "red", font = 2)
      next
    }

    dens <- get_density_colors(x, y, log_x = TRUE, log_y = TRUE)

    plot(x, y,
         pch = ".",
         col = dens,
         xlab = "",
         ylab = "",
         main = sample_name,
         xlim = c(100, 1e6),
         ylim = c(100, 6e6),
         cex.main = 1.2,
         cex.axis = 0.7,
         xaxt = "n",
         yaxt = "n",
         log = "xy")
    
    axis(1, at = c(100, 1000, 10000, 100000, 1000000),
         labels = c("100", "1K", "10K", "100K", "1M"), cex.axis = 0.5)
    axis(2, at = c(100, 1000, 10000, 100000, 1000000),
         labels = c("100", "1K", "10K", "100K", "1M"), cex.axis = 0.5)

    # Tint the HA-positive region (only if valid threshold)
    if(!is.null(ha_threshold) && length(ha_threshold) > 0 && !is.na(ha_threshold)) {
      rect(xleft = ha_threshold, ybottom = 100, xright = 1e7, ytop = 6e6,
           col = rgb(0.5, 0.7, 1, 0.25), border = NA)

      # Add threshold line
      abline(v = ha_threshold, col = "black", lwd = 1, lty = 2)
    }
    
    ha_values <- exprs(fcs_data)[, channels$HA]
    inside_gate <- if(!is.null(ha_threshold) && length(ha_threshold) > 0 && !is.na(ha_threshold)) {
      sum(ha_values >= ha_threshold)
    } else {
      0
    }
    pct <- 100 * inside_gate / nrow(fcs_data)
    
    # Determine sample type
    is_dox_minus <- grepl("Dox-", sample_name, ignore.case = TRUE)
    is_empty_vector <- grepl("Empty_Vector", sample_name, ignore.case = TRUE)
    is_dox_plus <- grepl("Dox\\+", sample_name, ignore.case = TRUE)
    
    # Flag based on sample type
    is_flagged <- FALSE
    flag_message <- ""
    
    if(is_dox_minus || is_empty_vector) {
      # Dox- and Empty Vector: flag if TOO HIGH (>5%)
      if(pct > 5) {
        is_flagged <- TRUE
        flag_message <- "HIGH\nHA+"
      }
    } else if(is_dox_plus) {
      # Dox+: flag if TOO LOW (<500 cells)
      if(inside_gate < 500) {
        is_flagged <- TRUE
        flag_message <- "LOW\nCELL COUNT"
      }
    }
    
    text_col <- ifelse(is_flagged, "red", "black")
    text(8e5, 200, sprintf("%.1f%%", pct), col = text_col, cex = 0.8, font = 2, pos = 2)
    
    if(is_flagged) {
      box(col = "red", lwd = 4)
      text(10000, 10000, flag_message, col = "red", cex = 1.5, font = 2)
    }
  }
  
  par(mfrow = c(1, 1), mar = c(5, 4, 4, 2), oma = c(0, 0, 0, 0))
}

## Gate 8: Final EdU vs HA Correlation ----

plot_edu_ha_correlation_single <- function(fcs_data, sample_name, ha_threshold, gates = GATES, channels = CHANNELS, show_sample_name = TRUE, edu_threshold = NULL) {
  # Apply Gates 1-6
  fcs_data <- apply_sequential_gates(fcs_data, up_to_gate = 4, gates = gates, channels = channels)

  # Gate 5: FxCycle quantile - read from gates
  fxcycle_gate <- gates$fxcycle_quantile
  fxcycle_values <- exprs(fcs_data)[, channels$FxCycle]
  fxcycle_limits <- quantile(fxcycle_values, probs = fxcycle_gate$probs, na.rm = TRUE)
  fxcycle_filter <- fxcycle_values >= fxcycle_limits[1] & fxcycle_values <= fxcycle_limits[2]
  fcs_data <- Subset(fcs_data, fxcycle_filter)

  # Gate 6: EdU + FxCycle - read from gates
  edu_gate <- gates$edu_fxcycle_sphase
  edu_values <- exprs(fcs_data)[, channels$EdU]
  edu_threshold_g6 <- quantile(edu_values, probs = edu_gate$edu_prob, na.rm = TRUE)

  fxcycle_values2 <- exprs(fcs_data)[, channels$FxCycle]
  fxcycle_bounds <- quantile(fxcycle_values2, probs = edu_gate$fxcycle_probs, na.rm = TRUE)

  edu_fxcycle_filter <- edu_values >= edu_threshold_g6 &
    fxcycle_values2 >= fxcycle_bounds[1] &
    fxcycle_values2 <= fxcycle_bounds[2]
  fcs_data <- Subset(fcs_data, edu_fxcycle_filter)

  # Gate 7: HA-positive (skip if edu_threshold provided - show quadrants instead)
  if(is.null(edu_threshold)) {
    ha_values <- exprs(fcs_data)[, channels$HA]
    ha_filter <- ha_values >= ha_threshold
    fcs_data <- Subset(fcs_data, ha_filter)
  }

  # Extract final data
  ha_final <- exprs(fcs_data)[, channels$HA]
  edu_final <- exprs(fcs_data)[, channels$EdU]

  # Check if we have enough cells to plot
  if(length(ha_final) < 10) {
    # Create empty plot with warning message
    plot(1, 1, type = "n", xlim = c(2, 6.5), ylim = c(4, 7),
         xlab = "log10(HA-A)", ylab = "log10(EdU-A)",
         main = sprintf("%s\nEdU vs HA Correlation", sample_name),
         xaxs = "i", yaxs = "i")
    text(4.25, 5.5, sprintf("Insufficient cells for plotting\n(n = %d)", length(ha_final)),
         cex = 1.2, col = "red", font = 2)

    return(invisible(list(
      sample_name = sample_name,
      correlation = NA,
      n_cells = length(ha_final),
      ha_log = numeric(0),
      edu_log = numeric(0)
    )))
  }

  # Log10 transform
  ha_log <- log10(ha_final + 1)
  edu_log <- log10(edu_final + 1)

  # Calculate correlation
  correlation <- cor(ha_log, edu_log, use = "complete.obs")

  # Calculate linear regression for trend line
  lm_fit <- lm(edu_log ~ ha_log)
  r_squared <- summary(lm_fit)$r.squared

  # Plot
  dens <- densCols(ha_log, edu_log, colramp = colorRampPalette(c("blue", "cyan", "yellow", "red")))

  par(mgp = c(3, 0.7, 0))

  # Create title with sample name and subtitle
  plot_title <- sprintf("%s\nEdU vs HA Correlation", sample_name)

  plot(ha_log, edu_log,
       pch = 16,
       cex = 0.5,
       col = dens,
       xlab = "log10(HA-A)",
       ylab = "log10(EdU-A)",
       main = plot_title,
       xlim = c(2, 6.5),
       ylim = c(4, 7),

       xaxs = "i",
       yaxs = "i")

  # Add threshold lines if edu_threshold is provided (quadrant mode)
  if(!is.null(edu_threshold)) {
    # Calculate quadrant populations (correct convention: Q1=top right, counterclockwise)
    ha_threshold_log <- log10(ha_threshold + 1)
    edu_threshold_log <- log10(edu_threshold + 1)

    total_cells <- length(ha_log)
    q1_ha_pos_edu_high <- sum(ha_log >= ha_threshold_log & edu_log >= edu_threshold_log)  # Top right
    q2_ha_neg_edu_high <- sum(ha_log < ha_threshold_log & edu_log >= edu_threshold_log)   # Top left
    q3_ha_neg_edu_low <- sum(ha_log < ha_threshold_log & edu_log < edu_threshold_log)     # Bottom left
    q4_ha_pos_edu_low <- sum(ha_log >= ha_threshold_log & edu_log < edu_threshold_log)    # Bottom right

    q1_pct <- (q1_ha_pos_edu_high / total_cells) * 100
    q2_pct <- (q2_ha_neg_edu_high / total_cells) * 100
    q3_pct <- (q3_ha_neg_edu_low / total_cells) * 100
    q4_pct <- (q4_ha_pos_edu_low / total_cells) * 100

    # Calculate strength ratio: Q4 / (Q1 + Q4)
    strength_ratio <- if((q1_pct + q4_pct) > 0) {
      q4_pct / (q1_pct + q4_pct)
    } else {
      NA_real_
    }

    # Add colored background for quadrants
    # Q1 (top right, HA+/EdU-high) = blue with transparency
    rect(ha_threshold_log, edu_threshold_log, 6.5, 7,
         col = rgb(0, 0, 1, alpha = 0.1), border = NA)

    # Q4 (bottom right, HA+/EdU-low) = red with transparency
    rect(ha_threshold_log, 4, 6.5, edu_threshold_log,
         col = rgb(1, 0, 0, alpha = 0.1), border = NA)

    # Re-plot points on top of colored rectangles
    points(ha_log, edu_log, pch = 16, cex = 0.5, col = dens)

    # Draw threshold lines (black)
    abline(v = ha_threshold_log, col = "black", lwd = 2, lty = 1)
    abline(h = edu_threshold_log, col = "black", lwd = 2, lty = 1)

    # Add quadrant percentages
    text(2.8, 6.5, sprintf("%.1f%%", q2_pct), col = "black", cex = 0.8, font = 2)  # Top left (Q2)
    text(5.5, 6.5, sprintf("%.1f%%", q1_pct), col = "blue", cex = 0.8, font = 2)   # Top right (Q1) - blue
    text(2.8, 4.5, sprintf("%.1f%%", q3_pct), col = "black", cex = 0.8, font = 2)  # Bottom left (Q3)
    text(5.5, 4.5, sprintf("%.1f%%", q4_pct), col = "red", cex = 0.8, font = 2)    # Bottom right (Q4) - red

    # Display strength ratio
    if(!is.na(strength_ratio)) {
      text(4.25, 4.2, sprintf("Strength Ratio = %.3f", strength_ratio),
           col = "black", cex = 1, font = 2, pos = 3)
    }
  } else {
    # Add regression line (only in non-quadrant mode)
    abline(lm_fit, col = "black", lwd = 1.5, lty=2)
  }

  # Add r and n label at top center (bold)
  text(4.25, 6.85, sprintf("r = %.3f, n = %s", correlation, format(length(ha_log), big.mark = ",")),
       cex = 1, font = 2, col = "black")

  # Add R² in top left
  text(2.15, 6.85, sprintf("R² = %.3f", r_squared), col = "black", cex = 0.9, font = 2, pos = 4)

  # Add correlation info at bottom right (keep for single plots)
  legend("bottomright",
         legend = c(sprintf("Pearson r = %.3f", correlation),
                    sprintf("n = %s cells", format(length(ha_log), big.mark = ","))),
         bty = "n",
         cex = 1)

  # Return correlation data
  invisible(list(
    sample_name = sample_name,
    correlation = correlation,
    n_cells = length(ha_log),
    ha_log = ha_log,
    edu_log = edu_log
  ))
}

plot_edu_ha_correlation_overview <- function(experiment, ha_threshold, gates = GATES, channels = CHANNELS) {
  # Count only Dox+ samples
  n_dox_plus <- sum(!grepl("Dox-", experiment$metadata$sample_name, ignore.case = TRUE))
  n_cols <- ceiling(sqrt(n_dox_plus))
  n_rows <- ceiling(n_dox_plus / n_cols)

  par(mfrow = c(n_rows, n_cols), mar = c(2, 2.5, 2, 0.5), oma = c(0, 0, 0, 0), mgp = c(3, 0.6, 0))

  # Read gate parameters
  fxcycle_gate <- gates$fxcycle_quantile
  edu_gate <- gates$edu_fxcycle_sphase

  for(i in seq_along(experiment$flowset)) {
    fcs_data <- experiment$flowset[[i]]
    sample_name <- experiment$metadata$sample_name[i]

    # Skip Dox- samples entirely
    is_dox_minus <- grepl("Dox-", sample_name, ignore.case = TRUE)
    if(is_dox_minus) {
      next
    }

    # Apply all gates
    fcs_data <- apply_sequential_gates(fcs_data, up_to_gate = 4, gates = gates, channels = channels)

    # Gate 5 - read from gates
    fxcycle_values <- exprs(fcs_data)[, channels$FxCycle]
    fxcycle_limits <- quantile(fxcycle_values, probs = fxcycle_gate$probs, na.rm = TRUE)
    fxcycle_filter <- fxcycle_values >= fxcycle_limits[1] & fxcycle_values <= fxcycle_limits[2]
    fcs_data <- Subset(fcs_data, fxcycle_filter)

    # Gate 6 - read from gates
    edu_values <- exprs(fcs_data)[, channels$EdU]
    edu_threshold <- quantile(edu_values, probs = edu_gate$edu_prob, na.rm = TRUE)

    fxcycle_values2 <- exprs(fcs_data)[, channels$FxCycle]
    fxcycle_bounds <- quantile(fxcycle_values2, probs = edu_gate$fxcycle_probs, na.rm = TRUE)
    
    edu_fxcycle_filter <- edu_values >= edu_threshold & 
      fxcycle_values2 >= fxcycle_bounds[1] & 
      fxcycle_values2 <= fxcycle_bounds[2]
    fcs_data <- Subset(fcs_data, edu_fxcycle_filter)
    
    # Gate 7: HA-positive
    ha_values <- exprs(fcs_data)[, channels$HA]
    ha_filter <- ha_values >= ha_threshold
    fcs_data <- Subset(fcs_data, ha_filter)
    
    # Get final data
    ha_final <- exprs(fcs_data)[, channels$HA]
    edu_final <- exprs(fcs_data)[, channels$EdU]
    
    # Log transform
    ha_log <- log10(ha_final + 1)
    edu_log <- log10(edu_final + 1)
    
    # Calculate correlation
    if(length(ha_log) > 10) {
      correlation <- cor(ha_log, edu_log, use = "complete.obs")
      lm_fit <- lm(edu_log ~ ha_log)
      r_squared <- summary(lm_fit)$r.squared

      dens <- densCols(ha_log, edu_log, colramp = colorRampPalette(c("blue", "cyan", "yellow", "red")))

      plot(ha_log, edu_log,
           pch = ".",
           col = dens,
           xlab = "",
           ylab = "",
           main = sprintf("%s\nEdU vs HA Correlation", sample_name),
           xlim = c(2, 6.5),
           ylim = c(4, 7),
           cex.main = 1.1,
           xaxs = "i",
           yaxs = "i")

      abline(lm_fit, col = "black", lwd = 1.5, lty =2)

      # Add r and n at top center
      text(4.25, 6.8, sprintf("r=%.3f, n=%s", correlation, format(length(ha_log), big.mark = ",")),
           col = "black", cex = 0.7, font = 2)

      # Add R² in top left
      text(2.15, 6.8, sprintf("R²=%.3f", r_squared), col = "black", cex = 0.7, font = 2, pos = 4)
      
      # Flag if low cell count or extreme correlation
      is_empty_vector <- grepl("Empty_Vector", sample_name, ignore.case = TRUE)
      is_flagged <- length(ha_log) < 500 && !is_empty_vector
      if(is_flagged) {
        box(col = "red", lwd = 4)
        text(5, 5.5, "LOW\nCELL COUNT", col = "red", cex = 1, font = 2)
      }
    } else {
      # Too few cells to plot
      plot.new()
      text(0.5, 0.5, "Too few\ncells", col = "red", cex = 1.5, font = 2)
      box(col = "red", lwd = 4)
    }
  }
  
  par(mfrow = c(1, 1), mar = c(5, 4, 4, 2), oma = c(0, 0, 0, 0))
}

## Table (Extract correlation results for all samples) ----

extract_correlations <- function(experiment, ha_threshold, gates = GATES, channels = CHANNELS) {
  results_list <- list()
  
  for(i in seq_along(experiment$flowset)) {
    fcs_data <- experiment$flowset[[i]]
    sample_name <- experiment$metadata$sample_name[i]
    well <- experiment$metadata$well[i]
    
    # Parse sample name components
    # Extract cell line (first 3 digits)
    cell_line <- str_extract(sample_name, "^\\d{3}")
    
    # Extract gene (SAMD9 or SAMD9L)
    gene <- str_extract(sample_name, "SAMD9L?")
    if(is.na(gene)) gene <- ""
    
    # Extract mutation (everything between gene and _Dox)
    # Remove everything up to and including the gene
    after_gene <- str_replace(sample_name, "^.*?SAMD9L?_", "")
    # Remove _Dox+ or _Dox- at the end
    mutation <- str_replace(after_gene, "_Dox[+-]$", "")
    if(is.na(mutation) || mutation == "") mutation <- ""
    
    # Skip Dox- samples
    is_dox_minus <- grepl("Dox-", sample_name, ignore.case = TRUE)
    
    # Skip Dox- samples
    is_dox_minus <- grepl("Dox-", sample_name, ignore.case = TRUE)
    if(is_dox_minus) {
      next
    }
    
    # Apply all gates
    fcs_data <- apply_sequential_gates(fcs_data, up_to_gate = 4, gates = gates, channels = channels)

    # Gate 5: FxCycle quantile - read from gates
    fxcycle_gate <- gates$fxcycle_quantile
    fxcycle_values <- exprs(fcs_data)[, channels$FxCycle]
    fxcycle_limits <- quantile(fxcycle_values, probs = fxcycle_gate$probs, na.rm = TRUE)
    fxcycle_filter <- fxcycle_values >= fxcycle_limits[1] & fxcycle_values <= fxcycle_limits[2]
    fcs_data <- Subset(fcs_data, fxcycle_filter)

    # Gate 6: EdU + FxCycle - read from gates
    edu_gate <- gates$edu_fxcycle_sphase
    edu_values <- exprs(fcs_data)[, channels$EdU]
    edu_threshold <- quantile(edu_values, probs = edu_gate$edu_prob, na.rm = TRUE)

    fxcycle_values2 <- exprs(fcs_data)[, channels$FxCycle]
    fxcycle_bounds <- quantile(fxcycle_values2, probs = edu_gate$fxcycle_probs, na.rm = TRUE)

    edu_fxcycle_filter <- edu_values >= edu_threshold &
      fxcycle_values2 >= fxcycle_bounds[1] &
      fxcycle_values2 <= fxcycle_bounds[2]
    fcs_data <- Subset(fcs_data, edu_fxcycle_filter)

    # Gate 7: HA-positive
    # Count cells before Gate 7 for HA+ percentage calculation
    n_before_ha_gate <- nrow(exprs(fcs_data))

    ha_values <- exprs(fcs_data)[, channels$HA]
    ha_filter <- ha_values >= ha_threshold
    fcs_data <- Subset(fcs_data, ha_filter)

    # Count cells after Gate 7
    n_after_ha_gate <- nrow(exprs(fcs_data))

    # Calculate HA+ percentage
    ha_pos_pct <- if(n_before_ha_gate > 0) {
      (n_after_ha_gate / n_before_ha_gate) * 100
    } else {
      NA_real_
    }

    # Extract final data
    ha_final <- exprs(fcs_data)[, channels$HA]
    edu_final <- exprs(fcs_data)[, channels$EdU]
    
    # Calculate correlation
    if(length(ha_final) > 10) {
      ha_log <- log10(ha_final + 1)
      edu_log <- log10(edu_final + 1)
      correlation <- cor(ha_log, edu_log, use = "complete.obs")
      n_cells <- length(ha_final)
    } else {
      correlation <- NA
      n_cells <- length(ha_final)
    }
    
    # Check for flags
    flags <- c()
    
    # Check if low cell count in final gate (excluding Empty Vector)
    is_empty_vector <- grepl("Empty_Vector", sample_name, ignore.case = TRUE)
    if(n_cells < 500 && !is_empty_vector) {
      flags <- c(flags, "Low final cell count")
    }
    
    # Check if correlation is NA (too few cells)
    if(is.na(correlation)) {
      flags <- c(flags, "Too few cells for correlation")
    }
    
    # Combine flags or leave blank
    notes <- if(length(flags) > 0) paste(flags, collapse = "; ") else ""
    
    # Store results
    results_list[[i]] <- data.frame(
      Experiment = experiment$experiment_name,
      Well = well,
      Sample = sample_name,
      Cell_line = cell_line,
      Gene = gene,
      Mutation = mutation,
      Correlation = correlation,
      HA_Pos_Pct = ha_pos_pct,
      N_cells = n_cells,
      Notes = notes,
      stringsAsFactors = FALSE
    )
  }

  # Combine all results
  results_df <- bind_rows(results_list)

  return(results_df)
}

## Quadrant Correlation Overview ----
plot_quadrant_correlation_overview <- function(experiment, gates = GATES, channels = CHANNELS) {
  # Count only Dox+ samples
  n_dox_plus <- sum(!grepl("Dox-", experiment$metadata$sample_name, ignore.case = TRUE))
  n_cols <- ceiling(sqrt(n_dox_plus))
  n_rows <- ceiling(n_dox_plus / n_cols)

  par(mfrow = c(n_rows, n_cols), mar = c(2, 2.5, 2, 0.5), oma = c(0, 0, 0, 0), mgp = c(3, 0.6, 0))

  for(i in seq_along(experiment$flowset)) {
    sample_name <- experiment$metadata$sample_name[i]

    # Skip Dox- samples entirely
    is_dox_minus <- grepl("Dox-", sample_name, ignore.case = TRUE)
    if(is_dox_minus) {
      next
    }

    # Find paired control
    control_idx <- find_paired_control(sample_name, experiment$metadata)
    
    if(is.null(control_idx)) {
      # No paired control - show error plot
      plot.new()
      text(0.5, 0.5, sprintf("No Dox- control\n%s", sample_name), 
           col = "red", cex = 0.8, font = 2)
      box(col = "red", lwd = 2)
      next
    }

    # Calculate quadrant thresholds
    control_fcs <- experiment$flowset[[control_idx]]
    test_fcs <- experiment$flowset[[i]]
    control_name <- experiment$metadata$sample_name[control_idx]

    tryCatch({
      quadrant_result <- calculate_quadrant_from_paired_control(
        control_fcs, test_fcs, control_name, sample_name,
        gates, channels
      )

      # Apply gates 1-6
      fcs_data <- apply_sequential_gates(test_fcs, up_to_gate = 4, gates = gates, channels = channels)

      # Gate 5
      fxcycle_gate <- gates$fxcycle_quantile
      fxcycle_values <- exprs(fcs_data)[, channels$FxCycle]
      fxcycle_limits <- quantile(fxcycle_values, probs = fxcycle_gate$probs, na.rm = TRUE)
      fxcycle_filter <- fxcycle_values >= fxcycle_limits[1] & fxcycle_values <= fxcycle_limits[2]
      fcs_data <- Subset(fcs_data, fxcycle_filter)

      # Gate 6
      edu_gate <- gates$edu_fxcycle_sphase
      edu_values <- exprs(fcs_data)[, channels$EdU]
      edu_threshold_g6 <- quantile(edu_values, probs = edu_gate$edu_prob, na.rm = TRUE)

      fxcycle_values2 <- exprs(fcs_data)[, channels$FxCycle]
      fxcycle_bounds <- quantile(fxcycle_values2, probs = edu_gate$fxcycle_probs, na.rm = TRUE)

      edu_fxcycle_filter <- edu_values >= edu_threshold_g6 &
        fxcycle_values2 >= fxcycle_bounds[1] &
        fxcycle_values2 <= fxcycle_bounds[2]
      fcs_data <- Subset(fcs_data, edu_fxcycle_filter)

      # Get final data (no Gate 7 filter for quadrants - show all cells)
      ha_final <- exprs(fcs_data)[, channels$HA]
      edu_final <- exprs(fcs_data)[, channels$EdU]

      # Check if enough cells
      if(length(ha_final) < 10) {
        plot.new()
        text(0.5, 0.5, sprintf("n = %d", length(ha_final)), col = "red", cex = 0.8, font = 2)
        box(col = "red", lwd = 2)
        next
      }

      # Log transform
      ha_log <- log10(ha_final + 1)
      edu_log <- log10(edu_final + 1)

      ha_threshold_log <- log10(quadrant_result$ha_threshold + 1)
      edu_threshold_log <- log10(quadrant_result$edu_threshold + 1)

      # Calculate quadrant populations
      total_cells <- length(ha_log)
      q1_ha_pos_edu_high <- sum(ha_log >= ha_threshold_log & edu_log >= edu_threshold_log)
      q2_ha_neg_edu_high <- sum(ha_log < ha_threshold_log & edu_log >= edu_threshold_log)
      q3_ha_neg_edu_low <- sum(ha_log < ha_threshold_log & edu_log < edu_threshold_log)
      q4_ha_pos_edu_low <- sum(ha_log >= ha_threshold_log & edu_log < edu_threshold_log)

      q1_pct <- (q1_ha_pos_edu_high / total_cells) * 100
      q2_pct <- (q2_ha_neg_edu_high / total_cells) * 100
      q3_pct <- (q3_ha_neg_edu_low / total_cells) * 100
      q4_pct <- (q4_ha_pos_edu_low / total_cells) * 100

      strength_ratio <- if((q1_pct + q4_pct) > 0) {
        q4_pct / (q1_pct + q4_pct)
      } else {
        NA_real_
      }

      # Calculate density colors
      dens <- densCols(ha_log, edu_log, colramp = colorRampPalette(c("blue", "cyan", "yellow", "red")))

      # Plot
      plot(ha_log, edu_log,
           pch = ".",
           col = dens,
           xlab = "",
           ylab = "",
           main = sprintf("%s\nStrength=%.3f", sample_name, strength_ratio),
           xlim = c(2, 6.5),
           ylim = c(4, 7),
           cex.main = 0.9,
           xaxs = "i",
           yaxs = "i")

      # Add colored quadrants
      rect(ha_threshold_log, edu_threshold_log, 6.5, 7,
           col = rgb(0, 0, 1, alpha = 0.1), border = NA)
      rect(ha_threshold_log, 4, 6.5, edu_threshold_log,
           col = rgb(1, 0, 0, alpha = 0.1), border = NA)

      # Re-plot points
      points(ha_log, edu_log, pch = ".", col = dens)

      # Threshold lines
      abline(v = ha_threshold_log, col = "black", lwd = 1, lty = 1)
      abline(h = edu_threshold_log, col = "black", lwd = 1, lty = 1)

      # Quadrant percentages
      text(2.8, 6.5, sprintf("%.1f%%", q2_pct), col = "black", cex = 0.6, font = 2)
      text(5.5, 6.5, sprintf("%.1f%%", q1_pct), col = "blue", cex = 0.6, font = 2)
      text(2.8, 4.5, sprintf("%.1f%%", q3_pct), col = "black", cex = 0.6, font = 2)
      text(5.5, 4.5, sprintf("%.1f%%", q4_pct), col = "red", cex = 0.6, font = 2)

    }, error = function(e) {
      plot.new()
      text(0.5, 0.5, sprintf("Error\n%s", e$message), col = "red", cex = 0.7)
      box(col = "red", lwd = 2)
    })
  }
}
