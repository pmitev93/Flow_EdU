# ==============================================================================
# GATE ADJUSTMENT MODULE
# Interactive gate editing functionality for Flow GUI
# ==============================================================================

library(shiny)
library(flowCore)

# ==============================================================================
# GATE STORAGE AND MANAGEMENT
# ==============================================================================

# Initialize gate storage
# This will hold custom gates that override the defaults
initialize_gate_storage <- function() {
  list(
    # Structure: experiment_name -> sample_index -> gate_name -> gate_matrix
    custom_gates = list(),
    # Track which scope each gate applies to: "sample", "experiment", "global"
    gate_scope = list(),
    # Store default gates for reset functionality
    default_gates = NULL
  )
}

# Get gate for a specific sample (checks custom gates first, then falls back to default)
get_gate_for_sample <- function(gate_storage, experiment_name, sample_idx, gate_name, default_gate) {
  # Check if there's a custom gate for this sample
  if(!is.null(gate_storage$custom_gates[[experiment_name]][[as.character(sample_idx)]][[gate_name]])) {
    return(gate_storage$custom_gates[[experiment_name]][[as.character(sample_idx)]][[gate_name]])
  }
  
  # Check if there's an experiment-level custom gate
  if(!is.null(gate_storage$custom_gates[[experiment_name]][["_experiment_"]][[gate_name]])) {
    return(gate_storage$custom_gates[[experiment_name]][["_experiment_"]][[gate_name]])
  }
  
  # Check if there's a global custom gate
  if(!is.null(gate_storage$custom_gates[["_global_"]][[gate_name]])) {
    return(gate_storage$custom_gates[["_global_"]][[gate_name]])
  }
  
  # Fall back to default
  return(default_gate)
}

# Save a custom gate
save_custom_gate <- function(gate_storage, experiment_name, sample_idx, gate_name, gate_matrix, scope = "sample") {
  if(scope == "global") {
    if(is.null(gate_storage$custom_gates[["_global_"]])) {
      gate_storage$custom_gates[["_global_"]] <- list()
    }
    gate_storage$custom_gates[["_global_"]][[gate_name]] <- gate_matrix
    gate_storage$gate_scope[[gate_name]] <- "global"
  } else if(scope == "experiment") {
    if(is.null(gate_storage$custom_gates[[experiment_name]])) {
      gate_storage$custom_gates[[experiment_name]] <- list()
    }
    if(is.null(gate_storage$custom_gates[[experiment_name]][["_experiment_"]])) {
      gate_storage$custom_gates[[experiment_name]][["_experiment_"]] <- list()
    }
    gate_storage$custom_gates[[experiment_name]][["_experiment_"]][[gate_name]] <- gate_matrix
    gate_storage$gate_scope[[paste0(experiment_name, "_", gate_name)]] <- "experiment"
  } else { # sample
    if(is.null(gate_storage$custom_gates[[experiment_name]])) {
      gate_storage$custom_gates[[experiment_name]] <- list()
    }
    if(is.null(gate_storage$custom_gates[[experiment_name]][[as.character(sample_idx)]])) {
      gate_storage$custom_gates[[experiment_name]][[as.character(sample_idx)]] <- list()
    }
    gate_storage$custom_gates[[experiment_name]][[as.character(sample_idx)]][[gate_name]] <- gate_matrix
    gate_storage$gate_scope[[paste0(experiment_name, "_", sample_idx, "_", gate_name)]] <- "sample"
  }
  
  return(gate_storage)
}

# Reset gate to default
reset_gate_to_default <- function(gate_storage, experiment_name, sample_idx, gate_name) {
  # Remove custom gate at all levels
  if(!is.null(gate_storage$custom_gates[[experiment_name]][[as.character(sample_idx)]])) {
    gate_storage$custom_gates[[experiment_name]][[as.character(sample_idx)]][[gate_name]] <- NULL
  }
  
  if(!is.null(gate_storage$custom_gates[[experiment_name]][["_experiment_"]])) {
    gate_storage$custom_gates[[experiment_name]][["_experiment_"]][[gate_name]] <- NULL
  }
  
  if(!is.null(gate_storage$custom_gates[["_global_"]])) {
    gate_storage$custom_gates[["_global_"]][[gate_name]] <- NULL
  }
  
  return(gate_storage)
}

# ==============================================================================
# GATE EDITING UI COMPONENTS
# ==============================================================================

# Create UI for gate editing panel
gate_edit_ui <- function(gate_name, gate_id) {
  tagList(
    h4(paste("Edit", gate_name)),
    fluidRow(
      column(6,
             actionButton(paste0("edit_", gate_id), "Edit Gate", 
                          class = "btn-primary", icon = icon("edit"))
      ),
      column(6,
             selectInput(paste0("scope_", gate_id), "Apply to:",
                         choices = c("This sample only" = "sample",
                                     "This experiment" = "experiment",
                                     "All experiments" = "global"),
                         selected = "sample")
      )
    ),
    br(),
    fluidRow(
      column(6,
             actionButton(paste0("reset_", gate_id), "Reset to Default", 
                          class = "btn-warning", icon = icon("undo"),
                          style = "width: 100%;")
      ),
      column(6,
             actionButton("deselect_vertex", "Deselect Vertex",
                          class = "btn-info", icon = icon("times-circle"),
                          style = "width: 100%;")
      )
    ),
    hr(),
    
    # Axis controls
    h5("Axis Limits"),
    fluidRow(
      column(6,
             numericInput(paste0("x_min_", gate_id), "X Min:", value = NULL),
             numericInput(paste0("x_max_", gate_id), "X Max:", value = NULL)
      ),
      column(6,
             numericInput(paste0("y_min_", gate_id), "Y Min:", value = NULL),
             numericInput(paste0("y_max_", gate_id), "Y Max:", value = NULL)
      )
    ),
    
    hr(),
    
    # EDITABLE VERTEX COORDINATES (appears when edit mode is active)
    conditionalPanel(
      condition = paste0("input.edit_", gate_id, " % 2 == 1"),
      wellPanel(
        style = "background-color: #fff3cd;",
        h5("✏️ Edit Mode Active"),
        p(strong("For polygon gates (1-2): Click near a vertex to select it, then click elsewhere to move it.")),
        p(strong("For rectangular gates (3-4): Adjust the corner coordinates below.")),
        hr(),
        # Dynamic numeric inputs will appear here
        uiOutput(paste0("vertex_editors_", gate_id)),
        hr(),
        fluidRow(
          column(6,
                 actionButton(paste0("apply_", gate_id), "Apply Changes", 
                              class = "btn-success btn-block", icon = icon("check"))
          ),
          column(6,
                 actionButton(paste0("cancel_", gate_id), "Cancel", 
                              class = "btn-secondary btn-block", icon = icon("times"))
          )
        )
      )
    ),
    
    hr(),
    actionButton("recalculate_all", "Recalculate Results", 
                 class = "btn-danger btn-lg", icon = icon("sync"))
  )
}

# ==============================================================================
# INTERACTIVE GATE PLOTTING
# ==============================================================================

# Plot gate with interactive editing capability
plot_gate_interactive <- function(fcs_data, gate_matrix, x_channel, y_channel, 
                                  gate_name, sample_name = "", edit_mode = FALSE, 
                                  selected_vertex = NULL, temp_gate = NULL,
                                  xlim = NULL, ylim = NULL,
                                  x_label = NULL, y_label = NULL,
                                  log_scale = "") {
  
  # Use temp gate if in edit mode, otherwise use provided gate
  current_gate <- if(!is.null(temp_gate)) temp_gate else gate_matrix
  
  # Extract data - keep full dataset for counting, subset for plotting
  x_data_full <- exprs(fcs_data)[, x_channel]
  y_data_full <- exprs(fcs_data)[, y_channel]
  total_cells <- length(x_data_full)
  
  # Check if we have any data
  if(total_cells == 0) {
    plot.new()
    text(0.5, 0.5, "No cells to display", cex = 1.5, col = "red")
    return()
  }
  
  # Subsample for plotting (max 10000 points for performance)
  if(total_cells > 10000) {
    sample_idx <- sample(total_cells, 10000)
    x_data <- x_data_full[sample_idx]
    y_data <- y_data_full[sample_idx]
  } else {
    x_data <- x_data_full
    y_data <- y_data_full
  }
  
  # For density calculation, use log-transformed data if log scale is applied
  x_for_dens <- if(grepl("x", log_scale)) log10(x_data + 1) else x_data
  y_for_dens <- if(grepl("y", log_scale)) log10(y_data + 1) else y_data
  
  # Calculate density for color mapping (with error handling)
  tryCatch({
    dens <- densCols(x_for_dens, y_for_dens, 
                     colramp = colorRampPalette(c("blue", "cyan", "yellow", "red")))
  }, error = function(e) {
    # If density calculation fails, use default color
    dens <- rep("blue", length(x_data))
  })
  
  # Auto-calculate axis limits if not provided
  if(is.null(xlim)) {
    xlim <- range(x_data, current_gate[, 1], na.rm = TRUE)
    xlim <- xlim + c(-0.05, 0.05) * diff(xlim)  # Add 5% padding
  }
  if(is.null(ylim)) {
    ylim <- range(y_data, current_gate[, 2], na.rm = TRUE)
    ylim <- ylim + c(-0.05, 0.05) * diff(ylim)  # Add 5% padding
  }
  
  # Use provided labels or default to channel names
  x_label <- if(!is.null(x_label)) x_label else x_channel
  y_label <- if(!is.null(y_label)) y_label else y_channel
  
  # Create plot with density colors
  plot(x_data, y_data, 
       pch = 16, cex = 0.5, col = dens,
       xlab = x_label, ylab = y_label,
       main = paste(gate_name, if(edit_mode) "- EDIT MODE" else "", 
                    if(sample_name != "") paste0("\n", sample_name) else ""),
       xlim = xlim, ylim = ylim,
       log = log_scale,
       xaxt = "n", yaxt = "n")
  
  # Custom axis formatting
  format_axis <- function(x) {
    ifelse(x == 0, "0", 
           ifelse(x >= 1e6, sprintf("%.0fM", x/1e6), 
                  ifelse(x >= 1e3, sprintf("%.0fK", x/1e3), sprintf("%.0f", x))))
  }
  
  # Add axes with nice formatting
  if(grepl("x", log_scale)) {
    # Log scale x-axis
    axis_vals <- c(100, 1000, 10000, 100000, 1000000)
    axis_labs <- c("100", "1K", "10K", "100K", "1M")
    axis(1, at = axis_vals, labels = axis_labs)
  } else {
    # Linear scale x-axis
    x_at <- pretty(xlim, n = 5)
    axis(1, at = x_at, labels = format_axis(x_at))
  }
  
  if(grepl("y", log_scale)) {
    # Log scale y-axis
    axis_vals <- c(100, 1000, 10000, 100000, 1000000)
    axis_labs <- c("100", "1K", "10K", "100K", "1M")
    axis(2, at = axis_vals, labels = axis_labs)
  } else {
    # Linear scale y-axis
    y_at <- pretty(ylim, n = 5)
    axis(2, at = y_at, labels = format_axis(y_at))
  }
  
  # Draw gate polygon
  polygon(current_gate[, 1], current_gate[, 2], 
          border = if(edit_mode) "red" else "blue", 
          lwd = if(edit_mode) 3 else 2,
          col = if(edit_mode) rgb(1, 0, 0, 0.1) else rgb(0, 0, 1, 0.1))
  
  # If in edit mode, show vertices
  if(edit_mode) {
    # Draw all vertices with numbers
    n_vertices <- nrow(current_gate)
    n_to_show <- if(all(current_gate[1,] == current_gate[n_vertices,])) n_vertices - 1 else n_vertices
    
    # Draw unselected vertices
    points(current_gate[1:n_to_show, 1], current_gate[1:n_to_show, 2], 
           pch = 21, cex = 2.5, bg = "yellow", col = "red", lwd = 2)
    
    # Highlight selected vertex in orange
    if(!is.null(selected_vertex) && selected_vertex <= n_to_show) {
      points(current_gate[selected_vertex, 1], current_gate[selected_vertex, 2],
             pch = 21, cex = 3, bg = "orange", col = "red", lwd = 3)
    }
    
    # Add vertex numbers
    text(current_gate[1:n_to_show, 1], current_gate[1:n_to_show, 2], 
         labels = 1:n_to_show,
         pos = 3, cex = 1, col = "red", font = 2)
  }
  
  # Count cells inside gate using FULL dataset (not subsampled)
  inside <- sum(point.in.polygon(x_data_full, y_data_full, 
                                 current_gate[, 1], 
                                 current_gate[, 2]) > 0)
  
  # Calculate percentage
  percent_inside <- (inside / total_cells) * 100
  
  # Add legend with count and percentage
  legend("bottomright", 
         legend = c(sprintf("Total: %s", format(total_cells, big.mark = ",")),
                    sprintf("Inside gate: %s (%.1f%%)", 
                            format(inside, big.mark = ","), percent_inside)),
         bty = "n", cex = 0.9)
  
  if(edit_mode) {
    legend("topleft",
           legend = "EDITING MODE",
           bty = "n", cex = 1, text.col = "red", text.font = 2)
  }
}

# Find nearest vertex to click position
find_nearest_vertex <- function(click_x, click_y, gate_matrix, plot_usr) {
  if(is.null(click_x) || is.null(click_y)) return(NULL)
  
  # Normalize coordinates to plot scale
  x_range <- plot_usr[2] - plot_usr[1]
  y_range <- plot_usr[4] - plot_usr[3]
  
  # Calculate distances to all vertices
  distances <- sqrt(((gate_matrix[, 1] - click_x) / x_range)^2 + 
                      ((gate_matrix[, 2] - click_y) / y_range)^2)
  
  # Find nearest (excluding last point if it's duplicate of first)
  n_vertices <- nrow(gate_matrix)
  if(all(gate_matrix[1, ] == gate_matrix[n_vertices, ])) {
    # Last point is duplicate, so consider n-1 vertices
    nearest <- which.min(distances[1:(n_vertices-1)])
  } else {
    nearest <- which.min(distances)
  }
  
  # Only return if click is reasonably close (within 5% of plot range)
  if(distances[nearest] < 0.05) {
    return(nearest)
  }
  
  return(NULL)
}

# Update vertex position
update_vertex_position <- function(gate_matrix, vertex_idx, new_x, new_y) {
  if(is.null(vertex_idx) || vertex_idx > nrow(gate_matrix)) return(gate_matrix)
  
  # Update the vertex
  gate_matrix[vertex_idx, 1] <- new_x
  gate_matrix[vertex_idx, 2] <- new_y
  
  # If first and last vertices are the same (closed polygon), update both
  if(vertex_idx == 1 && all(gate_matrix[1, ] == gate_matrix[nrow(gate_matrix), ])) {
    gate_matrix[nrow(gate_matrix), ] <- c(new_x, new_y)
  } else if(vertex_idx == nrow(gate_matrix) && all(gate_matrix[1, ] == gate_matrix[nrow(gate_matrix), ])) {
    gate_matrix[1, ] <- c(new_x, new_y)
  }
  
  return(gate_matrix)
}

# ==============================================================================
# EXPORT/IMPORT GATES
# ==============================================================================

# Export gates to RDS file
export_gates <- function(gate_storage, filepath) {
  saveRDS(gate_storage, file = filepath)
  return(TRUE)
}

# Import gates from RDS file
import_gates <- function(filepath) {
  if(!file.exists(filepath)) {
    warning("Gate file not found")
    return(NULL)
  }
  
  gate_storage <- readRDS(filepath)
  return(gate_storage)
}

# ==============================================================================
# GATE INFO DISPLAY
# ==============================================================================

# Format gate coordinates for display
format_gate_coordinates <- function(gate_matrix) {
  if(is.null(gate_matrix)) return("No gate defined")
  
  coords_text <- ""
  for(i in 1:nrow(gate_matrix)) {
    coords_text <- paste0(coords_text, 
                          sprintf("Vertex %d: (%.2e, %.2e)\n", 
                                  i, gate_matrix[i, 1], gate_matrix[i, 2]))
  }
  
  return(coords_text)
}

# ==============================================================================
# GATE VALIDATION
# ==============================================================================

# Validate that a gate is properly formed
validate_gate <- function(gate_matrix) {
  if(is.null(gate_matrix)) return(FALSE)
  if(!is.matrix(gate_matrix)) return(FALSE)
  if(ncol(gate_matrix) != 2) return(FALSE)
  if(nrow(gate_matrix) < 3) return(FALSE)
  
  # Check for any NA or infinite values
  if(any(is.na(gate_matrix)) || any(is.infinite(gate_matrix))) return(FALSE)
  
  return(TRUE)
}

cat("Gate adjustment module loaded successfully!\n")