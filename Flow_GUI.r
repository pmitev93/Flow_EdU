# ==============================================================================
# FLOW CYTOMETRY ANALYSIS - SHINY APP
# ==============================================================================

library(shiny)
library(flowCore)
library(tidyverse)
library(sp)
library(DT)
library(openxlsx)
library(sortable)
library(plotly)

# Source the master script to load all functions
# Make sure master script path is correct!
source("R_Flow_Claude.r")  # Uncomment and adjust path
source("plotting_functions.r")
source("gate_adjustment_module.r")  # NEW: Gate editing functionality
#source("Claude_testing_script.R") 
#source("Claude_testing_script2.R") 
#source("Final_analysis.R") 
#source("Analyze_everything.R") 

#To run it:
#shiny::runApp("/Users/petar.mitev/Library/CloudStorage/OneDrive-KarolinskaInstitutet/Experiments/Flow/R_New/Flow_GUI.r")
#or
#shiny::runApp("Flow_GUI.r")

# ==============================================================================
# UI
# ==============================================================================

ui <- fluidPage(
  titlePanel("The Mitev EdU Analysis Tool"),
  
  sidebarLayout(
    sidebarPanel(
      width = 4,  # Increased from 3 to give more space
      
      # Folder selection
      h4("1. Select Data Folder"),
      textInput("master_folder", "Master Folder Path:",
                value = "Experiments/"),
      fluidRow(
        column(6, actionButton("browse_folder", "Browse...", class = "btn-secondary btn-block")),
        column(6, actionButton("rescan_experiments", "Rescan Folder", class = "btn-info btn-block"))
      ),
      hr(),
      
      # Experiment selection with inline gating strategy
      h4("2. Select Experiments & Gating Strategy"),
      conditionalPanel(
        condition = "output.experiments_loaded",
        wellPanel(
          style = "background-color: #f8f9fa; padding: 8px;",
          fluidRow(
            column(5,
                   selectInput("default_gate_strategy", "Default Gates:",
                               choices = NULL, width = "100%")
            ),
            column(3,
                   actionButton("apply_default_gates", "Apply All",
                                class = "btn-sm btn-info btn-block",
                                style = "margin-top: 25px; font-size: 11px;")
            ),
            column(4,
                   actionButton("select_all", "Select All Exps",
                                class = "btn-sm btn-block",
                                style = "margin-top: 25px; font-size: 11px;")
            )
          )
        )
      ),
      uiOutput("experiment_selector_with_gates"),
      br(),
      actionButton("analyze_selected", "Analyze Selected Experiments",
                   class = "btn-success btn-block"),
      hr(),
      
      # Sample browser
      h4("3. Browse Samples"),
      selectInput("selected_experiment", "Experiment:", choices = NULL),
      selectInput("browse_gate_strategy", "Gating Strategy:", choices = NULL),
      selectInput("selected_sample", "Sample:", choices = NULL),
      hr(),
      
      # Download
      h4("4. Download Results"),
      downloadButton("download_results", "Download Excel")
    ),
    
    mainPanel(
      width = 8,  # Adjusted from 9 to match sidebar width=4
      
      tabsetPanel(
        id = "main_tabs",
        
        # Welcome tab
        tabPanel("Welcome",
                 h3("Welcome to the Mitev EdU Analysis Tool"),
                 p("This tool processes flow cytometry data with automated gating and correlation analysis."),
                 h4("How to use:"),
                 tags$ol(
                   tags$li("Experiments are automatically loaded from the folder path on startup"),
                   tags$li("Select which experiments to analyze"),
                   tags$li("(Optional) Choose different gating strategies for each experiment"),
                   tags$li("Click 'Analyze Selected Experiments' to process data"),
                   tags$li("View results in the 'Results Table' tab"),
                   tags$li("Browse individual gates in the 'Gate Inspection' tabs"),
                   tags$li("Download results as Excel file")
                 ),
                 p(strong("Tip:"), " Use the 'Rescan Folder' button if you change the experiments folder or add new data."),
                 hr(),
                 h4("Status:"),
                 verbatimTextOutput("status_text")
        ),
        
        # Results table
        tabPanel("Results Table",
                 h3("Correlation Results"),
                 DTOutput("results_table")
        ),
        
        # ADD THIS NEW TAB:
        tabPanel("Overview Plots",
                 h3("Experiment Overview"),
                 selectInput("overview_experiment", "Select Experiment:",
                             choices = NULL),
                 selectInput("overview_gate_strategy", "Gating Strategy:",
                             choices = NULL),
                 selectInput("overview_gate", "Select Gate:",
                             choices = c("Gate 1: Debris" = "gate1",
                                         "Gate 2: Singlets" = "gate2",
                                         "Gate 3: Live Cells" = "gate3",
                                         "Gate 4: S-phase Outliers" = "gate4",
                                         "Gate 5: FxCycle Quantile" = "gate5",
                                         "Gate 6: EdU + FxCycle" = "gate6",
                                         "Gate 7: HA-Positive" = "gate7",
                                         "Final: Correlation" = "correlation")),
                 plotOutput("overview_plot", height = "800px")
        ),
        
        # Gate Editing Tab
        tabPanel("Edit Gates",
                 h3("Interactive Gate Editor"),
                 p("Use this panel to adjust gates for your samples."),
                 
                 fluidRow(
                   column(4,
                          wellPanel(
                            h4("Select Gate to Edit"),
                            selectInput("gate_to_edit", "Gate:",
                                        choices = c("Gate 1: Debris" = "debris",
                                                    "Gate 2: Singlets" = "singlet",
                                                    "Gate 3: Live Cells" = "live_cells",
                                                    "Gate 4: S-phase Outliers" = "s_phase_outliers"),
                                        selected = "debris"),
                            hr(),
                            gate_edit_ui("Selected Gate", "main")
                          )
                   ),
                   column(8,
                          plotOutput("gate_edit_plot", height = "600px",
                                     click = "gate_plot_click",
                                     brush = brushOpts(id = "gate_plot_brush", 
                                                       resetOnNew = TRUE))
                   )
                 )
        ),
        
        # Gate inspection tabs
        tabPanel("Gate 1: Debris",
                 fluidRow(
                   column(2, actionButton("gate1_prev", "Previous", class = "btn-sm")),
                   column(8, h4(textOutput("current_sample_name"), align = "center")),
                   column(2, actionButton("gate1_next", "Next", class = "btn-sm", style = "float: right;"))
                 ),
                 plotOutput("gate1_plot", height = "600px")
        ),
        
        tabPanel("Gate 2: Singlets",
                 fluidRow(
                   column(2, actionButton("gate1_prev", "Previous", class = "btn-sm")),
                   column(8, h4(textOutput("current_sample_name"), align = "center")),
                   column(2, actionButton("gate1_next", "Next", class = "btn-sm", style = "float: right;"))
                 ),
                 plotOutput("gate2_plot", height = "600px")
        ),
        
        tabPanel("Gate 3: Live Cells",
                 fluidRow(
                   column(2, actionButton("gate1_prev", "Previous", class = "btn-sm")),
                   column(8, h4(textOutput("current_sample_name"), align = "center")),
                   column(2, actionButton("gate1_next", "Next", class = "btn-sm", style = "float: right;"))
                 ),
                 plotOutput("gate3_plot", height = "600px")
        ),
        
        tabPanel("Gate 4: S-phase Outliers",
                 fluidRow(
                   column(2, actionButton("gate1_prev", "Previous", class = "btn-sm")),
                   column(8, h4(textOutput("current_sample_name"), align = "center")),
                   column(2, actionButton("gate1_next", "Next", class = "btn-sm", style = "float: right;"))
                 ),
                 plotOutput("gate4_plot", height = "600px")
        ),
        
        tabPanel("Gate 5: FxCycle Quantile",
                 fluidRow(
                   column(2, actionButton("gate1_prev", "Previous", class = "btn-sm")),
                   column(8, h4(textOutput("current_sample_name"), align = "center")),
                   column(2, actionButton("gate1_next", "Next", class = "btn-sm", style = "float: right;"))
                 ),
                 plotOutput("gate5_plot", height = "600px")
        ),
        
        tabPanel("Gate 6: EdU + FxCycle",
                 fluidRow(
                   column(2, actionButton("gate1_prev", "Previous", class = "btn-sm")),
                   column(8, h4(textOutput("current_sample_name"), align = "center")),
                   column(2, actionButton("gate1_next", "Next", class = "btn-sm", style = "float: right;"))
                 ),
                 plotOutput("gate6_plot", height = "600px")
        ),
        
        tabPanel("Gate 7: HA-Positive",
                 fluidRow(
                   column(2, actionButton("gate1_prev", "Previous", class = "btn-sm")),
                   column(8, h4(textOutput("current_sample_name"), align = "center")),
                   column(2, actionButton("gate1_next", "Next", class = "btn-sm", style = "float: right;"))
                 ),
                 plotOutput("gate7_plot", height = "600px")
        ),
        
        tabPanel("Final: Correlation",
                 fluidRow(
                   column(2, actionButton("gate1_prev", "Previous", class = "btn-sm")),
                   column(8, h4(textOutput("current_sample_name"), align = "center")),
                   column(2, actionButton("gate1_next", "Next", class = "btn-sm", style = "float: right;"))
                 ),
                 plotOutput("correlation_plot", height = "600px")
        ),
        
        tabPanel("Multi-Sample Comparison",
                 h3("Compare Multiple Samples"),
                 
                 fluidRow(
                   column(12,
                          h4("Select Samples to Compare"),
                          DTOutput("comparison_sample_selector"),
                          actionButton("clear_selection", "Clear Selection", class = "btn-sm")
                   )
                 ),
                 
                 fluidRow(
                   column(12,
                          h4("Selected Samples"),
                          verbatimTextOutput("selected_samples_list"),
                          fluidRow(
                            column(3, actionButton("plot_comparison", "Generate Comparison Plot", 
                                                   class = "btn-primary")),
                            column(3, actionButton("reorder_samples", "Reorder Samples", 
                                                   class = "btn-secondary")),
                            column(3, downloadButton("download_prism", "Download for Prism",
                                                     class = "btn-success")),
                            column(3, selectInput("reference_group", "Reference for comparisons:",
                                                  choices = NULL))
                          )
                   )
                 ),
                 
                 fluidRow(
                   column(12,
                          h4("Correlation Comparison"),
                          plotOutput("comparison_plot", height = "600px")
                   )
                 ),
                 
                 hr(),
                 
                 fluidRow(
                   column(12,
                          h4("Statistical Analysis"),
                          verbatimTextOutput("stats_output")
                   )
                 )
        )
      )
    )
  )
)

# ==============================================================================
# SERVER
# ==============================================================================

server <- function(input, output, session) {
  
  # ==============================================================================
  # ANALYSIS CACHE SYSTEM
  # ==============================================================================
  
  # Ensure digest package is available for fingerprinting
  if(!require("digest", quietly = TRUE)) {
    install.packages("digest", repos = "http://cran.r-project.org")
    library(digest)
  }
  
  # Create cache directory if it doesn't exist
  CACHE_DIR <- "analysis_cache"
  if(!dir.exists(CACHE_DIR)) {
    dir.create(CACHE_DIR, recursive = TRUE)
  }
  
  # Generate unique fingerprint for gates (includes HA percentile method, not actual value)
  get_gate_fingerprint <- function(gates, ha_percentile = 0.98) {
    # Combine all gate coordinates into a single string
    gate_strings <- lapply(names(gates), function(gate_name) {
      gate <- gates[[gate_name]]
      if(is.matrix(gate)) {
        # For matrix gates, create string from coordinates
        paste(gate_name, paste(as.vector(gate), collapse = ","), sep = ":")
      } else if(is.list(gate)) {
        # For list gates (like quantile gates), use the parameters
        paste(gate_name, paste(unlist(gate), collapse = ","), sep = ":")
      } else {
        ""
      }
    })
    
    # Add HA percentile method (not the actual threshold value)
    gate_strings <- c(gate_strings, paste0("HA_percentile:", ha_percentile))
    
    # Create MD5 hash of the combined string
    combined <- paste(gate_strings, collapse = "|")
    digest::digest(combined, algo = "md5")
  }
  
  # Convert MD5 hash to human-readable gate ID
  get_readable_gate_id <- function(fingerprint) {
    # Check if this is a known gate ID already
    gate_map_file <- file.path(CACHE_DIR, "gate_id_map.rds")
    
    if(file.exists(gate_map_file)) {
      gate_map <- readRDS(gate_map_file)
    } else {
      gate_map <- list(
        fingerprints = character(0),
        ids = character(0),
        counter = 0
      )
    }
    
    # Check if we've seen this fingerprint before
    idx <- which(gate_map$fingerprints == fingerprint)
    
    if(length(idx) > 0) {
      cat(sprintf("Reusing existing Gate ID: %s for fingerprint %s\n", 
                  gate_map$ids[idx[1]], substr(fingerprint, 1, 8)))
      return(gate_map$ids[idx[1]])
    }
    
    # New fingerprint - assign new ID
    if(gate_map$counter == 0) {
      # First gate strategy gets "gdef" (gate default)
      new_id <- "gdef"
    } else {
      # Subsequent gates get gid1, gid2, etc.
      new_id <- sprintf("gid%d", gate_map$counter)
    }
    
    cat(sprintf("NEW Gate ID: %s for fingerprint %s\n", new_id, substr(fingerprint, 1, 8)))
    
    # Store mapping
    gate_map$fingerprints <- c(gate_map$fingerprints, fingerprint)
    gate_map$ids <- c(gate_map$ids, new_id)
    gate_map$counter <- gate_map$counter + 1
    
    saveRDS(gate_map, gate_map_file)
    
    return(new_id)
  }
  
  # Get cache file path for an experiment
  get_cache_path <- function(experiment_name, fingerprint, ha_threshold, gate_id = "gdef") {
    # Create subfolder for this gating strategy
    gate_dir <- file.path(CACHE_DIR, gate_id)
    if(!dir.exists(gate_dir)) {
      dir.create(gate_dir, recursive = TRUE)
    }

    # Include gate_id and HA threshold value in filename
    file.path(gate_dir, paste0(experiment_name, "_", gate_id, "_", fingerprint, "_",
                                sprintf("%.0f", ha_threshold), ".rds"))
  }
  
  # Save analysis results to cache
  save_to_cache <- function(experiment_name, results, ha_threshold, gates, ha_percentile = 0.98, gate_strategy_id = NULL) {
    fingerprint <- get_gate_fingerprint(gates, ha_percentile)

    # Use provided gate_strategy_id or fallback to auto-generated ID
    gate_id <- if(!is.null(gate_strategy_id)) {
      gate_strategy_id
    } else {
      get_readable_gate_id(fingerprint)
    }

    cache_file <- get_cache_path(experiment_name, fingerprint, ha_threshold, gate_id)

    cache_data <- list(
      results = results,
      ha_threshold = ha_threshold,
      ha_percentile = ha_percentile,
      fingerprint = fingerprint,
      gate_id = gate_id,
      gates = gates,
      timestamp = Sys.time()
    )

    saveRDS(cache_data, cache_file)
    cat(sprintf("Saved cache: %s (Gate ID: %s)\n", basename(cache_file), gate_id))

    # Also save human-readable gate details
    save_gate_details(gate_id, fingerprint, gates, ha_percentile)

    return(gate_id)
  }
  
  # Save gate details to human-readable text file
  save_gate_details <- function(gate_id, fingerprint, gates, ha_percentile) {
    # Create subfolder for this gating strategy
    gate_dir <- file.path(CACHE_DIR, gate_id)
    if(!dir.exists(gate_dir)) {
      dir.create(gate_dir, recursive = TRUE)
    }

    # Use gate ID for filename so same gates = same file
    details_file <- file.path(gate_dir, paste0("gate_", gate_id, ".txt"))
    
    # Only create if doesn't exist (same gates = same file)
    if(file.exists(details_file)) {
      return()
    }
    
    sink(details_file)
    cat("=", rep("=", 78), "\n", sep = "")
    cat("GATE STRATEGY DETAILS\n")
    cat("=", rep("=", 78), "\n", sep = "")
    cat(sprintf("Gate Strategy ID: %s\n", gate_id))
    cat(sprintf("Full Fingerprint: %s\n", fingerprint))
    cat(sprintf("Created: %s\n", Sys.time()))
    cat(sprintf("\nHA Threshold Method: %.0fth percentile of Empty Vector control\n", ha_percentile * 100))
    cat("NOTE: The actual HA threshold value varies by experiment.\n")
    cat("      All experiments using this gate strategy share the same Gate ID.\n")
    cat("\n")
    
    for(gate_name in names(gates)) {
      gate <- gates[[gate_name]]
      
      cat("-", rep("-", 78), "\n", sep = "")
      cat(sprintf("GATE: %s\n", gate_name))
      cat("-", rep("-", 78), "\n", sep = "")
      
      if(is.matrix(gate)) {
        cat("Type: Polygon\n")
        cat(sprintf("Vertices: %d\n", nrow(gate)))
        cat("\nCoordinates:\n")
        for(i in 1:nrow(gate)) {
          cat(sprintf("  Vertex %2d: X = %12.2f, Y = %12.2f\n", i, gate[i, 1], gate[i, 2]))
        }
      } else if(is.list(gate)) {
        cat(sprintf("Type: %s\n", gate$type))
        
        # Show all parameters in the list
        cat("\nParameters:\n")
        for(param_name in names(gate)) {
          if(param_name == "type") next  # Already shown
          
          param_value <- gate[[param_name]]
          
          if(is.numeric(param_value) && length(param_value) == 1) {
            cat(sprintf("  %s: %.6f\n", param_name, param_value))
          } else if(is.numeric(param_value) && length(param_value) > 1) {
            if(param_name == "probs") {
              cat(sprintf("  %s: %.2f%% to %.2f%%\n", param_name, 
                          param_value[1]*100, param_value[2]*100))
            } else {
              cat(sprintf("  %s: [%s]\n", param_name, 
                          paste(sprintf("%.4f", param_value), collapse = ", ")))
            }
          } else if(is.character(param_value)) {
            cat(sprintf("  %s: %s\n", param_name, param_value))
          } else {
            cat(sprintf("  %s: %s\n", param_name, paste(param_value, collapse = ", ")))
          }
        }
      }
      cat("\n")
    }
    
    cat("=", rep("=", 78), "\n", sep = "")
    sink()
    
    cat(sprintf("Saved gate strategy details: %s\n", details_file))
  }
  
  # Load analysis results from cache
  load_from_cache <- function(experiment_name, gates, ha_threshold, ha_percentile = 0.98, gate_id = "gdef") {
    fingerprint <- get_gate_fingerprint(gates, ha_percentile)
    cache_file <- get_cache_path(experiment_name, fingerprint, ha_threshold, gate_id)
    
    if(file.exists(cache_file)) {
      cache_data <- readRDS(cache_file)
      
      # Verify fingerprint matches
      if(cache_data$fingerprint == fingerprint) {
        gate_id <- if(!is.null(cache_data$gate_id)) cache_data$gate_id else get_readable_gate_id(fingerprint)
        cat(sprintf("Loaded from cache: %s (analyzed %s, Gate ID: %s)\n", 
                    experiment_name, cache_data$timestamp, gate_id))
        return(cache_data)
      } else {
        cat(sprintf("Cache fingerprint mismatch for %s\n", experiment_name))
        return(NULL)
      }
    }
    
    return(NULL)
  }
  
  # ==============================================================================
  # REACTIVE VALUES
  # ==============================================================================
  
  # Reactive values to store data
  rv <- reactiveValues(
    experiments = NULL,
    all_results = NULL,
    ha_thresholds = list(),
    experiment_gates = list(),  # Store which gates were used for each experiment (by experiment+gate_id)
    experiment_available_gates = list(),  # Store list of available gate IDs for each experiment
    gate_storage = initialize_gate_storage(),  # Store custom gates
    edit_mode = FALSE,  # Track if in edit mode
    temp_gate = NULL,  # Temporary gate during editing
    selected_vertex = NULL,  # Currently selected vertex
    available_gate_files = NULL,  # List of available gate strategy files
    experiment_gate_strategies = list(),  # Per-experiment gate strategy selection
    experiments_loaded = FALSE,  # Track if experiments have been loaded
    scan_trigger = 0,  # Trigger for rescanning experiments
    ui_refresh_trigger = 0  # Trigger for refreshing experiment UI
  )
  
  # Scan for available gate strategy files
  scan_gate_files <- function() {
    gate_dir <- "gate_definitions"
    if(!dir.exists(gate_dir)) {
      showNotification("gate_definitions folder not found!", type = "error")
      return(NULL)
    }

    gate_files <- list.files(gate_dir, pattern = "^gates_.*\\.r$", full.names = FALSE)

    if(length(gate_files) == 0) {
      showNotification("No gate files found in gate_definitions/", type = "warning")
      return(NULL)
    }

    # Extract gate IDs from filenames (gates_gdef.r -> gdef)
    gate_ids <- gsub("^gates_(.*)\\.r$", "\\1", gate_files)
    names(gate_files) <- gate_ids

    return(gate_files)
  }

  # Output variable for conditional panel
  output$experiments_loaded <- reactive({
    return(rv$experiments_loaded)
  })
  outputOptions(output, "experiments_loaded", suspendWhenHidden = FALSE)

  # Browse for folder
  observeEvent(input$browse_folder, {
    folder <- choose.dir(default = getwd(), caption = "Select Master Folder")
    if(!is.null(folder) && !is.na(folder)) {
      updateTextInput(session, "master_folder", value = folder)
    }
  })
  
  # Auto-scan for experiments on startup and when rescan is triggered
  observeEvent(rv$scan_trigger, ignoreInit = FALSE, {
    # Get folder path (isolate to prevent reactive dependency)
    master_folder <- isolate(input$master_folder)
    if(is.null(master_folder) || master_folder == "") {
      master_folder <- "Experiments/"
    }

    cat(sprintf("=== SCAN TRIGGERED (trigger=%d) ===\n", rv$scan_trigger))
    
    tryCatch({
      exp_folders <- list.dirs(master_folder, recursive = FALSE, full.names = TRUE)
      
      if(length(exp_folders) == 0) {
        showNotification("No experiment folders found!", type = "warning")
        return()
      }
      
      withProgress(message = 'Scanning experiments...', value = 0, {
        # Quick scan all experiments (just metadata, no FCS loading)
        all_metadata <- list()
        for(i in seq_along(exp_folders)) {
          incProgress(1/length(exp_folders))
          all_metadata[[i]] <- quick_scan_experiment(exp_folders[i])
        }
        
        # Combine all metadata
        rv$all_results <- bind_rows(all_metadata)
        
        # Store folder paths for later loading
        exp_names <- basename(exp_folders)
        rv$experiment_folders <- setNames(exp_folders, exp_names)
        
        # Update experiment selector with inline gate selection
        output$experiment_selector_with_gates <- renderUI({
          req(rv$available_gate_files)

          # React to UI refresh trigger to re-render after analysis
          rv$ui_refresh_trigger

          # Isolate exp_names to prevent unwanted re-rendering
          exp_names_local <- isolate(exp_names)

          # Check which experiments have cached analyses
          exp_info <- lapply(exp_names_local, function(exp_name) {
            # Scan all subfolders (gating strategies) for cache files
            pattern <- paste0("^", exp_name, "_.*\\.rds$")
            cache_files <- list.files(CACHE_DIR, pattern = pattern, full.names = TRUE, recursive = TRUE)

            if(length(cache_files) > 0) {
              # Collect all gate IDs for this experiment
              gate_ids <- character(0)
              for(cache_file in cache_files) {
                tryCatch({
                  cache_data <- readRDS(cache_file)
                  gate_id <- if(!is.null(cache_data$gate_id)) {
                    cache_data$gate_id
                  } else {
                    get_readable_gate_id(cache_data$fingerprint)
                  }
                  gate_ids <- c(gate_ids, gate_id)
                }, error = function(e) {
                  # Skip files that can't be read
                })
              }

              # Remove duplicates and sort
              gate_ids <- unique(gate_ids)
              gate_ids <- sort(gate_ids)

              return(list(analyzed = TRUE, gate_ids = gate_ids))
            } else {
              return(list(analyzed = FALSE, gate_ids = character(0)))
            }
          })
          names(exp_info) <- exp_names_local

          # Store available gate strategies for each experiment
          for(exp_name in exp_names_local) {
            rv$experiment_available_gates[[exp_name]] <- exp_info[[exp_name]]$gate_ids
          }

          # Get gate choices (isolate to prevent reactivity)
          gate_choices <- isolate({
            setNames(
              rv$available_gate_files,
              names(rv$available_gate_files)
            )
          })

          # Create rows with checkbox + experiment name + dropdown
          rows <- lapply(seq_along(exp_names_local), function(i) {
            exp_name <- exp_names_local[i]
            info <- exp_info[[exp_name]]

            # Current gate selection (isolate to prevent re-render on change)
            current_gate <- isolate(rv$experiment_gate_strategies[[exp_name]])
            if(is.null(current_gate)) current_gate <- "gates_gdef.r"

            # Create label with checkmark if analyzed
            label_text <- if(info$analyzed && length(info$gate_ids) > 0) {
              paste0(exp_name, " ✓ (", paste(info$gate_ids, collapse = ", "), ")")
            } else {
              exp_name
            }

            div(
              style = "margin-bottom: 5px; padding: 3px; border-bottom: 1px solid #eee; display: flex; align-items: center;",
              div(style = "width: 30px; flex-shrink: 0;",
                  checkboxInput(paste0("exp_check_", i), NULL,
                                value = FALSE, width = "100%")
              ),
              div(style = "flex: 1; min-width: 0; padding: 0 8px;",
                  p(style = "margin: 0; font-size: 13px; white-space: nowrap; overflow: hidden; text-overflow: ellipsis;",
                    title = label_text,  # Show full name on hover
                    strong(label_text))
              ),
              div(style = "width: 120px; flex-shrink: 0;",
                  selectInput(paste0("gate_strategy_", i), NULL,
                              choices = gate_choices,
                              selected = current_gate,
                              width = "100%")
              )
            )
          })

          div(
            style = "max-height: 400px; overflow-y: auto; border: 1px solid #ddd; padding: 10px;",
            rows
          )
        })

        # Scan for available gate files
        rv$available_gate_files <- scan_gate_files()

        # Initialize gate strategies with default (gdef)
        default_gate <- "gates_gdef.r"
        rv$experiment_gate_strategies <- setNames(
          rep(list(default_gate), length(exp_names)),
          exp_names
        )

        # Update default gate strategy dropdown
        if(!is.null(rv$available_gate_files)) {
          gate_choices <- setNames(
            rv$available_gate_files,
            names(rv$available_gate_files)
          )
          updateSelectInput(session, "default_gate_strategy",
                            choices = gate_choices,
                            selected = default_gate)
        }

        # Mark experiments as loaded
        rv$experiments_loaded <- TRUE

        showNotification(sprintf("Found %d experiments with %d samples!",
                                 length(exp_names), nrow(rv$all_results)),
                         type = "message")
      })
      
      # AUTO-LOAD CACHED ANALYSES
      withProgress(message = 'Loading cached analyses...', value = 0, {
        n_loaded <- 0
        
        for(i in seq_along(exp_names)) {
          exp_name <- exp_names[i]
          incProgress(1/length(exp_names), detail = exp_name)
          
          # Find most recent cache file for this experiment (scan all subfolders)
          pattern <- paste0("^", exp_name, "_.*\\.rds$")
          cache_files <- list.files(CACHE_DIR, pattern = pattern, full.names = TRUE, recursive = TRUE)
          
          if(length(cache_files) > 0) {
            # Use the most recent cache file
            cache_times <- file.mtime(cache_files)
            latest_cache <- cache_files[which.max(cache_times)]
            
            tryCatch({
              cache_data <- readRDS(latest_cache)
              
              # Add Gate_ID column if it doesn't exist
              if(!"Gate_ID" %in% names(rv$all_results)) {
                rv$all_results$Gate_ID <- NA_character_
              }
              
              # Update results table with cached data
              for(j in seq_len(nrow(cache_data$results))) {
                match_idx <- which(rv$all_results$Experiment == cache_data$results$Experiment[j] & 
                                     rv$all_results$Well == cache_data$results$Well[j])
                
                if(length(match_idx) > 0) {
                  rv$all_results$Correlation[match_idx] <- cache_data$results$Correlation[j]
                  rv$all_results$N_cells[match_idx] <- cache_data$results$N_cells[j]
                  rv$all_results$Notes[match_idx] <- cache_data$results$Notes[j]
                  rv$all_results$Gate_ID[match_idx] <- if(!is.null(cache_data$gate_id)) cache_data$gate_id else get_readable_gate_id(cache_data$fingerprint)
                }
              }
              
              # Store HA threshold
              rv$ha_thresholds[[exp_name]] <- cache_data$ha_threshold

              # Store gates used for this experiment+strategy combination
              if(!is.null(cache_data$gates)) {
                gate_id <- if(!is.null(cache_data$gate_id)) {
                  cache_data$gate_id
                } else {
                  get_readable_gate_id(cache_data$fingerprint)
                }
                composite_key <- paste0(exp_name, "::", gate_id)
                rv$experiment_gates[[composite_key]] <- cache_data$gates
              }

              # Also load the experiment FCS data for browsing
              if(is.null(rv$experiments)) {
                rv$experiments <- list()
              }
              if(is.null(rv$experiments[[exp_name]])) {
                exp_folder <- exp_folders[i]
                rv$experiments[[exp_name]] <- load_experiment(exp_folder)
              }

              n_loaded <- n_loaded + 1
              
            }, error = function(e) {
              cat(sprintf("Error loading cache for %s: %s\n", exp_name, e$message))
            })
          }
        }
        
        if(n_loaded > 0) {
          showNotification(sprintf("Auto-loaded %d cached analyses", n_loaded),
                           type = "message", duration = 5)

          # Update Browse Samples dropdown with loaded experiments
          if(!is.null(rv$experiments) && length(rv$experiments) > 0) {
            updateSelectInput(session, "selected_experiment",
                              choices = names(rv$experiments))
          }
        }
      })
      
    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error")
    })
  })
  
  # Rescan when button clicked
  observeEvent(input$rescan_experiments, {
    req(input$master_folder)

    # Clear existing data
    rv$experiments <- NULL
    rv$all_results <- NULL
    rv$experiments_loaded <- FALSE

    # Increment trigger to force re-scan
    rv$scan_trigger <- rv$scan_trigger + 1

    showNotification("Rescanning experiments folder...", type = "message", duration = 2)
  })
  
  # Select all experiments
  observeEvent(input$select_all, {
    req(rv$experiment_folders)

    exp_names <- names(rv$experiment_folders)
    for(i in seq_along(exp_names)) {
      updateCheckboxInput(session, paste0("exp_check_", i), value = TRUE)
    }
  })

  # Apply default gate strategy to all experiments
  observeEvent(input$apply_default_gates, {
    req(input$default_gate_strategy, rv$experiment_folders)

    exp_names <- names(rv$experiment_folders)
    default_gate <- input$default_gate_strategy

    # Update all experiment gate strategies
    rv$experiment_gate_strategies <- setNames(
      rep(list(default_gate), length(exp_names)),
      exp_names
    )

    # Update all gate dropdown selects
    for(i in seq_along(exp_names)) {
      updateSelectInput(session, paste0("gate_strategy_", i), selected = default_gate)
    }

    showNotification(sprintf("Applied %s to all experiments",
                             names(rv$available_gate_files)[rv$available_gate_files == default_gate]),
                     type = "message")
  })

  # Track individual gate strategy changes
  observe({
    req(rv$experiment_folders)

    exp_names <- names(rv$experiment_folders)

    lapply(seq_along(exp_names), function(i) {
      exp_name <- exp_names[i]
      input_id <- paste0("gate_strategy_", i)

      observeEvent(input[[input_id]], {
        rv$experiment_gate_strategies[[exp_name]] <- input[[input_id]]
      }, ignoreInit = TRUE)
    })
  })

  # Update sample selector when experiment changes
  observeEvent(input$selected_experiment, {
    req(input$selected_experiment)
    
    # Auto-load experiment if not already loaded
    if(is.null(rv$experiments[[input$selected_experiment]])) {
      req(rv$experiment_folders)
      
      exp_name <- input$selected_experiment
      exp_folder <- rv$experiment_folders[[exp_name]]
      
      if(!is.null(exp_folder)) {
        withProgress(message = paste('Loading', exp_name), value = 0.5, {
          if(is.null(rv$experiments)) {
            rv$experiments <- list()
          }
          rv$experiments[[exp_name]] <- load_experiment(exp_folder)
          
          showNotification(sprintf("Loaded %s for browsing", exp_name), 
                           type = "message", duration = 2)
        })
      }
    }
    
    req(rv$experiments[[input$selected_experiment]])
    
    exp <- rv$experiments[[input$selected_experiment]]
    samples <- exp$metadata$sample_name

    updateSelectInput(session, "selected_sample",
                      choices = setNames(seq_along(samples), samples))

    # Update gating strategy dropdown with available strategies for this experiment
    exp_name <- input$selected_experiment
    available_gates <- rv$experiment_available_gates[[exp_name]]
    if(is.null(available_gates) || length(available_gates) == 0) {
      available_gates <- "gdef"  # Default if none available
    }
    updateSelectInput(session, "browse_gate_strategy",
                      choices = available_gates,
                      selected = available_gates[1])
  })
  
  # Analyze selected experiments (load data now)
  observeEvent(input$analyze_selected, {
    req(rv$experiment_folders)

    # Get selected experiments from checkboxes
    exp_names <- names(rv$experiment_folders)
    selected_exps <- c()
    for(i in seq_along(exp_names)) {
      if(isTRUE(input[[paste0("exp_check_", i)]])) {
        selected_exps <- c(selected_exps, exp_names[i])
      }
    }

    if(length(selected_exps) == 0) {
      showNotification("Please select at least one experiment to analyze", type = "warning")
      return()
    }

    withProgress(message = 'Loading and analyzing experiments...', value = 0, {

      # Load only selected experiments
      if(is.null(rv$experiments)) {
        rv$experiments <- list()
      }

      all_results <- list()
      n_exp <- length(selected_exps)
      n_from_cache <- 0
      n_analyzed <- 0

      for(exp_idx in seq_along(selected_exps)) {
        exp_name <- selected_exps[exp_idx]
        exp_folder <- rv$experiment_folders[[exp_name]]

        incProgress(1/(n_exp*2), detail = sprintf("Loading %s", exp_name))

        # Load the appropriate gate strategy for this experiment
        gate_file <- rv$experiment_gate_strategies[[exp_name]]
        if(is.null(gate_file)) gate_file <- "gates_gdef.r"

        gate_path <- file.path("gate_definitions", gate_file)
        if(!file.exists(gate_path)) {
          showNotification(sprintf("Gate file not found: %s. Using default.", gate_file),
                           type = "warning")
          gate_path <- "gate_definitions/gates_gdef.r"
        }

        # Load gates from the selected file
        GATES_env <- new.env()
        tryCatch({
          source(gate_path, local = GATES_env)
          GATES_selected <- GATES_env$GATES
          GATE_STRATEGY_selected <- GATES_env$GATE_STRATEGY
        }, error = function(e) {
          showNotification(sprintf("Error loading %s: %s", gate_file, e$message),
                           type = "error")
          GATES_selected <<- GATES  # Fallback to global GATES
          GATE_STRATEGY_selected <<- NULL
        })

        # Load experiment if not already loaded
        if(is.null(rv$experiments[[exp_name]])) {
          rv$experiments[[exp_name]] <- load_experiment(exp_folder)
        }

        exp <- rv$experiments[[exp_name]]

        # Find control for HA threshold
        control_idx <- find_control_sample(exp$metadata, "Empty_Vector_Dox-")

        if(is.null(control_idx)) {
          showNotification(sprintf("No control found for %s", exp_name),
                           type = "warning")
          next
        }

        control_fcs <- exp$flowset[[control_idx]]
        control_name <- exp$metadata$sample_name[control_idx]
        control_result <- calculate_ha_threshold_from_control(control_fcs, control_name,
                                                               gates = GATES_selected,
                                                               channels = CHANNELS)
        ha_threshold <- control_result$threshold

        # Determine gate ID for this analysis
        gate_id_for_cache <- if(!is.null(GATE_STRATEGY_selected$id)) {
          GATE_STRATEGY_selected$id
        } else {
          "gdef"
        }

        # Check cache first
        cache_data <- load_from_cache(exp_name, GATES_selected, ha_threshold,
                                       gate_id = gate_id_for_cache)

        if(!is.null(cache_data)) {
          # Use cached results
          incProgress(1/(n_exp*2), detail = sprintf("Loaded %s from cache", exp_name))
          rv$ha_thresholds[[exp_name]] <- cache_data$ha_threshold
          exp_results <- cache_data$results

          # Use gate strategy ID (prefer from cache, fallback to GATE_STRATEGY from file)
          gate_id <- if(!is.null(cache_data$gate_id)) {
            cache_data$gate_id
          } else if(!is.null(GATE_STRATEGY_selected$id)) {
            GATE_STRATEGY_selected$id
          } else {
            "gdef"
          }

          # Store gates used for this experiment+strategy combination (prefer from cache, fallback to selected)
          composite_key <- paste0(exp_name, "::", gate_id)
          rv$experiment_gates[[composite_key]] <- if(!is.null(cache_data$gates)) {
            cache_data$gates
          } else {
            GATES_selected
          }

          exp_results$Gate_ID <- gate_id

          n_from_cache <- n_from_cache + 1
        } else {
          # Run analysis
          incProgress(1/(n_exp*2), detail = sprintf("Analyzing %s", exp_name))

          rv$ha_thresholds[[exp_name]] <- ha_threshold
          exp_results <- extract_correlations(exp, ha_threshold,
                                               gates = GATES_selected,
                                               channels = CHANNELS)

          # Save to cache with gate strategy ID
          gate_id <- if(!is.null(GATE_STRATEGY_selected$id)) {
            GATE_STRATEGY_selected$id
          } else {
            "gdef"  # Default fallback
          }

          save_to_cache(exp_name, exp_results, ha_threshold, GATES_selected,
                        gate_strategy_id = gate_id)

          # Store gates used for this experiment+strategy combination
          composite_key <- paste0(exp_name, "::", gate_id)
          rv$experiment_gates[[composite_key]] <- GATES_selected

          # Add Gate_ID to results
          exp_results$Gate_ID <- gate_id

          n_analyzed <- n_analyzed + 1
        }

        all_results[[i]] <- exp_results
      }
      
      # Update sample browser with loaded experiments
      updateSelectInput(session, "selected_experiment", 
                        choices = names(rv$experiments))
      
      # Combine new results
      new_results <- bind_rows(all_results)
      
      # Update existing results table with analyzed data
      for(i in seq_len(nrow(new_results))) {
        # Ensure Gate_ID column exists
        if(!"Gate_ID" %in% names(rv$all_results)) {
          rv$all_results$Gate_ID <- NA_character_
        }

        # Find matching row in all_results (must match Experiment, Well, AND Gate_ID)
        match_idx <- which(rv$all_results$Experiment == new_results$Experiment[i] &
                             rv$all_results$Well == new_results$Well[i] &
                             rv$all_results$Gate_ID == new_results$Gate_ID[i])

        if(length(match_idx) > 0) {
          # Update existing row with this gate strategy
          rv$all_results$Correlation[match_idx] <- new_results$Correlation[i]
          rv$all_results$N_cells[match_idx] <- new_results$N_cells[i]
          rv$all_results$Notes[match_idx] <- new_results$Notes[i]
        } else {
          # This is a new gate strategy for this experiment/well - add a new row
          rv$all_results <- bind_rows(rv$all_results, new_results[i, ])
        }
      }
      
      # Trigger UI refresh to show updated checkmarks/gate IDs
      rv$ui_refresh_trigger <- rv$ui_refresh_trigger + 1
      
      showNotification(sprintf("Analysis complete! %d samples processed (%d from cache, %d newly analyzed).", 
                               nrow(new_results), n_from_cache, n_analyzed), 
                       type = "message", duration = 5)
    })
  })
  
  # Display results table with row selection
  output$results_table <- renderDT({
    req(rv$all_results)
    
    # Format correlation to 4 decimal places
    display_data <- rv$all_results
    if("Correlation" %in% names(display_data)) {
      display_data$Correlation <- ifelse(
        display_data$Correlation == "Not analyzed",
        "Not analyzed",
        sprintf("%.4f", as.numeric(display_data$Correlation))
      )
    }
    
    datatable(display_data, 
              selection = 'single',  # Enable single row selection
              options = list(pageLength = 25, scrollX = TRUE),
              filter = 'top')
  })
  
  # When a row is clicked, load that sample
  observeEvent(input$results_table_rows_selected, {
    req(rv$all_results, input$results_table_rows_selected)
    
    selected_row <- input$results_table_rows_selected
    selected_data <- rv$all_results[selected_row, ]
    
    exp_name <- selected_data$Experiment
    well <- selected_data$Well
    
    # Check if experiment is loaded
    if(!exp_name %in% names(rv$experiments)) {
      showNotification("Please analyze this experiment first to view individual gates", 
                       type = "warning")
      return()
    }
    
    # Find the sample index in the experiment
    exp <- rv$experiments[[exp_name]]
    sample_idx <- which(exp$metadata$well == well)
    
    if(length(sample_idx) == 0) {
      showNotification("Sample not found in loaded experiment", type = "error")
      return()
    }
    
    # Update the selectors to show this sample
    updateSelectInput(session, "selected_experiment", selected = exp_name)
    updateSelectInput(session, "selected_sample", selected = as.character(sample_idx))
    
    showNotification(sprintf("Loaded: %s from %s", selected_data$Sample, exp_name), 
                     type = "message")
  })
  
  # Update overview experiment selector when experiments are loaded
  observe({
    req(rv$experiments)
    updateSelectInput(session, "overview_experiment",
                      choices = names(rv$experiments))
  })

  # Update overview gating strategy selector when experiment changes
  observeEvent(input$overview_experiment, {
    req(input$overview_experiment)
    exp_name <- input$overview_experiment
    available_gates <- rv$experiment_available_gates[[exp_name]]
    if(is.null(available_gates) || length(available_gates) == 0) {
      available_gates <- "gdef"  # Default if none available
    }
    updateSelectInput(session, "overview_gate_strategy",
                      choices = available_gates,
                      selected = available_gates[1])
  })

  # Render overview plot based on selected gate
  output$overview_plot <- renderPlot({
    req(rv$experiments, input$overview_experiment, input$overview_gate, input$overview_gate_strategy)

    exp <- rv$experiments[[input$overview_experiment]]
    exp_name <- input$overview_experiment

    # Use gates for this experiment+strategy combination
    composite_key <- paste0(exp_name, "::", input$overview_gate_strategy)
    gates_to_use <- if(!is.null(rv$experiment_gates[[composite_key]])) {
      rv$experiment_gates[[composite_key]]
    } else {
      GATES
    }

    if(input$overview_gate == "gate1") {
      plot_debris_gate_overview(exp, gates = gates_to_use)

    } else if(input$overview_gate == "gate2") {
      plot_singlet_gate_overview(exp, gates = gates_to_use)

    } else if(input$overview_gate == "gate3") {
      plot_live_gate_overview(exp, gates = gates_to_use)

    } else if(input$overview_gate == "gate4") {
      plot_sphase_outlier_gate_overview(exp, gates = gates_to_use)

    } else if(input$overview_gate == "gate5") {
      plot_fxcycle_quantile_gate_overview(exp, gates = gates_to_use)

    } else if(input$overview_gate == "gate6") {
      plot_edu_fxcycle_gate_overview(exp, gates = gates_to_use)

    } else if(input$overview_gate == "gate7") {
      # Calculate HA threshold
      control_idx <- find_control_sample(exp$metadata, "Empty_Vector_Dox-")
      if(is.null(control_idx)) {
        plot.new()
        text(0.5, 0.5, "No control sample found", cex = 2)
        return()
      }
      control_fcs <- exp$flowset[[control_idx]]
      control_name <- exp$metadata$sample_name[control_idx]
      control_result <- calculate_ha_threshold_from_control(control_fcs, control_name,
                                                             gates = gates_to_use,
                                                             channels = CHANNELS)
      ha_threshold <- control_result$threshold

      plot_ha_gate_overview(exp, ha_threshold, gates = gates_to_use)

    } else if(input$overview_gate == "correlation") {
      # Calculate HA threshold
      control_idx <- find_control_sample(exp$metadata, "Empty_Vector_Dox-")
      if(is.null(control_idx)) {
        plot.new()
        text(0.5, 0.5, "No control sample found", cex = 2)
        return()
      }
      control_fcs <- exp$flowset[[control_idx]]
      control_name <- exp$metadata$sample_name[control_idx]
      control_result <- calculate_ha_threshold_from_control(control_fcs, control_name,
                                                             gates = gates_to_use,
                                                             channels = CHANNELS)
      ha_threshold <- control_result$threshold

      plot_edu_ha_correlation_overview(exp, ha_threshold, gates = gates_to_use)
    }
  })
  
  # Gate plots
  
  # Display current sample name
  output$current_sample_name <- renderText({
    req(rv$experiments, input$selected_experiment, input$selected_sample)
    exp <- rv$experiments[[input$selected_experiment]]
    idx <- as.numeric(input$selected_sample)
    exp$metadata$sample_name[idx]
  })
  
  # Navigation functions for all gates
  navigate_sample <- function(direction) {
    req(rv$experiments, input$selected_experiment, input$selected_sample)
    exp <- rv$experiments[[input$selected_experiment]]
    current_idx <- as.numeric(input$selected_sample)
    n_samples <- length(exp$flowset)
    
    if(direction == "next") {
      new_idx <- if(current_idx < n_samples) current_idx + 1 else 1
    } else {  # previous
      new_idx <- if(current_idx > 1) current_idx - 1 else n_samples
    }
    
    updateSelectInput(session, "selected_sample", selected = as.character(new_idx))
  }
  
  # Gate 1 navigation
  observeEvent(input$gate1_prev, { navigate_sample("prev") })
  observeEvent(input$gate1_next, { navigate_sample("next") })
  
  # Gate 2 navigation
  observeEvent(input$gate2_prev, { navigate_sample("prev") })
  observeEvent(input$gate2_next, { navigate_sample("next") })
  
  # Gate 3 navigation
  observeEvent(input$gate3_prev, { navigate_sample("prev") })
  observeEvent(input$gate3_next, { navigate_sample("next") })
  
  # Gate 4 navigation
  observeEvent(input$gate4_prev, { navigate_sample("prev") })
  observeEvent(input$gate4_next, { navigate_sample("next") })
  
  # Gate 5 navigation
  observeEvent(input$gate5_prev, { navigate_sample("prev") })
  observeEvent(input$gate5_next, { navigate_sample("next") })
  
  # Gate 6 navigation
  observeEvent(input$gate6_prev, { navigate_sample("prev") })
  observeEvent(input$gate6_next, { navigate_sample("next") })
  
  # Gate 7 navigation
  observeEvent(input$gate7_prev, { navigate_sample("prev") })
  observeEvent(input$gate7_next, { navigate_sample("next") })
  
  # Correlation navigation
  observeEvent(input$correlation_prev, { navigate_sample("prev") })
  observeEvent(input$correlation_next, { navigate_sample("next") })
  
  # ==============================================================================
  # GATE EDITING FUNCTIONALITY
  # ==============================================================================
  
  # Store default gates on startup (run once)
  observeEvent(rv$experiments, once = TRUE, {
    isolate({
      if(is.null(rv$gate_storage$default_gates)) {
        rv$gate_storage$default_gates <- GATES
      }
    })
  })
  
  # Enter edit mode (not a toggle - always enter when clicked)
  observeEvent(input$edit_main, {
    # Always enter edit mode when button is clicked
    req(input$selected_experiment, input$selected_sample, input$gate_to_edit)
    
    exp_name <- input$selected_experiment
    sample_idx <- as.numeric(input$selected_sample)
    gate_name <- input$gate_to_edit
    
    # Get default gate
    default_gate <- GATES[[gate_name]]
    
    # Get current gate (custom or default)
    current_gate <- get_gate_for_sample(rv$gate_storage, exp_name, sample_idx, 
                                        gate_name, default_gate)
    
    # Store as temp gate for editing
    rv$temp_gate <- current_gate
    rv$selected_vertex <- NULL
    rv$edit_mode <- TRUE
    
    # Debug notification
    showNotification(sprintf("Edit mode ON - %d vertices loaded", nrow(current_gate)), 
                     type = "message", duration = 2)
  })
  
  # Generate dynamic UI for editing vertex coordinates
  output$vertex_editors_main <- renderUI({
    req(rv$edit_mode, rv$temp_gate)
    
    gate_name <- input$gate_to_edit
    n_vertices <- nrow(rv$temp_gate)
    
    # Check if this is a rectangular gate (Gates 3 and 4)
    is_rectangle <- gate_name %in% c("live_cells", "s_phase_outliers")
    
    cat(sprintf("Rendering vertex editors for %s (rectangle: %s)\n", 
                gate_name, is_rectangle))
    
    if(is_rectangle) {
      # For rectangles, only show bottom-left and top-right corners
      # Assuming standard rectangle: (x1,y1), (x2,y1), (x2,y2), (x1,y2), (x1,y1)
      tagList(
        h5("Rectangle Editing", style = "color: #0066cc;"),
        p("Edit the two diagonal corners. The other corners will be calculated automatically."),
        hr(),
        fluidRow(
          column(12, h6("Bottom-Left Corner", style = "margin-bottom: 5px;")),
          column(6,
                 numericInput("rect_x_min", "X Min:", 
                              value = min(rv$temp_gate[,1]),
                              width = "100%")
          ),
          column(6,
                 numericInput("rect_y_min", "Y Min:", 
                              value = min(rv$temp_gate[,2]),
                              width = "100%")
          )
        ),
        fluidRow(
          column(12, h6("Top-Right Corner", style = "margin-bottom: 5px; margin-top: 10px;")),
          column(6,
                 numericInput("rect_x_max", "X Max:", 
                              value = max(rv$temp_gate[,1]),
                              width = "100%")
          ),
          column(6,
                 numericInput("rect_y_max", "Y Max:", 
                              value = max(rv$temp_gate[,2]),
                              width = "100%")
          )
        )
      )
    } else {
      # For polygons (Gates 1 and 2), show all vertices
      # Don't show last vertex if it's a duplicate of first (closed polygon)
      n_to_show <- if(all(rv$temp_gate[1,] == rv$temp_gate[n_vertices,])) {
        n_vertices - 1
      } else {
        n_vertices
      }
      
      vertex_inputs <- lapply(1:n_to_show, function(i) {
        fluidRow(
          column(12,
                 h6(paste("Vertex", i), style = "margin-bottom: 5px;")
          ),
          column(6,
                 numericInput(paste0("vertex_", i, "_x_main"), "X:", 
                              value = rv$temp_gate[i, 1],
                              width = "100%")
          ),
          column(6,
                 numericInput(paste0("vertex_", i, "_y_main"), "Y:", 
                              value = rv$temp_gate[i, 2],
                              width = "100%")
          )
        )
      })
      
      do.call(tagList, vertex_inputs)
    }
  })
  
  # Update temp_gate when any vertex input changes
  observe({
    req(rv$edit_mode, rv$temp_gate)
    
    gate_name <- input$gate_to_edit
    is_rectangle <- gate_name %in% c("live_cells", "s_phase_outliers")
    
    tryCatch({
      if(is_rectangle) {
        # Handle rectangular gates using corner inputs
        x_min <- input$rect_x_min
        y_min <- input$rect_y_min
        x_max <- input$rect_x_max
        y_max <- input$rect_y_max
        
        if(!is.null(x_min) && !is.null(y_min) && !is.null(x_max) && !is.null(y_max) &&
           !is.na(x_min) && !is.na(y_min) && !is.na(x_max) && !is.na(y_max)) {
          
          # Reconstruct rectangle from corners
          new_gate <- matrix(c(
            x_min, y_min,
            x_max, y_min,
            x_max, y_max,
            x_min, y_max,
            x_min, y_min  # Close the rectangle
          ), ncol = 2, byrow = TRUE)
          
          cat(sprintf("Rectangle updated: X[%.2e, %.2e], Y[%.2e, %.2e]\n", 
                      x_min, x_max, y_min, y_max))
          
          rv$temp_gate <- new_gate
        }
      } else {
        # Handle polygon gates (original code)
        n_vertices <- nrow(rv$temp_gate)
        n_to_show <- if(all(rv$temp_gate[1,] == rv$temp_gate[n_vertices,])) {
          n_vertices - 1
        } else {
          n_vertices
        }
        
        new_gate <- rv$temp_gate
        changes_made <- FALSE
        
        for(i in 1:n_to_show) {
          x_input <- input[[paste0("vertex_", i, "_x_main")]]
          y_input <- input[[paste0("vertex_", i, "_y_main")]]
          
          if(!is.null(x_input) && !is.null(y_input) && 
             !is.na(x_input) && !is.na(y_input)) {
            
            if(!is.numeric(x_input) || !is.numeric(y_input)) {
              cat(sprintf("Warning: Non-numeric input for vertex %d\n", i))
              next
            }
            
            if(is.infinite(x_input) || is.infinite(y_input)) {
              cat(sprintf("Warning: Invalid value (Inf) for vertex %d\n", i))
              next
            }
            
            x_changed <- !isTRUE(all.equal(new_gate[i, 1], x_input))
            y_changed <- !isTRUE(all.equal(new_gate[i, 2], y_input))
            
            if(x_changed || y_changed) {
              cat(sprintf("Vertex %d changed: (%.2e, %.2e) -> (%.2e, %.2e)\n", 
                          i, new_gate[i, 1], new_gate[i, 2], x_input, y_input))
              changes_made <- TRUE
            }
            
            new_gate[i, 1] <- x_input
            new_gate[i, 2] <- y_input
            
            if(i == 1 && n_to_show < n_vertices) {
              new_gate[n_vertices, 1] <- x_input
              new_gate[n_vertices, 2] <- y_input
            }
          }
        }
        
        if(changes_made) {
          cat("Updating rv$temp_gate with new values\n")
        }
        
        rv$temp_gate <- new_gate
      }
      
    }, error = function(e) {
      cat("ERROR in vertex update observer:", e$message, "\n")
      showNotification(paste("Error updating gate:", e$message), 
                       type = "error", duration = 5)
    })
  })
  
  # Remove the old click-based vertex selection (not working reliably)
  # observeEvent(input$gate_plot_click, ...) - REMOVED
  
  # Click on plot to select vertex (for polygon gates only)
  observeEvent(input$gate_plot_click, {
    req(rv$edit_mode, rv$temp_gate)
    
    gate_name <- input$gate_to_edit
    is_polygon <- gate_name %in% c("debris", "singlet")
    
    if(!is_polygon) return()  # Only for polygon gates
    
    click <- input$gate_plot_click
    if(is.null(click)) return()
    
    # Find nearest vertex
    distances <- sqrt((rv$temp_gate[, 1] - click$x)^2 + (rv$temp_gate[, 2] - click$y)^2)
    
    # Get plot range for threshold calculation
    x_range <- diff(par("usr")[1:2])
    y_range <- diff(par("usr")[3:4])
    threshold <- 0.05 * sqrt(x_range^2 + y_range^2)
    
    nearest_idx <- which.min(distances)
    
    # Exclude duplicate last vertex if it exists
    n_vertices <- nrow(rv$temp_gate)
    if(all(rv$temp_gate[1,] == rv$temp_gate[n_vertices,]) && nearest_idx == n_vertices) {
      nearest_idx <- 1
    }
    
    if(distances[nearest_idx] < threshold) {
      rv$selected_vertex <- nearest_idx
      cat(sprintf("Vertex %d selected\n", nearest_idx))
    } else {
      # If a vertex is selected, move it to click location
      if(!is.null(rv$selected_vertex)) {
        rv$temp_gate[rv$selected_vertex, 1] <- click$x
        rv$temp_gate[rv$selected_vertex, 2] <- click$y
        
        # Update duplicate last vertex if needed
        if(rv$selected_vertex == 1 && all(rv$temp_gate[1,] == rv$temp_gate[n_vertices,])) {
          rv$temp_gate[n_vertices, 1] <- click$x
          rv$temp_gate[n_vertices, 2] <- click$y
        }
        
        cat(sprintf("Vertex %d moved to (%.2e, %.2e)\n", rv$selected_vertex, click$x, click$y))
      }
    }
  })
  
  # Deselect vertex button
  observeEvent(input$deselect_vertex, {
    rv$selected_vertex <- NULL
    cat("Vertex deselected\n")
  })
  
  # Apply gate changes
  observeEvent(input$apply_main, {
    req(rv$temp_gate, input$selected_experiment, input$selected_sample, input$gate_to_edit)
    
    cat("=== APPLY BUTTON CLICKED ===\n")
    cat("Gate name:", input$gate_to_edit, "\n")
    cat("Experiment:", input$selected_experiment, "\n")
    cat("Sample:", input$selected_sample, "\n")
    cat("Scope:", input$scope_main, "\n")
    cat("Temp gate dimensions:", nrow(rv$temp_gate), "x", ncol(rv$temp_gate), "\n")
    cat("First vertex:", rv$temp_gate[1,1], ",", rv$temp_gate[1,2], "\n")
    
    # Validate gate
    if(!validate_gate(rv$temp_gate)) {
      showNotification("Invalid gate! Please check coordinates.", type = "error")
      cat("VALIDATION FAILED\n")
      return()
    }
    
    cat("Validation passed, saving gate...\n")
    
    # Save the gate
    rv$gate_storage <- save_custom_gate(
      rv$gate_storage,
      input$selected_experiment,
      as.numeric(input$selected_sample),
      input$gate_to_edit,
      rv$temp_gate,
      input$scope_main
    )
    
    cat("Gate saved successfully\n")
    
    # Exit edit mode
    rv$edit_mode <- FALSE
    rv$temp_gate <- NULL
    rv$selected_vertex <- NULL
    
    showNotification("Gate updated successfully! Click 'Recalculate Results' to apply changes.", 
                     type = "message", duration = 5)
  })
  
  # Cancel editing
  observeEvent(input$cancel_main, {
    rv$edit_mode <- FALSE
    rv$temp_gate <- NULL
    rv$selected_vertex <- NULL
  })
  
  # Reset gate to default
  observeEvent(input$reset_main, {
    req(input$selected_experiment, input$selected_sample, input$gate_to_edit)
    
    rv$gate_storage <- reset_gate_to_default(
      rv$gate_storage,
      input$selected_experiment,
      as.numeric(input$selected_sample),
      input$gate_to_edit
    )
    
    showNotification("Gate reset to default. Click 'Recalculate Results' to apply changes.", 
                     type = "message", duration = 5)
  })
  
  # Recalculate button (placeholder for now)
  observeEvent(input$recalculate_all, {
    showNotification("Recalculation functionality coming in Phase 3!", 
                     type = "message", duration = 3)
  })
  
  # Set default axis limits based on selected gate
  observeEvent(input$gate_to_edit, {
    # Define default limits for each gate (matching the plotting functions)
    defaults <- switch(input$gate_to_edit,
                       "debris" = list(x_min = 0, x_max = 20e6, y_min = 0, y_max = 20e6),
                       "singlet" = list(x_min = 0, x_max = 15e6, y_min = 0, y_max = 4e6),
                       "live_cells" = list(x_min = 100, x_max = 1e6, y_min = 0, y_max = 15e6),
                       "s_phase_outliers" = list(x_min = 0, x_max = 12e6, y_min = 100, y_max = 6e6),
                       list(x_min = NULL, x_max = NULL, y_min = NULL, y_max = NULL)
    )
    
    # Update the input fields
    updateNumericInput(session, "x_min_main", value = defaults$x_min)
    updateNumericInput(session, "x_max_main", value = defaults$x_max)
    updateNumericInput(session, "y_min_main", value = defaults$y_min)
    updateNumericInput(session, "y_max_main", value = defaults$y_max)
  })
  
  # Render gate edit plot
  output$gate_edit_plot <- renderPlot({
    # Require all necessary inputs
    req(rv$experiments)
    req(input$selected_experiment)
    req(input$selected_sample)
    req(input$gate_to_edit)
    
    # Get experiment - return NULL if not found
    exp <- rv$experiments[[input$selected_experiment]]
    req(exp)
    
    idx <- as.numeric(input$selected_sample)
    fcs_data <- exp$flowset[[idx]]
    req(fcs_data)
    
    # Get sample name
    sample_name <- exp$metadata$sample_name[idx]
    
    gate_name <- input$gate_to_edit
    
    # ==============================================================
    # SEQUENTIAL GATING: Apply previous gates first
    # ==============================================================
    current_data <- fcs_data
    
    # Gate order: debris -> singlet -> live_cells -> s_phase_outliers
    gate_sequence <- c("debris", "singlet", "live_cells", "s_phase_outliers")
    current_gate_idx <- which(gate_sequence == gate_name)
    
    cat(sprintf("\n=== GATE EDIT PLOT: %s ===\n", gate_name))
    cat(sprintf("Starting cells: %d\n", nrow(exprs(current_data))))
    
    # Apply all previous gates
    if(current_gate_idx > 1) {
      for(i in 1:(current_gate_idx - 1)) {
        prev_gate_name <- gate_sequence[i]
        
        # Get the gate (custom or default)
        default_prev_gate <- GATES[[prev_gate_name]]
        prev_gate <- isolate({
          get_gate_for_sample(rv$gate_storage, input$selected_experiment, 
                              idx, prev_gate_name, default_prev_gate)
        })
        
        # Determine channels for this previous gate
        prev_channels <- switch(prev_gate_name,
                                "debris" = c("FSC-A", "SSC-A"),
                                "singlet" = c("FSC-A", "FSC-H"),
                                "live_cells" = c("FL5-A", "SSC-A"),
                                "s_phase_outliers" = c("FL6-A", "FL1-A")
        )
        
        # Apply the gate filter
        prev_filter <- point.in.polygon(
          exprs(current_data)[, prev_channels[1]],
          exprs(current_data)[, prev_channels[2]],
          prev_gate[, 1],
          prev_gate[, 2]
        ) > 0
        
        current_data <- Subset(current_data, prev_filter)
        cat(sprintf("After gate %d (%s): %d cells remaining\n", 
                    i, prev_gate_name, nrow(exprs(current_data))))
      }
    }
    
    # Now current_data contains only cells that passed all previous gates
    # ==============================================================
    
    cat(sprintf("Final cells for %s gate: %d\n", gate_name, nrow(exprs(current_data))))
    
    # Check if we have any cells left
    if(nrow(exprs(current_data)) == 0) {
      # No cells remaining - show empty plot with error message
      plot.new()
      text(0.5, 0.5, 
           paste0("No cells remaining after previous gates!\n\n",
                  "Previous gates may be filtering out all cells.\n",
                  "Try editing earlier gates or click 'Reset to Default'."),
           cex = 1.2, col = "red")
      return()
    }
    
    # Get current gate (custom or default) - isolate to prevent reactivity loops
    default_gate <- GATES[[gate_name]]
    current_gate <- isolate({
      get_gate_for_sample(rv$gate_storage, input$selected_experiment, 
                          idx, gate_name, default_gate)
    })
    
    # Determine channels and labels based on gate
    channel_info <- switch(gate_name,
                           "debris" = list(
                             x_channel = "FSC-A", y_channel = "SSC-A",
                             x_label = "FSC-A", y_label = "SSC-A",
                             log_scale = ""
                           ),
                           "singlet" = list(
                             x_channel = "FSC-A", y_channel = "FSC-H",
                             x_label = "FSC-A", y_label = "FSC-H",
                             log_scale = ""
                           ),
                           "live_cells" = list(
                             x_channel = "FL5-A", y_channel = "SSC-A",
                             x_label = "DCM-A", y_label = "SSC-A",
                             log_scale = "x"
                           ),
                           "s_phase_outliers" = list(
                             x_channel = "FL6-A", y_channel = "FL1-A",
                             x_label = "FxCycle-A", y_label = "EdU-A",
                             log_scale = "y"
                           )
    )
    
    # Show data range for debugging
    x_data_all <- exprs(current_data)[, channel_info$x_channel]
    y_data_all <- exprs(current_data)[, channel_info$y_channel]
    cat(sprintf("Data range for %s:\n", gate_name))
    cat(sprintf("  X (%s): %.2e to %.2e\n", 
                channel_info$x_label, min(x_data_all), max(x_data_all)))
    cat(sprintf("  Y (%s): %.2e to %.2e\n", 
                channel_info$y_label, min(y_data_all), max(y_data_all)))
    cat(sprintf("Gate coordinates:\n"))
    cat(sprintf("  X range: %.2e to %.2e\n", 
                min(current_gate[,1]), max(current_gate[,1])))
    cat(sprintf("  Y range: %.2e to %.2e\n", 
                min(current_gate[,2]), max(current_gate[,2])))
    
    
    # Get axis limits from inputs (NULL if not set)
    xlim <- if(!is.null(input$x_min_main) && !is.null(input$x_max_main)) {
      c(input$x_min_main, input$x_max_main)
    } else {
      NULL
    }
    
    ylim <- if(!is.null(input$y_min_main) && !is.null(input$y_max_main)) {
      c(input$y_min_main, input$y_max_main)
    } else {
      NULL
    }
    
    # Plot with interactive editing if in edit mode
    if(rv$edit_mode && !is.null(rv$temp_gate)) {
      plot_gate_interactive(current_data, current_gate, 
                            channel_info$x_channel, channel_info$y_channel,
                            gate_name, sample_name, edit_mode = TRUE, 
                            selected_vertex = rv$selected_vertex,
                            temp_gate = rv$temp_gate,
                            xlim = xlim, ylim = ylim,
                            x_label = channel_info$x_label,
                            y_label = channel_info$y_label,
                            log_scale = channel_info$log_scale)
    } else {
      plot_gate_interactive(current_data, current_gate, 
                            channel_info$x_channel, channel_info$y_channel,
                            gate_name, sample_name, edit_mode = FALSE,
                            xlim = xlim, ylim = ylim,
                            x_label = channel_info$x_label,
                            y_label = channel_info$y_label,
                            log_scale = channel_info$log_scale)
    }
  })
  
  # Display current gate coordinates
  output$coords_main <- renderText({
    if(rv$edit_mode && !is.null(rv$temp_gate)) {
      format_gate_coordinates(rv$temp_gate)
    } else if(!is.null(input$selected_experiment) && !is.null(input$selected_sample) && 
              !is.null(input$gate_to_edit)) {
      default_gate <- GATES[[input$gate_to_edit]]
      current_gate <- isolate({
        get_gate_for_sample(rv$gate_storage, input$selected_experiment,
                            as.numeric(input$selected_sample),
                            input$gate_to_edit, default_gate)
      })
      format_gate_coordinates(current_gate)
    } else {
      "Select a sample to view coordinates"
    }
  })
  
  # ==============================================================================
  # EXISTING GATE PLOT OUTPUTS
  # ==============================================================================
  
  output$gate1_plot <- renderPlot({
    req(rv$experiments, input$selected_experiment, input$selected_sample, input$browse_gate_strategy)
    exp <- rv$experiments[[input$selected_experiment]]
    exp_name <- input$selected_experiment
    idx <- as.numeric(input$selected_sample)

    # Use gates for this experiment+strategy combination
    composite_key <- paste0(exp_name, "::", input$browse_gate_strategy)
    gates_to_use <- if(!is.null(rv$experiment_gates[[composite_key]])) {
      rv$experiment_gates[[composite_key]]
    } else {
      GATES
    }

    plot_debris_gate_single(exp$flowset[[idx]], exp$metadata$sample_name[idx], gates = gates_to_use)
  })

  output$gate2_plot <- renderPlot({
    req(rv$experiments, input$selected_experiment, input$selected_sample, input$browse_gate_strategy)
    exp <- rv$experiments[[input$selected_experiment]]
    exp_name <- input$selected_experiment
    idx <- as.numeric(input$selected_sample)

    # Use gates for this experiment+strategy combination
    composite_key <- paste0(exp_name, "::", input$browse_gate_strategy)
    gates_to_use <- if(!is.null(rv$experiment_gates[[composite_key]])) {
      rv$experiment_gates[[composite_key]]
    } else {
      GATES
    }

    plot_singlet_gate_single(exp$flowset[[idx]], exp$metadata$sample_name[idx], gates = gates_to_use)
  })

  output$gate3_plot <- renderPlot({
    req(rv$experiments, input$selected_experiment, input$selected_sample, input$browse_gate_strategy)
    exp <- rv$experiments[[input$selected_experiment]]
    exp_name <- input$selected_experiment
    idx <- as.numeric(input$selected_sample)

    # Use gates for this experiment+strategy combination
    composite_key <- paste0(exp_name, "::", input$browse_gate_strategy)
    gates_to_use <- if(!is.null(rv$experiment_gates[[composite_key]])) {
      rv$experiment_gates[[composite_key]]
    } else {
      GATES
    }

    plot_live_gate_single(exp$flowset[[idx]], exp$metadata$sample_name[idx], gates = gates_to_use)
  })

  output$gate4_plot <- renderPlot({
    req(rv$experiments, input$selected_experiment, input$selected_sample, input$browse_gate_strategy)
    exp <- rv$experiments[[input$selected_experiment]]
    exp_name <- input$selected_experiment
    idx <- as.numeric(input$selected_sample)

    # Use gates for this experiment+strategy combination
    composite_key <- paste0(exp_name, "::", input$browse_gate_strategy)
    gates_to_use <- if(!is.null(rv$experiment_gates[[composite_key]])) {
      rv$experiment_gates[[composite_key]]
    } else {
      GATES
    }

    plot_sphase_outlier_gate_single(exp$flowset[[idx]], exp$metadata$sample_name[idx], gates = gates_to_use)
  })
  
  output$gate5_plot <- renderPlot({
    req(rv$experiments, input$selected_experiment, input$selected_sample, input$browse_gate_strategy)
    exp <- rv$experiments[[input$selected_experiment]]
    exp_name <- input$selected_experiment
    idx <- as.numeric(input$selected_sample)

    # Use gates for this experiment+strategy combination
    composite_key <- paste0(exp_name, "::", input$browse_gate_strategy)
    gates_to_use <- if(!is.null(rv$experiment_gates[[composite_key]])) {
      rv$experiment_gates[[composite_key]]
    } else {
      GATES
    }

    plot_fxcycle_quantile_gate_single(exp$flowset[[idx]], exp$metadata$sample_name[idx], gates = gates_to_use)
  })

  output$gate6_plot <- renderPlot({
    req(rv$experiments, input$selected_experiment, input$selected_sample, input$browse_gate_strategy)
    exp <- rv$experiments[[input$selected_experiment]]
    exp_name <- input$selected_experiment
    idx <- as.numeric(input$selected_sample)

    # Use gates for this experiment+strategy combination
    composite_key <- paste0(exp_name, "::", input$browse_gate_strategy)
    gates_to_use <- if(!is.null(rv$experiment_gates[[composite_key]])) {
      rv$experiment_gates[[composite_key]]
    } else {
      GATES
    }

    plot_edu_fxcycle_gate_single(exp$flowset[[idx]], exp$metadata$sample_name[idx], gates = gates_to_use)
  })
  
  output$gate7_plot <- renderPlot({
    req(rv$experiments, input$selected_experiment, input$selected_sample, input$browse_gate_strategy)
    exp <- rv$experiments[[input$selected_experiment]]
    exp_name <- input$selected_experiment
    idx <- as.numeric(input$selected_sample)

    # Use gates for this experiment+strategy combination
    composite_key <- paste0(exp_name, "::", input$browse_gate_strategy)
    gates_to_use <- if(!is.null(rv$experiment_gates[[composite_key]])) {
      rv$experiment_gates[[composite_key]]
    } else {
      GATES
    }

    # Calculate HA threshold for this experiment
    control_idx <- find_control_sample(exp$metadata, "Empty_Vector_Dox-")
    if(is.null(control_idx)) {
      plot.new()
      text(0.5, 0.5, "No control sample found", cex = 1.5)
      return()
    }

    control_fcs <- exp$flowset[[control_idx]]
    control_name <- exp$metadata$sample_name[control_idx]
    control_result <- calculate_ha_threshold_from_control(control_fcs, control_name,
                                                           gates = gates_to_use,
                                                           channels = CHANNELS)
    ha_threshold <- control_result$threshold

    plot_ha_gate_single(exp$flowset[[idx]], exp$metadata$sample_name[idx], ha_threshold,
                        gates = gates_to_use)
  })
  
  output$correlation_plot <- renderPlot({
    req(rv$experiments, input$selected_experiment, input$selected_sample, input$browse_gate_strategy)
    exp <- rv$experiments[[input$selected_experiment]]
    exp_name <- input$selected_experiment
    idx <- as.numeric(input$selected_sample)

    # Use gates for this experiment+strategy combination
    composite_key <- paste0(exp_name, "::", input$browse_gate_strategy)
    gates_to_use <- if(!is.null(rv$experiment_gates[[composite_key]])) {
      rv$experiment_gates[[composite_key]]
    } else {
      GATES
    }

    # Calculate HA threshold for this experiment
    control_idx <- find_control_sample(exp$metadata, "Empty_Vector_Dox-")
    if(is.null(control_idx)) {
      plot.new()
      text(0.5, 0.5, "No control sample found", cex = 1.5)
      return()
    }

    control_fcs <- exp$flowset[[control_idx]]
    control_name <- exp$metadata$sample_name[control_idx]
    control_result <- calculate_ha_threshold_from_control(control_fcs, control_name,
                                                           gates = gates_to_use,
                                                           channels = CHANNELS)
    ha_threshold <- control_result$threshold

    plot_edu_ha_correlation_single(exp$flowset[[idx]], exp$metadata$sample_name[idx], ha_threshold,
                                     gates = gates_to_use)
  })
  
  # Status text
  output$status_text <- renderText({
    status <- "Ready"
    if(!is.null(rv$experiments)) {
      status <- sprintf("%s\nLoaded %d experiments", status, length(rv$experiments))
    }
    if(!is.null(rv$all_results)) {
      status <- sprintf("%s\nAnalyzed %d samples", status, nrow(rv$all_results))
    }
    status
  })
  
  # Download handler
  output$download_results <- downloadHandler(
    filename = function() {
      paste0("flow_analysis_results_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".xlsx")
    },
    content = function(file) {
      req(rv$all_results)
      write.xlsx(rv$all_results, file)
    }
  )
  
  # Multi-sample comparison
  
  # Create reactive value to store selected samples for comparison
  rv$comparison_samples <- reactiveVal(data.frame())
  
  # Display sample selector table (only analyzed samples)
  output$comparison_sample_selector <- renderDT({
    req(rv$all_results)
    
    # Filter to only analyzed samples (those with numeric correlation)
    analyzed <- rv$all_results[rv$all_results$Correlation != "Not analyzed", ]
    
    if(nrow(analyzed) == 0) {
      return(data.frame(Message = "No analyzed samples available. Please analyze experiments first."))
    }
    
    # Show relevant columns
    cols_to_show <- c("Experiment", "Sample", "Cell_line", "Gene",
                      "Mutation", "Correlation", "N_cells", "Gate_ID")
    # Add Notes if it exists
    if("Notes" %in% names(analyzed)) {
      cols_to_show <- c(cols_to_show, "Notes")
    }
    display_data <- analyzed[, cols_to_show]

    # Format correlation to 4 decimal places
    display_data$Correlation <- sprintf("%.4f", as.numeric(display_data$Correlation))

    # Ensure all columns are character/factor for proper filtering
    display_data$Gate_ID <- as.character(display_data$Gate_ID)
    display_data$Experiment <- as.character(display_data$Experiment)
    if("Notes" %in% names(display_data)) {
      display_data$Notes <- as.character(display_data$Notes)
    }

    datatable(display_data,
              selection = 'multiple',
              options = list(
                pageLength = 10,
                scrollX = TRUE,
                search = list(regex = FALSE, caseInsensitive = TRUE)
              ),
              filter = list(position = 'top', clear = FALSE))
  })
  
  # Store previous selection to detect additions vs removals
  rv$prev_selected_rows <- reactiveVal(integer(0))
  
  # Auto-select samples from same cell line (but allow deselection)
  observeEvent(input$comparison_sample_selector_rows_selected, {
    req(rv$all_results)
    
    analyzed <- rv$all_results[rv$all_results$Correlation != "Not analyzed", ]
    selected_rows <- input$comparison_sample_selector_rows_selected
    prev_rows <- rv$prev_selected_rows()
    
    if(length(selected_rows) == 0) {
      rv$prev_selected_rows(integer(0))
      return()
    }
    
    # Determine if user is adding or removing selections
    newly_added <- setdiff(selected_rows, prev_rows)
    
    # Only auto-select if rows were added (not removed)
    if(length(newly_added) > 0) {
      # Get cell lines of newly added samples
      new_cell_lines <- unique(analyzed$Cell_line[newly_added])
      
      # Find all rows with matching cell lines
      matching_rows <- which(analyzed$Cell_line %in% new_cell_lines)
      
      # Combine with existing selection
      all_selected <- unique(c(selected_rows, matching_rows))
      
      # Update selection if different
      if(!identical(sort(all_selected), sort(selected_rows))) {
        dataTableProxy('comparison_sample_selector') %>% 
          selectRows(all_selected)
        rv$prev_selected_rows(all_selected)
      } else {
        rv$prev_selected_rows(selected_rows)
      }
    } else {
      # User is deselecting - just update the stored selection
      rv$prev_selected_rows(selected_rows)
    }
  })
  
  # Show selected samples
  output$selected_samples_list <- renderText({
    req(input$comparison_sample_selector_rows_selected)
    
    analyzed <- rv$all_results[rv$all_results$Correlation != "Not analyzed", ]
    selected_rows <- input$comparison_sample_selector_rows_selected
    selected_data <- analyzed[selected_rows, ]
    
    # Group by cell line
    grouped <- selected_data %>%
      group_by(Cell_line, Mutation) %>%
      summarize(n = n(), .groups = 'drop') %>%
      mutate(
        Cell_line = as.character(Cell_line),
        Mutation = as.character(Mutation)
      ) %>%
      arrange(Cell_line)
    
    paste(sprintf("%s (%s): %d replicates",
                  grouped$Cell_line,
                  grouped$Mutation,
                  grouped$n),
          collapse = "\n")
  })
  
  # Clear selection
  observeEvent(input$clear_selection, {
    dataTableProxy('comparison_sample_selector') %>% 
      selectRows(NULL)
  })
  
  # Generate comparison plot
  observeEvent(input$plot_comparison, {
    req(input$comparison_sample_selector_rows_selected)
    
    analyzed <- rv$all_results[rv$all_results$Correlation != "Not analyzed", ]
    selected_rows <- input$comparison_sample_selector_rows_selected
    selected_data <- analyzed[selected_rows, ]
    
    # Store for plotting
    rv$comparison_samples(selected_data)
    
    # Update reference group choices
    unique_lines <- unique(selected_data$Cell_line)
    # Properly match each cell line to its mutation
    choices_labels <- sapply(unique_lines, function(cl) {
      mut <- selected_data$Mutation[selected_data$Cell_line == cl][1]
      paste0(mut, " (#", cl, ")")
    })
    choices <- setNames(unique_lines, choices_labels)
    updateSelectInput(session, "reference_group", 
                      choices = choices,
                      selected = unique_lines[1])
  })
  
  
  # Render comparison plot
  output$comparison_plot <- renderPlot({
    req(rv$comparison_samples())
    
    plot_data <- rv$comparison_samples()
    
    if(nrow(plot_data) == 0) {
      plot.new()
      text(0.5, 0.5, "Select samples and click 'Generate Comparison Plot'", cex = 1.5)
      return()
    }
    
    # Convert correlation to numeric
    plot_data$Correlation <- as.numeric(plot_data$Correlation)
    
    # Group by cell line to get means
    plot_summary <- plot_data %>%
      group_by(Cell_line, Mutation, Gene) %>%
      summarize(
        mean_corr = mean(Correlation, na.rm = TRUE),
        sd_corr = ifelse(n() > 1, sd(Correlation, na.rm = TRUE), 0),  # Set SD to 0 for single replicates
        n = n(),
        .groups = 'drop'
      ) %>%
      mutate(
        Cell_line = as.character(Cell_line),
        Mutation = as.character(Mutation),
        Gene = as.character(Gene)
      )
    
    # Sort: use custom order if available, otherwise WT first then by correlation
    if(!is.null(rv$sample_order())) {
      # Apply custom order
      plot_summary$Cell_line <- factor(plot_summary$Cell_line, levels = rv$sample_order())
      plot_summary <- plot_summary %>% arrange(Cell_line)
    } else {
      # Default: WT first, others by correlation
      wt_rows <- grep("^WT", plot_summary$Mutation, ignore.case = TRUE)
      if(length(wt_rows) > 0) {
        wt_data <- plot_summary[wt_rows, ]
        other_data <- plot_summary[-wt_rows, ] %>% arrange(mean_corr)
        plot_summary <- bind_rows(wt_data, other_data)
      } else {
        plot_summary <- plot_summary %>% arrange(mean_corr)
      }
    }
    
    # Create labels - just mutation name (no cell line number)
    plot_summary$label <- plot_summary$Mutation
    
    # Assign lighter colors to experiments for points
    exp_colors <- rainbow(length(unique(plot_data$Experiment)), alpha = 0.5)
    names(exp_colors) <- unique(plot_data$Experiment)
    
    # Create the plot with adaptive top margin based on number of experiments
    n_experiments <- length(unique(plot_data$Experiment))
    top_margin <- 8 + ceiling(n_experiments * 0.7)  # Scales with number of experiments
    par(mar = c(8, 4, top_margin, 4))
    
    # Fixed narrow bars, add empty space on right when few samples
    n_samples <- nrow(plot_summary)
    bar_width <- 0.8  # Fixed narrow width
    bar_space <- 1  # Fixed moderate spacing
    
    # Calculate total x-range: expand to fixed width regardless of sample count
    bars_width <- n_samples * (bar_width + bar_space)
    x_max <- max(bars_width, 15)  # Always extend to at least 15 units
    
    bp <- barplot(plot_summary$mean_corr,
                  names.arg = plot_summary$label,
                  las = 2,
                  xlim = c(0, x_max),  # Fixed x-range
                  ylim = c(-0.8, 
                           max(c(0, plot_data$Correlation, plot_summary$mean_corr)) + 0.3),
                  ylab = "EdU vs HA Correlation (r)",
                  main = "",
                  col = "lightgrey",
                  border = "black",
                  cex.names = 0.8,
                  width = bar_width,
                  space = bar_space)
    
    # Add title with space for legend
    title(main = "Multi-Sample Correlation Comparison", line = 3.5)
    
    # Add error bars (SD) - only for samples with SD > 0
    for(i in seq_len(nrow(plot_summary))) {
      if(!is.na(plot_summary$sd_corr[i]) && plot_summary$sd_corr[i] > 0) {
        arrows(bp[i], plot_summary$mean_corr[i] - plot_summary$sd_corr[i],
               bp[i], plot_summary$mean_corr[i] + plot_summary$sd_corr[i],
               angle = 90, code = 3, length = 0.1)
      }
    }
    
    # Add horizontal line at 0
    abline(h = 0, lty = 2, col = "gray40")
    
    # Add individual data points
    for(i in seq_len(nrow(plot_summary))) {
      cell_line <- plot_summary$Cell_line[i]
      cell_data <- plot_data[plot_data$Cell_line == cell_line, ]
      
      # Add jitter to x position for visibility
      x_pos <- bp[i] + runif(nrow(cell_data), -0.15, 0.15)
      
      for(j in seq_len(nrow(cell_data))) {
        exp_name <- cell_data$Experiment[j]
        points(x_pos[j], cell_data$Correlation[j], 
               pch = 21, cex = 2, 
               bg = exp_colors[exp_name], 
               col = "black", lwd = 1.5)
      }
    }
    
    # Add statistical significance using Holm-Sidak correction (all at same height)
    if(nrow(plot_summary) > 1) {
      # Prepare data for matched analysis
      # If there are multiple replicates per experiment/cell_line, take the mean
      wide_data <- plot_data %>%
        select(Experiment, Cell_line, Correlation) %>%
        pivot_wider(names_from = Cell_line, values_from = Correlation, values_fn = mean)
      
      long_data <- wide_data %>%
        pivot_longer(cols = -Experiment, names_to = "Cell_line", values_to = "Correlation") %>%
        filter(!is.na(Correlation))
      
      # Reference group - use selected reference or default to first
      ref_group <- if(!is.null(input$reference_group) && input$reference_group != "") {
        input$reference_group
      } else {
        plot_summary$Cell_line[1]
      }
      ref_data <- long_data %>% filter(Cell_line == ref_group)
      
      # Perform pairwise t-tests
      p_values <- c()
      test_groups <- c()
      
      for(i in 1:nrow(plot_summary)) {
        test_group <- plot_summary$Cell_line[i]

        # Skip the reference group itself
        if(test_group == ref_group) next

        test_data <- long_data %>% filter(Cell_line == test_group)

        # Match by experiment - only keep experiments that have both samples
        merged <- merge(ref_data, test_data, by = "Experiment", suffixes = c("_ref", "_test"))

        # Remove any rows with NA correlations
        merged <- merged %>% filter(!is.na(Correlation_ref) & !is.na(Correlation_test))

        if(nrow(merged) >= 2) {
          # Paired t-test (only for experiments with both samples)
          tryCatch({
            t_result <- t.test(merged$Correlation_ref, merged$Correlation_test, paired = TRUE)
            p_values <- c(p_values, t_result$p.value)
            test_groups <- c(test_groups, test_group)
          }, error = function(e) {
            # If t-test fails, mark as NA
            p_values <<- c(p_values, NA)
            test_groups <<- c(test_groups, test_group)
          })
        } else {
          # Not enough matched pairs for comparison
          p_values <- c(p_values, NA)
          test_groups <- c(test_groups, test_group)
        }
      }
      
      # Apply Holm-Sidak correction
      p_adjusted <- p.adjust(p_values, method = "holm")
      
      # Find the maximum y position for all bars
      max_y <- max(c(plot_summary$mean_corr + plot_summary$sd_corr, 
                     plot_data$Correlation))
      sig_y_pos <- max_y + 0.12
      
      # Display results
      for(i in seq_along(test_groups)) {
        plot_idx <- which(plot_summary$Cell_line == test_groups[i])
        
        # Only display if we have a valid p-value (skip n/a cases)
        if(!is.na(p_adjusted[i])) {
          stars <- if(p_adjusted[i] < 0.001) "***"
          else if(p_adjusted[i] < 0.01) "**"
          else if(p_adjusted[i] < 0.05) "*"
          else "ns"
          
          # Format p-value
          if(p_adjusted[i] < 0.001) {
            p_text <- sprintf("(p<0.001)")
          } else {
            p_text <- sprintf("(p=%.3f)", p_adjusted[i])
          }
          
          # Show stars/ns at consistent height
          text(bp[plot_idx], sig_y_pos, stars, cex = 1.2, font = 2, col = "black")
          
          # Show p-value below stars
          text(bp[plot_idx], sig_y_pos - 0.06, p_text, cex = 0.7, col = "black")
        }
      }
    }
    
    # Add legend above the plot (anchored to right edge)
    par(xpd = TRUE)
    legend(x = par("usr")[2], y = par("usr")[4] + 0.25,  # Right edge of plot area
           legend = names(exp_colors),
           pch = 21, pt.bg = exp_colors, pt.cex = 1.3, col = "black",
           title = "Experiments",
           ncol = 1,
           cex = 0.7,
           bg = "white",
           box.lty = 1,
           xjust = 1,
           yjust = 0)
    par(xpd = FALSE)
    
  })
  
  
  # Store custom order
  rv$sample_order <- reactiveVal(NULL)
  
  # Reorder samples with drag-and-drop
  observeEvent(input$reorder_samples, {
    req(rv$comparison_samples())
    
    plot_data <- rv$comparison_samples()
    plot_data$Correlation <- as.numeric(plot_data$Correlation)
    
    # Get current order
    current_order <- plot_data %>%
      group_by(Cell_line, Mutation) %>%
      summarize(mean_corr = mean(Correlation), .groups = 'drop') %>%
      arrange(mean_corr)
    
    # WT first
    wt_rows <- grep("^WT", current_order$Mutation, ignore.case = TRUE)
    if(length(wt_rows) > 0) {
      wt_data <- current_order[wt_rows, ]
      other_data <- current_order[-wt_rows, ]
      current_order <- bind_rows(wt_data, other_data)
    }
    
    # Create labels for sortable list
    sample_labels <- setNames(
      current_order$Mutation,
      current_order$Cell_line
    )
    
    showModal(modalDialog(
      title = "Reorder Samples (Drag to Reorder)",
      size = "m",
      p(strong("Drag samples to reorder. WT will remain first.")),
      rank_list(
        text = "Sample Order:",
        labels = sample_labels,
        input_id = "sample_rank_list"
      ),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("apply_reorder", "Apply Order", class = "btn-primary")
      )
    ))
  })
  
  # Apply custom order from drag-and-drop
  observeEvent(input$apply_reorder, {
    req(input$sample_rank_list)
    
    # The rank_list returns cell_line values in the new order
    new_order <- input$sample_rank_list
    
    rv$sample_order(new_order)
    
    removeModal()
    showNotification("Custom order applied! Click 'Generate Comparison Plot' to update.", 
                     type = "message")
  })
  
  # Apply custom order
  observeEvent(input$apply_reorder, {
    req(input$reorder_sequence, rv$comparison_samples())
    
    plot_data <- rv$comparison_samples()
    plot_data$Correlation <- as.numeric(plot_data$Correlation)
    
    current_order <- plot_data %>%
      group_by(Cell_line, Mutation) %>%
      summarize(mean_corr = mean(Correlation), .groups = 'drop') %>%
      arrange(mean_corr)
    
    # WT first
    wt_rows <- grep("^WT", current_order$Mutation, ignore.case = TRUE)
    if(length(wt_rows) > 0) {
      wt_data <- current_order[wt_rows, ]
      other_data <- current_order[-wt_rows, ]
      current_order <- bind_rows(wt_data, other_data)
    }
    
    # Parse the sequence
    tryCatch({
      sequence <- as.integer(strsplit(input$reorder_sequence, ",")[[1]])
      
      if(length(sequence) != nrow(current_order) || 
         !all(sort(sequence) == seq_len(nrow(current_order)))) {
        showNotification("Invalid sequence! Must include all numbers from 1 to the number of samples.", 
                         type = "error")
        return()
      }
      
      # Apply the new order
      new_order <- current_order$Cell_line[sequence]
      rv$sample_order(new_order)
      
      removeModal()
      showNotification("Custom order applied! Click 'Generate Comparison Plot' to update.", 
                       type = "message")
      
    }, error = function(e) {
      showNotification("Error parsing order sequence. Use format: 1,2,3,4", type = "error")
    })
  })
  
  # Statistical analysis (matching GraphPad approach)
  output$stats_output <- renderText({
    req(rv$comparison_samples())
    
    plot_data <- rv$comparison_samples()
    
    if(nrow(plot_data) < 2) {
      return("Select at least 2 samples for statistical comparison")
    }
    
    plot_data$Correlation <- as.numeric(plot_data$Correlation)
    
    # Prepare grouped data (for reference, not displayed)
    grouped <- plot_data %>%
      group_by(Cell_line, Mutation) %>%
      summarize(
        mean_corr = mean(Correlation, na.rm = TRUE),
        sd_corr = ifelse(n() > 1, sd(Correlation, na.rm = TRUE), 0),
        n = n(),
        .groups = 'drop'
      ) %>%
      mutate(
        Cell_line = as.character(Cell_line),
        Mutation = as.character(Mutation)
      )
    
    stats_text <- ""
    
    # Repeated measures one-way ANOVA
    if(length(unique(plot_data$Cell_line)) > 1) {
      
      # Check if we have matched data (same experiments across groups)
      experiments_per_group <- plot_data %>%
        group_by(Cell_line) %>%
        summarize(exps = list(unique(Experiment)), .groups = 'drop')
      
      # Prepare data for RM-ANOVA
      # Need balanced design with same experiments
      # If there are multiple replicates per experiment/cell_line, take the mean
      wide_data <- plot_data %>%
        select(Experiment, Cell_line, Correlation) %>%
        pivot_wider(names_from = Cell_line, values_from = Correlation, values_fn = mean)
      
      # Check if data is suitable for RM-ANOVA
      n_complete <- sum(complete.cases(wide_data))
      
      if(n_complete >= 2) {
        # Perform repeated measures ANOVA
        stats_text <- paste0(stats_text, "\n=== Repeated Measures One-Way ANOVA ===\n")
        stats_text <- paste0(stats_text, "(Samples matched by experiment, assuming Gaussian distribution)\n\n")
        
        # Convert back to long format for aov
        long_data <- wide_data %>%
          pivot_longer(cols = -Experiment, names_to = "Cell_line", values_to = "Correlation") %>%
          filter(!is.na(Correlation))
        
        # Perform RM-ANOVA
        aov_model <- aov(Correlation ~ Cell_line + Error(Experiment/Cell_line), data = long_data)
        aov_summary <- summary(aov_model)
        
        # Extract F-statistic and p-value
        # The within-subjects effect is in the second element
        within_effect <- aov_summary$`Error: Experiment:Cell_line`[[1]]
        
        if(!is.null(within_effect) && nrow(within_effect) > 0) {
          f_stat <- within_effect["Cell_line", "F value"]
          p_val <- within_effect["Cell_line", "Pr(>F)"]
          df1 <- within_effect["Cell_line", "Df"]
          df2 <- within_effect["Residuals", "Df"]
          
          stats_text <- paste0(stats_text, sprintf("F(%d, %d) = %.3f, p = %.4f\n", 
                                                   df1, df2, f_stat, p_val))
          
          if(p_val < 0.05) {
            stats_text <- paste0(stats_text, "Result: Significant difference detected (p < 0.05)\n\n")
          } else {
            stats_text <- paste0(stats_text, "Result: No significant difference (p >= 0.05)\n\n")
          }
          
          # Always show pairwise comparisons (regardless of ANOVA result)
          stats_text <- paste0(stats_text, "\n=== Pairwise Comparisons (Holm-Sidak) ===\n")
          
          # Get reference group - use selected reference or default to first
          ref_group <- if(!is.null(input$reference_group) && input$reference_group != "") {
            input$reference_group
          } else {
            grouped %>% slice(1) %>% pull(Cell_line)
          }

          ref_mut <- grouped %>% filter(Cell_line == ref_group) %>% slice(1) %>% pull(Mutation)
          stats_text <- paste0(stats_text, sprintf("(Compared to %s #%s)\n\n", as.character(ref_mut), as.character(ref_group)))
          
          ref_data <- long_data %>% filter(Cell_line == ref_group)
          
          # Perform pairwise t-tests
          p_values <- c()
          test_info <- list()
          
          for(i in 1:nrow(grouped)) {
            test_group <- grouped %>% slice(i) %>% pull(Cell_line)

            # Skip if this is the reference group
            if(test_group == ref_group) next

            test_data <- long_data %>% filter(Cell_line == test_group)
            test_mut <- grouped %>% slice(i) %>% pull(Mutation)

            # Match by experiment - only keep experiments with both samples
            merged <- merge(ref_data, test_data, by = "Experiment", suffixes = c("_ref", "_test"))

            # Remove any rows with NA correlations
            merged <- merged %>% filter(!is.na(Correlation_ref) & !is.na(Correlation_test))

            if(nrow(merged) >= 2) {
              # Paired t-test
              tryCatch({
                t_result <- t.test(merged$Correlation_ref, merged$Correlation_test, paired = TRUE)
                p_values <- c(p_values, t_result$p.value)
                test_info[[length(test_info) + 1]] <- list(
                  test_group = as.character(test_group),
                  test_mut = as.character(test_mut),
                  n_pairs = nrow(merged)
                )
              }, error = function(e) {
                # Skip if t-test fails
              })
            }
          }
          
          # Apply Holm-Sidak correction
          if(length(p_values) > 0) {
            p_adjusted <- p.adjust(p_values, method = "holm")
            
            for(i in seq_along(test_info)) {
              sig <- if(p_adjusted[i] < 0.001) "***"
              else if(p_adjusted[i] < 0.01) "**"
              else if(p_adjusted[i] < 0.05) "*"
              else "ns"
              
              stats_text <- paste0(stats_text,
                                   sprintf("%s (#%s) vs %s (#%s): p = %.4f, %s\n",
                                           as.character(test_info[[i]]$test_mut),
                                           as.character(test_info[[i]]$test_group),
                                           as.character(ref_mut),
                                           as.character(ref_group),
                                           p_adjusted[i],
                                           sig))
            }
          }
          
        }
        
      } else {
        stats_text <- paste0(stats_text, "\n=== Statistical Analysis ===\n")
        stats_text <- paste0(stats_text, "Note: Insufficient matched data for repeated measures ANOVA.\n")
        stats_text <- paste0(stats_text, sprintf("Only %d experiments have data for all groups.\n", n_complete))
        stats_text <- paste0(stats_text, "Need at least 2 complete sets for RM-ANOVA.\n")
      }
    }
    
    return(stats_text)
  })
  
  # ==============================================================================
  # DOWNLOAD FOR PRISM
  # ==============================================================================
  
  output$download_prism <- downloadHandler(
    filename = function() {
      paste0("correlation_data_", format(Sys.Date(), "%Y%m%d"), ".xlsx")
    },
    content = function(file) {
      req(rv$comparison_samples())
      
      plot_data <- rv$comparison_samples()
      
      if(nrow(plot_data) == 0) {
        showNotification("No data selected for download", type = "error")
        return()
      }
      
      unique_groups <- plot_data %>%
        group_by(Cell_line, Mutation) %>%
        summarize(.groups = 'drop') %>%
        arrange(Cell_line)
      
      col_names <- paste0(unique_groups$Mutation, " (#", unique_groups$Cell_line, ")")
      
      unique_experiments <- unique(plot_data$Experiment)
      
      output_df <- data.frame(Experiment = unique_experiments, stringsAsFactors = FALSE)
      
      low_cell_warnings <- list()
      
      for(i in 1:nrow(unique_groups)) {
        cell_line <- unique_groups$Cell_line[i]
        mutation <- unique_groups$Mutation[i]
        col_name <- col_names[i]
        
        cell_data <- plot_data %>%
          filter(Cell_line == cell_line) %>%
          select(Experiment, Correlation, N_cells)
        
        values <- sapply(unique_experiments, function(exp) {
          match_idx <- which(cell_data$Experiment == exp)
          if(length(match_idx) > 0) {
            corr <- as.numeric(cell_data$Correlation[match_idx[1]])
            n_cells <- as.numeric(cell_data$N_cells[match_idx[1]])
            
            is_empty_vector <- grepl("Empty.?Vector", mutation, ignore.case = TRUE)
            if(!is.na(n_cells) && n_cells < 500 && !is_empty_vector) {
              if(is.null(low_cell_warnings[[exp]])) {
                low_cell_warnings[[exp]] <- character(0)
              }
              low_cell_warnings[[exp]] <- c(low_cell_warnings[[exp]], 
                                            paste0(mutation, " (#", cell_line, ")"))
            }
            
            if(!is.na(corr)) {
              return(round(corr, 4))
            } else {
              return(NA_real_)
            }
          } else {
            return(NA_real_)
          }
        })
        
        output_df[[col_name]] <- values
      }
      
      warnings_col <- sapply(unique_experiments, function(exp) {
        if(!is.null(low_cell_warnings[[exp]]) && length(low_cell_warnings[[exp]]) > 0) {
          return(paste(low_cell_warnings[[exp]], collapse = "; "))
        } else {
          return("")
        }
      })
      
      output_df$Low_Cell_Warning <- warnings_col
      
      if(!require("writexl", quietly = TRUE)) {
        install.packages("writexl", repos = "http://cran.r-project.org")
        library(writexl)
      }
      
      notes_df <- data.frame(
        Note = c("Correlation values are reported to 4 decimal places",
                 "Each row represents one experiment",
                 "Each column represents one cell line",
                 "Low_Cell_Warning column lists cell lines with N_cells < 500",
                 paste0("Data exported: ", Sys.time())),
        stringsAsFactors = FALSE
      )
      
      write_xlsx(list(
        "Correlation_Data" = output_df,
        "Notes" = notes_df
      ), path = file)
      
      showNotification("Data downloaded successfully!", type = "message", duration = 3)
    }
  )
}

# ==============================================================================
# RUN APP
# ==============================================================================

shinyApp(ui = ui, server = server)