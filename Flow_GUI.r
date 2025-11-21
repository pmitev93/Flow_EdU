# ==============================================================================
# FLOW CYTOMETRY ANALYSIS - SHINY APP
# ==============================================================================

library(shiny)
library(shinyjs)
library(flowCore)
library(tidyverse)
library(sp)
library(DT)
library(openxlsx)
library(sortable)
library(plotly)
library(future)
library(future.apply)

# Set up parallel processing plan
# Use multisession to load experiments in parallel
plan(multisession, workers = availableCores() - 1)

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
  useShinyjs(),
  titlePanel("The MITEV EdU Analysis Tool"),

  # Add custom CSS for sidebar collapse
  tags$head(
    tags$style(HTML("
      #sidebar_panel {
        transition: all 0.3s;
      }
      .sidebar-collapsed {
        display: none !important;
        width: 0 !important;
      }
      .main-expanded {
        width: 100% !important;
        max-width: 100% !important;
      }
      #toggle_sidebar {
        position: fixed;
        left: 0;
        top: 100px;
        z-index: 1000;
        border-radius: 0 5px 5px 0;
      }
      #toggle_sidebar.collapsed {
        left: 0;
      }

      /* Gating Strategy Creator panel toggles */
      .creator-panel-collapsed {
        display: none !important;
        width: 0 !important;
      }
      .creator-plot-width-12 {
        width: 100% !important;
        max-width: 100% !important;
      }
      .creator-plot-width-9 {
        width: 75% !important;
        max-width: 75% !important;
      }
      .creator-plot-width-8 {
        width: 66.66666667% !important;
        max-width: 66.66666667% !important;
      }
    "))
  ),

  # Toggle button (always visible)
  actionButton("toggle_sidebar",
               HTML("<i class='glyphicon glyphicon-chevron-left'></i>"),
               class = "btn-primary btn-sm",
               style = "position: fixed; left: 0; top: 100px; z-index: 1000; border-radius: 0 5px 5px 0;"),

  sidebarLayout(
    sidebarPanel(
      id = "sidebar_panel",
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
      id = "main_panel",
      width = 8,  # Adjusted from 9 to match sidebar width=4

      tabsetPanel(
        id = "main_tabs",
        
        # Welcome tab
        tabPanel("Welcome",
                 h3("Welcome to the Multivariate Identification and Tracking of EdU-incorporating Variants (MITEV) EdU Analysis Tool"),
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
                 fluidRow(
                   column(4, checkboxInput("show_dox_minus", "Show Dox- samples", value = FALSE)),
                   column(4, checkboxInput("show_unanalyzed", "Show unanalyzed samples", value = FALSE))
                 ),
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

        # Sample Overview Tab
        tabPanel("Sample Overview",
                 h3("All Gates for Single Sample"),
                 selectInput("sample_overview_experiment", "Select Experiment:",
                             choices = NULL),
                 selectInput("sample_overview_gate_strategy", "Gating Strategy:",
                             choices = NULL),
                 selectInput("sample_overview_sample", "Select Sample:",
                             choices = NULL),
                 plotOutput("sample_overview_plot", height = "1200px")
        ),

        # Gating Strategy Creator Tab
        tabPanel("Gating Strategy Creator",
                 h3("Create or Modify Gating Strategies"),
                 p("Load an existing strategy, modify parameters, and save as a new strategy file."),

                 # Toggle buttons for creator panels
                 div(style = "margin-bottom: 10px;",
                     actionButton("toggle_creator_load",
                                  HTML("<i class='glyphicon glyphicon-chevron-left'></i> Load/Sample"),
                                  class = "btn-primary btn-sm",
                                  style = "margin-right: 5px;"),
                     actionButton("toggle_creator_edit",
                                  HTML("<i class='glyphicon glyphicon-chevron-left'></i> Edit Gates"),
                                  class = "btn-primary btn-sm")
                 ),

                 fluidRow(
                   # Left panel: Load/Sample/Save controls
                   column(3,
                          id = "creator_load_panel",
                          wellPanel(
                            h4("1. Load Base Strategy"),
                            selectInput("creator_base_strategy", "Base Strategy:",
                                        choices = NULL),
                            actionButton("creator_load", "Load Strategy",
                                        class = "btn-primary btn-block"),
                            hr(),

                            h4("2. Select Sample for Preview"),
                            selectInput("creator_experiment", "Experiment:",
                                       choices = NULL),
                            selectInput("creator_sample", "Sample:",
                                       choices = NULL),
                            hr(),

                            h4("3. Save New Strategy"),
                            textInput("creator_new_id", "Strategy ID:",
                                     placeholder = "e.g., gid18"),
                            textInput("creator_new_name", "Strategy Name:",
                                     placeholder = "e.g., My Custom Strategy"),
                            textAreaInput("creator_new_desc", "Description:",
                                         placeholder = "Describe your changes...",
                                         rows = 3),
                            actionButton("creator_save", "Save as New Strategy",
                                        class = "btn-success btn-block")
                          )
                   ),

                   # Middle panel: All gate editing controls (scrollable)
                   column(4,
                          id = "creator_edit_panel",
                          wellPanel(
                            style = "overflow-y: auto; max-height: 900px;",
                            h4("Edit Gate Parameters"),
                            p(style = "font-size: 12px; color: #666;", "Click gate names to expand/collapse"),

                            # Gate 1: Debris (collapsible)
                            tags$a(href = "#collapse_debris", `data-toggle` = "collapse",
                                   h5(HTML("<i class='glyphicon glyphicon-chevron-down'></i> Gate 1: Debris (FSC-A, SSC-A)"),
                                      style = "color: #337ab7; margin-top: 15px; cursor: pointer;")),
                            tags$div(id = "collapse_debris", class = "collapse",
                                    uiOutput("creator_debris_ui")),
                            hr(),

                            # Gate 2: Singlets (collapsible)
                            tags$a(href = "#collapse_singlet", `data-toggle` = "collapse",
                                   h5(HTML("<i class='glyphicon glyphicon-chevron-down'></i> Gate 2: Singlets (FSC-A, FSC-H)"),
                                      style = "color: #337ab7; cursor: pointer;")),
                            tags$div(id = "collapse_singlet", class = "collapse",
                                    uiOutput("creator_singlet_ui")),
                            hr(),

                            # Gate 3: Live Cells (collapsible)
                            tags$a(href = "#collapse_live", `data-toggle` = "collapse",
                                   h5(HTML("<i class='glyphicon glyphicon-chevron-down'></i> Gate 3: Live Cells (DCM-A, SSC-A)"),
                                      style = "color: #337ab7; cursor: pointer;")),
                            tags$div(id = "collapse_live", class = "collapse",
                                    uiOutput("creator_live_ui")),
                            hr(),

                            # Gate 4: S-phase Outliers (collapsible)
                            tags$a(href = "#collapse_sphase", `data-toggle` = "collapse",
                                   h5(HTML("<i class='glyphicon glyphicon-chevron-down'></i> Gate 4: S-phase Outliers (FxCycle-A, EdU-A)"),
                                      style = "color: #337ab7; cursor: pointer;")),
                            tags$div(id = "collapse_sphase", class = "collapse",
                                    uiOutput("creator_sphase_ui")),
                            hr(),

                            # Gate 5: FxCycle Quantile (collapsible, default open)
                            tags$a(href = "#collapse_fxcycle", `data-toggle` = "collapse",
                                   h5(HTML("<i class='glyphicon glyphicon-chevron-down'></i> Gate 5: FxCycle Quantile"),
                                      style = "color: #337ab7; cursor: pointer;")),
                            tags$div(id = "collapse_fxcycle", class = "collapse in",
                                    uiOutput("creator_fxcycle_ui")),
                            hr(),

                            # Gate 6: EdU + FxCycle (collapsible, default open)
                            tags$a(href = "#collapse_edu_fxcycle", `data-toggle` = "collapse",
                                   h5(HTML("<i class='glyphicon glyphicon-chevron-down'></i> Gate 6: EdU + FxCycle"),
                                      style = "color: #337ab7; cursor: pointer;")),
                            tags$div(id = "collapse_edu_fxcycle", class = "collapse in",
                                    uiOutput("creator_edu_fxcycle_ui")),
                            hr(),

                            # Gate 7: HA Positive / Quadrant (collapsible, default open)
                            uiOutput("creator_gate7_header"),
                            tags$div(id = "collapse_ha", class = "collapse in",
                                    uiOutput("creator_ha_ui"))
                          )
                   ),

                   # Right panel: Plot viewer with toggle
                   column(5,
                          id = "creator_plot_panel",
                          wellPanel(
                            radioButtons("creator_plot_mode", "View Mode:",
                                        choices = c("All Gates" = "all",
                                                   "Individual Plot" = "single"),
                                        selected = "all",
                                        inline = TRUE),
                            conditionalPanel(
                              condition = "input.creator_plot_mode == 'single'",
                              selectInput("creator_single_plot", "Select Plot:",
                                         choices = c("Gate 1: Debris" = "gate1",
                                                    "Gate 2: Singlets" = "gate2",
                                                    "Gate 3: Live Cells" = "gate3",
                                                    "Gate 4: S-phase Outliers" = "gate4",
                                                    "Gate 5: FxCycle Quantile" = "gate5",
                                                    "Gate 6: EdU + FxCycle" = "gate6",
                                                    "Gate 7" = "gate7",
                                                    "Final Correlation" = "correlation"),
                                         selected = "correlation")
                            ),
                            conditionalPanel(
                              condition = "input.creator_plot_mode == 'all'",
                              h4(textOutput("creator_preview_sample_title"), align = "center", style = "font-weight: bold; margin-bottom: 10px;"),
                              plotOutput("creator_preview_plot", height = "900px")
                            ),
                            conditionalPanel(
                              condition = "input.creator_plot_mode == 'single'",
                              plotOutput("creator_single_plot_output", height = "500px")
                            )
                          )
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
                   column(4, checkboxInput("msc_show_dox_minus", "Show Dox- samples", value = FALSE)),
                   column(4, checkboxInput("msc_show_unanalyzed", "Show unanalyzed samples", value = FALSE))
                 ),
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
  # SIDEBAR TOGGLE
  # ==============================================================================

  # Track sidebar state
  sidebar_visible <- reactiveVal(TRUE)

  observeEvent(input$toggle_sidebar, {
    sidebar_visible(!sidebar_visible())
    if(sidebar_visible()) {
      # Show sidebar
      shinyjs::removeClass(id = "sidebar_panel", class = "sidebar-collapsed")
      shinyjs::removeClass(id = "main_panel", class = "main-expanded")
      shinyjs::html(id = "toggle_sidebar",
                   html = "<i class='glyphicon glyphicon-chevron-left'></i>")
    } else {
      # Hide sidebar
      shinyjs::addClass(id = "sidebar_panel", class = "sidebar-collapsed")
      shinyjs::addClass(id = "main_panel", class = "main-expanded")
      shinyjs::html(id = "toggle_sidebar",
                   html = "<i class='glyphicon glyphicon-chevron-right'></i>")
    }
  })

  # ==============================================================================
  # GATING STRATEGY CREATOR PANEL TOGGLES
  # ==============================================================================

  # Track creator panel states
  creator_load_visible <- reactiveVal(TRUE)
  creator_edit_visible <- reactiveVal(TRUE)

  # Function to update plot panel width based on visible panels
  update_creator_plot_width <- function() {
    load_vis <- creator_load_visible()
    edit_vis <- creator_edit_visible()

    # Remove all width classes first
    shinyjs::removeClass(id = "creator_plot_panel", class = "creator-plot-width-12")
    shinyjs::removeClass(id = "creator_plot_panel", class = "creator-plot-width-9")
    shinyjs::removeClass(id = "creator_plot_panel", class = "creator-plot-width-8")

    # Add appropriate width class based on visible panels
    if(!load_vis && !edit_vis) {
      # Both hidden - full width (12 cols)
      shinyjs::addClass(id = "creator_plot_panel", class = "creator-plot-width-12")
    } else if(load_vis && !edit_vis) {
      # Only load visible - 9 cols for plot
      shinyjs::addClass(id = "creator_plot_panel", class = "creator-plot-width-9")
    } else if(!load_vis && edit_vis) {
      # Only edit visible - 8 cols for plot
      shinyjs::addClass(id = "creator_plot_panel", class = "creator-plot-width-8")
    }
    # If both visible, default width (5 cols) is used
  }

  # Toggle Load/Sample panel
  observeEvent(input$toggle_creator_load, {
    creator_load_visible(!creator_load_visible())
    if(creator_load_visible()) {
      # Show panel
      shinyjs::removeClass(id = "creator_load_panel", class = "creator-panel-collapsed")
      shinyjs::html(id = "toggle_creator_load",
                   html = "<i class='glyphicon glyphicon-chevron-left'></i> Load/Sample")
    } else {
      # Hide panel
      shinyjs::addClass(id = "creator_load_panel", class = "creator-panel-collapsed")
      shinyjs::html(id = "toggle_creator_load",
                   html = "<i class='glyphicon glyphicon-chevron-right'></i> Load/Sample")
    }
    update_creator_plot_width()
  })

  # Toggle Edit Gates panel
  observeEvent(input$toggle_creator_edit, {
    creator_edit_visible(!creator_edit_visible())
    if(creator_edit_visible()) {
      # Show panel
      shinyjs::removeClass(id = "creator_edit_panel", class = "creator-panel-collapsed")
      shinyjs::html(id = "toggle_creator_edit",
                   html = "<i class='glyphicon glyphicon-chevron-left'></i> Edit Gates")
    } else {
      # Hide panel
      shinyjs::addClass(id = "creator_edit_panel", class = "creator-panel-collapsed")
      shinyjs::html(id = "toggle_creator_edit",
                   html = "<i class='glyphicon glyphicon-chevron-right'></i> Edit Gates")
    }
    update_creator_plot_width()
  })

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

  # User preferences file
  PREFS_FILE <- file.path(CACHE_DIR, "user_preferences.rds")

  # Save user preference for gating strategy per experiment
  save_gating_preference <- function(experiment_name, gating_strategy) {
    prefs <- list()
    if(file.exists(PREFS_FILE)) {
      prefs <- readRDS(PREFS_FILE)
    }

    if(is.null(prefs$gating_strategies)) {
      prefs$gating_strategies <- list()
    }

    prefs$gating_strategies[[experiment_name]] <- gating_strategy
    saveRDS(prefs, PREFS_FILE)
  }

  # Load saved gating strategy preference for an experiment
  load_gating_preference <- function(experiment_name) {
    if(!file.exists(PREFS_FILE)) {
      return(NULL)
    }

    prefs <- readRDS(PREFS_FILE)
    if(is.null(prefs$gating_strategies)) {
      return(NULL)
    }

    return(prefs$gating_strategies[[experiment_name]])
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
  save_to_cache <- function(experiment_name, results, ha_threshold, gates, ha_percentile = 0.98, gate_strategy_id = NULL, gate_strategy = NULL) {
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
      gate_strategy = gate_strategy,  # Store GATE_STRATEGY metadata
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
    gate_strategies = list(),  # Store GATE_STRATEGY metadata for each gate_id
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
          cat("\n=== EXPERIMENT SELECTOR DEBUG ===\n")
          exp_info <- lapply(exp_names_local, function(exp_name) {
            # Scan all subfolders (gating strategies) for cache files
            pattern <- paste0("^", exp_name, "_.*\\.rds$")
            cache_files <- list.files(CACHE_DIR, pattern = pattern, full.names = TRUE, recursive = TRUE)

            cat(sprintf("Experiment: %s, found %d cache files\n", exp_name, length(cache_files)))

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
                  cat(sprintf("  File: %s -> Gate ID: %s\n", basename(cache_file), gate_id))
                  gate_ids <- c(gate_ids, gate_id)
                }, error = function(e) {
                  # Skip files that can't be read
                  cat(sprintf("  ERROR reading %s: %s\n", basename(cache_file), e$message))
                })
              }

              # Remove duplicates and sort
              gate_ids <- unique(gate_ids)
              gate_ids <- sort(gate_ids)

              cat(sprintf("  Final gate IDs for %s: %s\n", exp_name, paste(gate_ids, collapse = ", ")))

              return(list(analyzed = TRUE, gate_ids = gate_ids))
            } else {
              return(list(analyzed = FALSE, gate_ids = character(0)))
            }
          })
          names(exp_info) <- exp_names_local
          cat("=== END EXPERIMENT SELECTOR DEBUG ===\n\n")

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

        cat("\n=== AUTO-LOAD DEBUG ===\n")
        cat(sprintf("Initial rv$all_results has %d rows\n", nrow(rv$all_results)))

        # Convert Correlation and N_cells to proper types to avoid bind_rows errors
        # (quick_scan creates them as character, but cache has numeric)
        if("Correlation" %in% names(rv$all_results)) {
          rv$all_results$Correlation <- ifelse(
            rv$all_results$Correlation == "Not analyzed",
            NA_real_,
            as.numeric(rv$all_results$Correlation)
          )
          cat("Converted Correlation column to numeric\n")
        }
        if("N_cells" %in% names(rv$all_results)) {
          rv$all_results$N_cells <- ifelse(
            rv$all_results$N_cells == "Not analyzed",
            NA_integer_,
            as.integer(rv$all_results$N_cells)
          )
          cat("Converted N_cells column to numeric\n")
        }

        # First pass: Load cached results and identify experiments that need FCS loading
        experiments_to_load <- list()

        for(i in seq_along(exp_names)) {
          exp_name <- exp_names[i]
          incProgress(1/length(exp_names), detail = exp_name)

          # Find ALL cache files for this experiment (scan all subfolders for all gate strategies)
          pattern <- paste0("^", exp_name, "_.*\\.rds$")
          cache_files <- list.files(CACHE_DIR, pattern = pattern, full.names = TRUE, recursive = TRUE)

          cat(sprintf("\nExperiment: %s\n", exp_name))
          cat(sprintf("  Found %d cache files: %s\n", length(cache_files),
                      paste(basename(cache_files), collapse = ", ")))

          if(length(cache_files) > 0) {
            # Load ALL cache files (one for each gate strategy)
            for(cache_file in cache_files) {
              tryCatch({
                cache_data <- readRDS(cache_file)

                # Add Gate_ID column if it doesn't exist
                if(!"Gate_ID" %in% names(rv$all_results)) {
                  rv$all_results$Gate_ID <- NA_character_
                  cat("  Added Gate_ID column to rv$all_results\n")
                }

                gate_id <- if(!is.null(cache_data$gate_id)) {
                  cache_data$gate_id
                } else {
                  get_readable_gate_id(cache_data$fingerprint)
                }

                cat(sprintf("  Loading cache file: %s (Gate ID: %s, %d wells)\n",
                            basename(cache_file), gate_id, nrow(cache_data$results)))

                # Update results table with cached data
                n_updated <- 0
                n_added <- 0
                for(j in seq_len(nrow(cache_data$results))) {
                  # First try to match on Experiment, Well, AND Gate_ID
                  match_idx <- which(rv$all_results$Experiment == cache_data$results$Experiment[j] &
                                       rv$all_results$Well == cache_data$results$Well[j] &
                                       !is.na(rv$all_results$Gate_ID) &
                                       rv$all_results$Gate_ID == gate_id)

                  # If no match found, try to find an unanalyzed row (Gate_ID is NA)
                  if(length(match_idx) == 0) {
                    match_idx <- which(rv$all_results$Experiment == cache_data$results$Experiment[j] &
                                         rv$all_results$Well == cache_data$results$Well[j] &
                                         is.na(rv$all_results$Gate_ID))
                  }

                  if(length(match_idx) > 0) {
                    # Update existing row (take first match if multiple)
                    match_idx <- match_idx[1]
                    rv$all_results$Correlation[match_idx] <- cache_data$results$Correlation[j]
                    rv$all_results$N_cells[match_idx] <- cache_data$results$N_cells[j]
                    rv$all_results$Notes[match_idx] <- cache_data$results$Notes[j]
                    rv$all_results$Gate_ID[match_idx] <- gate_id

                    # Update Strength_Ratio if it exists in cached results
                    if("Strength_Ratio" %in% names(cache_data$results)) {
                      if(!"Strength_Ratio" %in% names(rv$all_results)) {
                        rv$all_results$Strength_Ratio <- NA_real_
                      }
                      rv$all_results$Strength_Ratio[match_idx] <- cache_data$results$Strength_Ratio[j]
                    }

                    # Update HA_Pos_Pct if it exists in cached results
                    if("HA_Pos_Pct" %in% names(cache_data$results)) {
                      if(!"HA_Pos_Pct" %in% names(rv$all_results)) {
                        rv$all_results$HA_Pos_Pct <- NA_real_
                      }
                      rv$all_results$HA_Pos_Pct[match_idx] <- cache_data$results$HA_Pos_Pct[j]
                    }

                    n_updated <- n_updated + 1
                  } else {
                    # This is a second/third gate strategy for same well - add new row
                    new_row <- cache_data$results[j, ]
                    new_row$Gate_ID <- gate_id

                    # Ensure Strength_Ratio column exists before binding
                    if("Strength_Ratio" %in% names(cache_data$results) && !"Strength_Ratio" %in% names(rv$all_results)) {
                      rv$all_results$Strength_Ratio <- NA_real_
                    }

                    # Ensure HA_Pos_Pct column exists before binding
                    if("HA_Pos_Pct" %in% names(cache_data$results) && !"HA_Pos_Pct" %in% names(rv$all_results)) {
                      rv$all_results$HA_Pos_Pct <- NA_real_
                    }

                    rv$all_results <- bind_rows(rv$all_results, new_row)
                    n_added <- n_added + 1
                  }
                }

                cat(sprintf("    Updated %d rows, added %d rows\n", n_updated, n_added))

                # Store HA threshold
                rv$ha_thresholds[[exp_name]] <- cache_data$ha_threshold

                # Store gates used for this experiment+strategy combination
                if(!is.null(cache_data$gates)) {
                  composite_key <- paste0(exp_name, "::", gate_id)
                  rv$experiment_gates[[composite_key]] <- cache_data$gates
                  cat(sprintf("    Stored gates with key: %s\n", composite_key))
                }

                # Store gate strategy metadata
                if(!is.null(cache_data$gate_strategy)) {
                  gate_strategy_key <- paste0("GATE_STRATEGY_", gate_id)
                  rv$gate_strategies[[gate_strategy_key]] <- cache_data$gate_strategy
                  cat(sprintf("    Stored gate_strategy with key: %s\n", gate_strategy_key))
                }

                # Mark this experiment for FCS loading (only once per experiment)
                if(!exp_name %in% names(experiments_to_load)) {
                  experiments_to_load[[exp_name]] <- exp_folders[i]
                }

                n_loaded <- n_loaded + 1

              }, error = function(e) {
                cat(sprintf("  ERROR loading cache for %s: %s\n", basename(cache_file), e$message))
              })
            }
          }
        }

        # Second pass: Load FCS data in parallel for all experiments
        if(length(experiments_to_load) > 0) {
          cat(sprintf("\n=== PARALLEL LOADING %d EXPERIMENTS ===\n", length(experiments_to_load)))

          # Initialize experiments list if needed
          if(is.null(rv$experiments)) {
            rv$experiments <- list()
          }

          # Load experiments in parallel
          loaded_experiments <- future_lapply(experiments_to_load, function(exp_folder) {
            load_experiment(exp_folder)
          }, future.seed = TRUE)

          # Store loaded experiments in reactive values
          names(loaded_experiments) <- names(experiments_to_load)
          rv$experiments <- c(rv$experiments, loaded_experiments)

          cat(sprintf("=== PARALLEL LOADING COMPLETE ===\n"))
        }

        cat(sprintf("\nFinal rv$all_results has %d rows\n", nrow(rv$all_results)))
        cat(sprintf("Gate_ID distribution:\n"))
        print(table(rv$all_results$Gate_ID, useNA = "ifany"))
        cat("=== END AUTO-LOAD DEBUG ===\n\n")
        
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

    # Try to load saved preference for this experiment
    saved_strategy <- load_gating_preference(exp_name)
    default_strategy <- available_gates[1]

    # Use saved strategy if it's still available, otherwise use first available
    if(!is.null(saved_strategy) && saved_strategy %in% available_gates) {
      default_strategy <- saved_strategy
    }

    updateSelectInput(session, "browse_gate_strategy",
                      choices = available_gates,
                      selected = default_strategy)
  })

  # Save gating strategy preference when user changes it
  observeEvent(input$browse_gate_strategy, {
    req(input$selected_experiment, input$browse_gate_strategy)

    exp_name <- input$selected_experiment
    gate_strategy <- input$browse_gate_strategy

    # Save the preference
    save_gating_preference(exp_name, gate_strategy)
  }, ignoreInit = TRUE)  # Don't save on initial load, only on user changes

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
        # Use filename-based ID to avoid collisions (e.g., quadrant vs quadrant2)
        gate_id_from_file <- gsub("^gates_(.*)\\.r$", "\\1", gate_file)

        gate_id_for_cache <- if(!is.null(GATE_STRATEGY_selected$id)) {
          # Prefer filename-based ID to avoid collisions
          gate_id_from_file
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

          # Use gate strategy ID (prefer from cache, fallback to filename-based ID)
          gate_id <- if(!is.null(cache_data$gate_id)) {
            cache_data$gate_id
          } else {
            # For old cache files, use filename-based ID to avoid collisions
            gate_id_from_file
          }

          # Store gates used for this experiment+strategy combination (prefer from cache, fallback to selected)
          composite_key <- paste0(exp_name, "::", gate_id)
          rv$experiment_gates[[composite_key]] <- if(!is.null(cache_data$gates)) {
            cache_data$gates
          } else {
            GATES_selected
          }

          # Store gate strategy metadata (prefer from cache, fallback to loaded from file)
          gate_strategy_key <- paste0("GATE_STRATEGY_", gate_id)
          if(!is.null(cache_data$gate_strategy)) {
            rv$gate_strategies[[gate_strategy_key]] <- cache_data$gate_strategy
            cat(sprintf("  Stored gate_strategy from cache with key: %s\n", gate_strategy_key))
          } else if(!is.null(GATE_STRATEGY_selected)) {
            rv$gate_strategies[[gate_strategy_key]] <- GATE_STRATEGY_selected
            cat(sprintf("  Stored gate_strategy from file with key: %s\n", gate_strategy_key))
          } else {
            cat(sprintf("  WARNING: No gate_strategy to store for key: %s\n", gate_strategy_key))
          }

          exp_results$Gate_ID <- gate_id

          n_from_cache <- n_from_cache + 1
        } else {
          # Run analysis
          incProgress(1/(n_exp*2), detail = sprintf("Analyzing %s", exp_name))

          rv$ha_thresholds[[exp_name]] <- ha_threshold

          # Check if using quadrant strategy
          use_quadrant <- !is.null(GATE_STRATEGY_selected$analysis_type) &&
                          GATE_STRATEGY_selected$analysis_type == "quadrant_ratio"

          # Use appropriate analysis function
          # For quadrant analysis, pass NULL for ha_threshold to skip Gate 7 filtering
          exp_results <- extract_correlations_with_quadrants(exp,
                                                             if(use_quadrant) NULL else ha_threshold,
                                                             gates = GATES_selected,
                                                             channels = CHANNELS,
                                                             use_quadrant = use_quadrant)

          # Save to cache with gate strategy ID
          # Use same filename-based ID as above to avoid collisions
          gate_id <- gate_id_from_file

          save_to_cache(exp_name, exp_results, ha_threshold, GATES_selected,
                        gate_strategy_id = gate_id,
                        gate_strategy = GATE_STRATEGY_selected)

          # Store gates used for this experiment+strategy combination
          composite_key <- paste0(exp_name, "::", gate_id)
          rv$experiment_gates[[composite_key]] <- GATES_selected

          # Store gate strategy metadata
          gate_strategy_key <- paste0("GATE_STRATEGY_", gate_id)
          if(!is.null(GATE_STRATEGY_selected)) {
            rv$gate_strategies[[gate_strategy_key]] <- GATE_STRATEGY_selected
            cat(sprintf("  Stored gate_strategy from newly analyzed data with key: %s (analysis_type: %s)\n",
                        gate_strategy_key,
                        if(!is.null(GATE_STRATEGY_selected$analysis_type)) GATE_STRATEGY_selected$analysis_type else "NULL"))
          } else {
            cat(sprintf("  WARNING: GATE_STRATEGY_selected is NULL for key: %s\n", gate_strategy_key))
          }

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

        # First try to match on Experiment, Well, AND Gate_ID
        match_idx <- which(rv$all_results$Experiment == new_results$Experiment[i] &
                             rv$all_results$Well == new_results$Well[i] &
                             !is.na(rv$all_results$Gate_ID) &
                             rv$all_results$Gate_ID == new_results$Gate_ID[i])

        # If no match found, try to find an unanalyzed row (Gate_ID is NA)
        if(length(match_idx) == 0) {
          match_idx <- which(rv$all_results$Experiment == new_results$Experiment[i] &
                               rv$all_results$Well == new_results$Well[i] &
                               is.na(rv$all_results$Gate_ID))
        }

        if(length(match_idx) > 0) {
          # Update existing row (take first match if multiple)
          match_idx <- match_idx[1]
          rv$all_results$Correlation[match_idx] <- new_results$Correlation[i]
          rv$all_results$N_cells[match_idx] <- new_results$N_cells[i]
          rv$all_results$Notes[match_idx] <- new_results$Notes[i]
          rv$all_results$Gate_ID[match_idx] <- new_results$Gate_ID[i]

          # Update Strength_Ratio column if it exists in new_results
          if("Strength_Ratio" %in% names(new_results)) {
            if(!"Strength_Ratio" %in% names(rv$all_results)) {
              rv$all_results$Strength_Ratio <- NA_real_
            }
            rv$all_results$Strength_Ratio[match_idx] <- new_results$Strength_Ratio[i]
          }

          # Update HA_Pos_Pct column if it exists in new_results
          if("HA_Pos_Pct" %in% names(new_results)) {
            if(!"HA_Pos_Pct" %in% names(rv$all_results)) {
              rv$all_results$HA_Pos_Pct <- NA_real_
            }
            rv$all_results$HA_Pos_Pct[match_idx] <- new_results$HA_Pos_Pct[i]
          }
        } else {
          # This is a second/third gate strategy for same well - add new row
          # Ensure Strength_Ratio column exists in rv$all_results before binding
          if("Strength_Ratio" %in% names(new_results) && !"Strength_Ratio" %in% names(rv$all_results)) {
            rv$all_results$Strength_Ratio <- NA_real_
          }

          # Ensure HA_Pos_Pct column exists in rv$all_results before binding
          if("HA_Pos_Pct" %in% names(new_results) && !"HA_Pos_Pct" %in% names(rv$all_results)) {
            rv$all_results$HA_Pos_Pct <- NA_real_
          }
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

    # React to ui_refresh_trigger to update after analysis
    rv$ui_refresh_trigger

    # Start with all results
    display_data <- rv$all_results

    # Filter out Dox- samples unless checkbox is checked
    if(!isTRUE(input$show_dox_minus)) {
      display_data <- display_data[!grepl("Dox-", display_data$Sample, ignore.case = TRUE), ]
    }

    # Filter out unanalyzed samples unless checkbox is checked
    if(!isTRUE(input$show_unanalyzed)) {
      display_data <- display_data[!is.na(display_data$Correlation) & display_data$Correlation != "Not analyzed", ]
    }

    if(nrow(display_data) == 0) {
      return(datatable(data.frame(Message = "No samples match current filters")))
    }

    # Format correlation to 4 decimal places (handle both character and numeric)
    if("Correlation" %in% names(display_data)) {
      display_data$Correlation <- ifelse(
        is.na(display_data$Correlation) | display_data$Correlation == "Not analyzed",
        "Not analyzed",
        sprintf("%.4f", as.numeric(display_data$Correlation))
      )
    }

    # Format ratio to 4 decimal places
    if("Strength_Ratio" %in% names(display_data)) {
      display_data$Strength_Ratio <- ifelse(
        is.na(display_data$Strength_Ratio),
        "",
        sprintf("%.4f", as.numeric(display_data$Strength_Ratio))
      )
    }

    # Format HA_Pos_Pct to 2 decimal places
    if("HA_Pos_Pct" %in% names(display_data)) {
      display_data$HA_Pos_Pct <- ifelse(
        is.na(display_data$HA_Pos_Pct),
        "",
        sprintf("%.2f%%", as.numeric(display_data$HA_Pos_Pct))
      )
    }

    # Ensure columns are character/numeric for proper filtering
    if("N_cells" %in% names(display_data)) {
      display_data$N_cells <- ifelse(
        is.na(display_data$N_cells) | display_data$N_cells == "Not analyzed",
        "Not analyzed",
        as.character(display_data$N_cells)
      )
    }
    if("Notes" %in% names(display_data)) {
      display_data$Notes <- as.character(display_data$Notes)
    }
    if("Gate_ID" %in% names(display_data)) {
      display_data$Gate_ID <- as.character(display_data$Gate_ID)
    }

    # Calculate summary statistics per Cell_line
    if("Correlation" %in% names(rv$all_results) && "Cell_line" %in% names(rv$all_results)) {
      # Get numeric correlations for summary
      numeric_corr <- as.numeric(rv$all_results$Correlation)
      valid_data <- rv$all_results[!is.na(numeric_corr), ]
      valid_data$Correlation_num <- as.numeric(valid_data$Correlation)

      if(nrow(valid_data) > 0) {
        # Calculate stats per cell line (all experiments)
        all_stats <- valid_data %>%
          group_by(Cell_line) %>%
          summarize(
            Avg_All = mean(Correlation_num, na.rm = TRUE),
            SD_All = ifelse(n() > 1, sd(Correlation_num, na.rm = TRUE), NA),
            n_all = n(),
            .groups = 'drop'
          )

        # Calculate stats per cell line and week (same week experiments)
        week_stats <- NULL
        if("Year_Week" %in% names(valid_data)) {
          week_stats <- valid_data %>%
            group_by(Cell_line, Year_Week) %>%
            summarize(
              Avg_Week = mean(Correlation_num, na.rm = TRUE),
              SD_Week = ifelse(n() > 1, sd(Correlation_num, na.rm = TRUE), NA),
              n_week = n(),
              .groups = 'drop'
            )
        }

        # Add 4 columns to display_data
        display_data$Avg_All <- ""
        display_data$SD_All <- ""
        display_data$Avg_Week <- ""
        display_data$SD_Week <- ""

        # Populate per-row stats
        for(i in seq_len(nrow(display_data))) {
          cell_line <- display_data$Cell_line[i]
          # All experiments stats
          cell_stats <- all_stats[all_stats$Cell_line == cell_line, ]
          if(nrow(cell_stats) > 0) {
            display_data$Avg_All[i] <- sprintf("%.4f", cell_stats$Avg_All[1])
            if(!is.na(cell_stats$SD_All[1])) {
              display_data$SD_All[i] <- sprintf("%.4f", cell_stats$SD_All[1])
            }
          }
          # Same week stats
          if(!is.null(week_stats) && "Year_Week" %in% names(display_data)) {
            year_week <- display_data$Year_Week[i]
            wk_stats <- week_stats[week_stats$Cell_line == cell_line & week_stats$Year_Week == year_week, ]
            if(nrow(wk_stats) > 0) {
              display_data$Avg_Week[i] <- sprintf("%.4f", wk_stats$Avg_Week[1])
              if(!is.na(wk_stats$SD_Week[1])) {
                display_data$SD_Week[i] <- sprintf("%.4f", wk_stats$SD_Week[1])
              }
            }
          }
        }
      }
    }

    datatable(display_data,
              selection = 'single',  # Enable single row selection
              options = list(
                pageLength = 25,
                scrollX = TRUE,
                search = list(regex = FALSE, caseInsensitive = TRUE)
              ),
              filter = list(position = 'top', clear = FALSE))
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

    # Try to load saved preference for this experiment
    saved_strategy <- load_gating_preference(exp_name)
    default_strategy <- available_gates[1]

    # Use saved strategy if it's still available, otherwise use first available
    if(!is.null(saved_strategy) && saved_strategy %in% available_gates) {
      default_strategy <- saved_strategy
    }

    updateSelectInput(session, "overview_gate_strategy",
                      choices = available_gates,
                      selected = default_strategy)
  })

  # Save overview gating strategy preference when user changes it
  observeEvent(input$overview_gate_strategy, {
    req(input$overview_experiment, input$overview_gate_strategy)

    exp_name <- input$overview_experiment
    gate_strategy <- input$overview_gate_strategy

    # Save the preference
    save_gating_preference(exp_name, gate_strategy)
  }, ignoreInit = TRUE)  # Don't save on initial load, only on user changes

  # Update sample overview experiment selector when experiments are loaded
  observe({
    req(rv$experiments)
    updateSelectInput(session, "sample_overview_experiment",
                      choices = names(rv$experiments))
  })

  # Update sample overview gating strategy selector when experiment changes
  observeEvent(input$sample_overview_experiment, {
    req(input$sample_overview_experiment)
    exp_name <- input$sample_overview_experiment
    available_gates <- rv$experiment_available_gates[[exp_name]]
    if(is.null(available_gates) || length(available_gates) == 0) {
      available_gates <- "gdef"  # Default if none available
    }

    # Try to load saved preference for this experiment
    saved_strategy <- load_gating_preference(exp_name)
    default_strategy <- available_gates[1]

    # Use saved strategy if it's still available, otherwise use first available
    if(!is.null(saved_strategy) && saved_strategy %in% available_gates) {
      default_strategy <- saved_strategy
    }

    updateSelectInput(session, "sample_overview_gate_strategy",
                      choices = available_gates,
                      selected = default_strategy)
  })

  # Update sample overview sample selector when experiment or strategy changes
  observeEvent(c(input$sample_overview_experiment, input$sample_overview_gate_strategy), {
    req(input$sample_overview_experiment, rv$experiments)
    exp <- rv$experiments[[input$sample_overview_experiment]]
    if(!is.null(exp)) {
      sample_choices <- setNames(seq_along(exp$metadata$sample_name),
                                  exp$metadata$sample_name)
      updateSelectInput(session, "sample_overview_sample",
                        choices = sample_choices)
    }
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
      # Get strategy metadata
      gate_strategy_key <- paste0("GATE_STRATEGY_", input$overview_gate_strategy)

      # Debug: Show what's available in rv$gate_strategies
      cat(sprintf("\n=== Checking rv$gate_strategies ===\n"))
      cat(sprintf("Looking for key: %s\n", gate_strategy_key))
      cat(sprintf("Available keys in rv$gate_strategies: %s\n",
                  if(length(names(rv$gate_strategies)) > 0) paste(names(rv$gate_strategies), collapse=", ") else "NONE"))

      gate_strategy <- if(!is.null(rv$gate_strategies[[gate_strategy_key]])) {
        rv$gate_strategies[[gate_strategy_key]]
      } else {
        NULL
      }

      # If gate strategy not in memory, try to load from cache
      if(is.null(gate_strategy)) {
        cat(sprintf("Gate strategy NULL - attempting to load from cache for '%s'\n", input$overview_gate_strategy))
        cache_dir <- file.path("analysis_cache", exp_name)
        cat(sprintf("  Cache dir: %s (exists: %s)\n", cache_dir, dir.exists(cache_dir)))

        if(dir.exists(cache_dir)) {
          cache_files <- list.files(cache_dir, pattern = "\\.rds$", full.names = TRUE)
          cat(sprintf("  Found %d cache files\n", length(cache_files)))

          for(cache_file in cache_files) {
            tryCatch({
              cache_data <- readRDS(cache_file)
              cat(sprintf("    Checking %s: gate_id=%s (looking for '%s')\n",
                          basename(cache_file),
                          if(!is.null(cache_data$gate_id)) cache_data$gate_id else "NULL",
                          input$overview_gate_strategy))

              if(!is.null(cache_data$gate_id) && cache_data$gate_id == input$overview_gate_strategy) {
                # Found matching cache file, extract gate strategy AND gates
                cat(sprintf("    MATCH FOUND! gate_strategy present: %s, gates present: %s\n",
                            !is.null(cache_data$gate_strategy), !is.null(cache_data$gates)))

                if(!is.null(cache_data$gate_strategy)) {
                  gate_strategy <- cache_data$gate_strategy
                  rv$gate_strategies[[gate_strategy_key]] <- gate_strategy

                  # Also load the gates themselves
                  if(!is.null(cache_data$gates)) {
                    composite_key <- paste0(exp_name, "::", input$overview_gate_strategy)
                    rv$experiment_gates[[composite_key]] <- cache_data$gates
                    gates_to_use <- cache_data$gates
                    cat(sprintf("    ✓ Loaded gates and strategy '%s' from cache\n", input$overview_gate_strategy))
                    cat(sprintf("    gates_to_use$quadrant present: %s\n", !is.null(gates_to_use$quadrant)))
                  } else {
                    cat(sprintf("    WARNING: gates data missing from cache\n"))
                  }
                  break
                }
              }
            }, error = function(e) {
              cat(sprintf("    ERROR reading %s: %s\n", basename(cache_file), e$message))
            })
          }
        }
      } else {
        cat(sprintf("Gate strategy already in memory for '%s'\n", input$overview_gate_strategy))
      }

      # Re-evaluate gates_to_use after cache loading to ensure we have the latest gates
      composite_key <- paste0(exp_name, "::", input$overview_gate_strategy)
      if(!is.null(rv$experiment_gates[[composite_key]])) {
        gates_to_use <- rv$experiment_gates[[composite_key]]
        cat(sprintf("Updated gates_to_use from rv$experiment_gates for '%s'\n", input$overview_gate_strategy))
      }

      # Debug output
      cat(sprintf("\n=== Gate 7 Overview Debug ===\n"))
      cat(sprintf("overview_gate_strategy: %s\n", input$overview_gate_strategy))
      cat(sprintf("gate_strategy_key: %s\n", gate_strategy_key))
      cat(sprintf("gate_strategy: %s\n", if(!is.null(gate_strategy)) "PRESENT" else "NULL"))
      if(!is.null(gate_strategy)) {
        cat(sprintf("  analysis_type: %s\n", gate_strategy$analysis_type))
      }
      cat(sprintf("gates$quadrant: %s\n", if(!is.null(gates_to_use$quadrant)) "PRESENT" else "NULL"))

      # Check if using quadrant strategy
      use_quadrant <- !is.null(gate_strategy$analysis_type) &&
                      gate_strategy$analysis_type == "quadrant_ratio" &&
                      !is.null(gates_to_use$quadrant)

      cat(sprintf("use_quadrant: %s\n", use_quadrant))

      if(use_quadrant) {
        # Quadrant strategy: show quadrant plots for all Dox+ samples
        plot_quadrant_correlation_overview(exp, gates = gates_to_use, channels = CHANNELS)
      } else {
        # Old strategy: Calculate HA threshold
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
      }

    } else if(input$overview_gate == "correlation") {
      # Get strategy metadata
      gate_strategy_key <- paste0("GATE_STRATEGY_", input$overview_gate_strategy)
      gate_strategy <- if(!is.null(rv$gate_strategies[[gate_strategy_key]])) {
        rv$gate_strategies[[gate_strategy_key]]
      } else {
        NULL
      }

      # If gate strategy not in memory, try to load from cache
      if(is.null(gate_strategy)) {
        cache_dir <- file.path("analysis_cache", exp_name)
        if(dir.exists(cache_dir)) {
          cache_files <- list.files(cache_dir, pattern = "\\.rds$", full.names = TRUE)
          for(cache_file in cache_files) {
            tryCatch({
              cache_data <- readRDS(cache_file)
              if(!is.null(cache_data$gate_id) && cache_data$gate_id == input$overview_gate_strategy) {
                # Found matching cache file, extract gate strategy AND gates
                if(!is.null(cache_data$gate_strategy)) {
                  gate_strategy <- cache_data$gate_strategy
                  rv$gate_strategies[[gate_strategy_key]] <- gate_strategy

                  # Also load the gates themselves
                  if(!is.null(cache_data$gates)) {
                    composite_key <- paste0(exp_name, "::", input$overview_gate_strategy)
                    rv$experiment_gates[[composite_key]] <- cache_data$gates
                    gates_to_use <- cache_data$gates
                    cat(sprintf("Loaded gates and strategy '%s' from cache\n", input$overview_gate_strategy))
                  }
                  break
                }
              }
            }, error = function(e) {
              # Skip files that can't be read
            })
          }
        }
      }

      # Check if using quadrant strategy
      use_quadrant <- !is.null(gate_strategy$analysis_type) &&
                      gate_strategy$analysis_type == "quadrant_ratio" &&
                      !is.null(gates_to_use$quadrant)

      if(use_quadrant) {
        # Quadrant strategy: show quadrant plots for all Dox+ samples
        plot_quadrant_correlation_overview(exp, gates = gates_to_use, channels = CHANNELS)
      } else {
        # Old strategy: use global Empty_Vector control
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
    }
  })

  # Render sample overview plot (all gates for one sample)
  output$sample_overview_plot <- renderPlot({
    req(rv$experiments, input$sample_overview_experiment,
        input$sample_overview_gate_strategy, input$sample_overview_sample)

    exp <- rv$experiments[[input$sample_overview_experiment]]
    exp_name <- input$sample_overview_experiment
    idx <- as.numeric(input$sample_overview_sample)
    sample_name <- exp$metadata$sample_name[idx]
    fcs <- exp$flowset[[idx]]

    # Use gates for this experiment+strategy combination
    composite_key <- paste0(exp_name, "::", input$sample_overview_gate_strategy)
    gates_to_use <- if(!is.null(rv$experiment_gates[[composite_key]])) {
      rv$experiment_gates[[composite_key]]
    } else {
      GATES
    }

    # Get strategy metadata
    gate_strategy_key <- paste0("GATE_STRATEGY_", input$sample_overview_gate_strategy)
    gate_strategy <- if(!is.null(rv$gate_strategies[[gate_strategy_key]])) {
      rv$gate_strategies[[gate_strategy_key]]
    } else {
      NULL
    }

    # If gate strategy not in memory, try to load from cache
    if(is.null(gate_strategy)) {
      cache_dir <- file.path("analysis_cache", exp_name)
      if(dir.exists(cache_dir)) {
        cache_files <- list.files(cache_dir, pattern = "\\.rds$", full.names = TRUE)
        for(cache_file in cache_files) {
          tryCatch({
            cache_data <- readRDS(cache_file)
            if(!is.null(cache_data$gate_id) && cache_data$gate_id == input$sample_overview_gate_strategy) {
              # Found matching cache file, extract gate strategy AND gates
              if(!is.null(cache_data$gate_strategy)) {
                gate_strategy <- cache_data$gate_strategy
                rv$gate_strategies[[gate_strategy_key]] <- gate_strategy

                # Also load the gates themselves
                if(!is.null(cache_data$gates)) {
                  composite_key <- paste0(exp_name, "::", input$sample_overview_gate_strategy)
                  rv$experiment_gates[[composite_key]] <- cache_data$gates
                  gates_to_use <- cache_data$gates
                  cat(sprintf("Loaded gates and strategy '%s' from cache\n", input$sample_overview_gate_strategy))
                }
                break
              }
            }
          }, error = function(e) {
            # Skip files that can't be read
          })
        }
      }
    }

    # Re-evaluate gates_to_use after cache loading to ensure we have the latest gates
    composite_key <- paste0(exp_name, "::", input$sample_overview_gate_strategy)
    if(!is.null(rv$experiment_gates[[composite_key]])) {
      gates_to_use <- rv$experiment_gates[[composite_key]]
    }

    # Check if using quadrant strategy
    use_quadrant <- !is.null(gate_strategy$analysis_type) &&
                    gate_strategy$analysis_type == "quadrant_ratio" &&
                    !is.null(gates_to_use$quadrant)

    # Calculate thresholds based on strategy
    ha_threshold <- NULL
    edu_threshold <- NULL

    if(use_quadrant) {
      # Quadrant strategy: find paired control
      control_idx <- find_paired_control(sample_name, exp$metadata)
      if(!is.null(control_idx)) {
        control_fcs <- exp$flowset[[control_idx]]
        control_name <- exp$metadata$sample_name[control_idx]

        tryCatch({
          quadrant_result <- calculate_quadrant_from_paired_control(
            control_fcs, fcs, control_name, sample_name,
            gates_to_use, CHANNELS
          )
          ha_threshold <- quadrant_result$ha_threshold
          edu_threshold <- quadrant_result$edu_threshold
        }, error = function(e) {
          cat(sprintf("Error calculating quadrant thresholds: %s\n", e$message))
        })
      }
    } else {
      # Old strategy: use global Empty_Vector control
      control_idx <- find_control_sample(exp$metadata, "Empty_Vector_Dox-")
      if(!is.null(control_idx)) {
        control_fcs <- exp$flowset[[control_idx]]
        control_name <- exp$metadata$sample_name[control_idx]
        control_result <- calculate_ha_threshold_from_control(control_fcs, control_name,
                                                               gates = gates_to_use,
                                                               channels = CHANNELS)
        ha_threshold <- control_result$threshold
      }
    }

    # Set up multi-panel layout
    if(!is.null(ha_threshold)) {
      # 4x2 grid for 8 plots
      par(mfrow = c(4, 2), mar = c(5, 4, 3, 1))

      plot_debris_gate_single(fcs, sample_name, gates = gates_to_use)
      plot_singlet_gate_single(fcs, sample_name, gates = gates_to_use)
      plot_live_gate_single(fcs, sample_name, gates = gates_to_use)
      plot_sphase_outlier_gate_single(fcs, sample_name, gates = gates_to_use)
      plot_fxcycle_quantile_gate_single(fcs, sample_name, gates = gates_to_use)
      plot_edu_fxcycle_gate_single(fcs, sample_name, gates = gates_to_use)

      if(use_quadrant) {
        # Show quadrant plot for Gate 7
        plot_edu_ha_correlation_single(fcs, sample_name, ha_threshold,
                                       gates = gates_to_use, channels = CHANNELS,
                                       show_sample_name = TRUE,
                                       edu_threshold = edu_threshold)
        # Show quadrant plot again (same as Gate 7 for quadrant strategy)
        plot_edu_ha_correlation_single(fcs, sample_name, ha_threshold,
                                       gates = gates_to_use, channels = CHANNELS,
                                       show_sample_name = TRUE,
                                       edu_threshold = edu_threshold)
      } else {
        # Show traditional Gate 7 and correlation plot
        plot_ha_gate_single(fcs, sample_name, ha_threshold, gates = gates_to_use)
        plot_edu_ha_correlation_single(fcs, sample_name, ha_threshold,
                                       gates = gates_to_use, channels = CHANNELS)
      }
    } else {
      # 3x2 grid for 6 plots
      par(mfrow = c(3, 2), mar = c(5, 4, 3, 1))

      plot_debris_gate_single(fcs, sample_name, gates = gates_to_use)
      plot_singlet_gate_single(fcs, sample_name, gates = gates_to_use)
      plot_live_gate_single(fcs, sample_name, gates = gates_to_use)
      plot_sphase_outlier_gate_single(fcs, sample_name, gates = gates_to_use)
      plot_fxcycle_quantile_gate_single(fcs, sample_name, gates = gates_to_use)
      plot_edu_fxcycle_gate_single(fcs, sample_name, gates = gates_to_use)
    }
  })

  # ============================================================================
  # GATING STRATEGY CREATOR
  # ============================================================================

  # Reactive values for creator
  creator_rv <- reactiveValues(
    current_gates = NULL,
    current_strategy = NULL
  )

  # Populate base strategy dropdown
  observe({
    req(rv$available_gate_files)
    updateSelectInput(session, "creator_base_strategy",
                      choices = rv$available_gate_files,
                      selected = "gates_gdef.r")
  })

  # Populate experiment dropdown
  observe({
    req(rv$experiments)
    updateSelectInput(session, "creator_experiment",
                      choices = names(rv$experiments))
  })

  # Update sample selector when experiment changes
  observeEvent(input$creator_experiment, {
    req(input$creator_experiment, rv$experiments)
    exp <- rv$experiments[[input$creator_experiment]]
    if(!is.null(exp)) {
      sample_choices <- setNames(seq_along(exp$metadata$sample_name),
                                  exp$metadata$sample_name)
      updateSelectInput(session, "creator_sample",
                        choices = sample_choices)
    }
  })

  # Load strategy when button clicked
  observeEvent(input$creator_load, {
    req(input$creator_base_strategy)

    gate_file <- input$creator_base_strategy
    gate_path <- file.path("gate_definitions", gate_file)

    if(!file.exists(gate_path)) {
      showNotification("Gate file not found!", type = "error")
      return()
    }

    # Load gates from file
    GATES_env <- new.env()
    tryCatch({
      source(gate_path, local = GATES_env)
      creator_rv$current_gates <- GATES_env$GATES
      creator_rv$current_strategy <- GATES_env$GATE_STRATEGY

      showNotification("Strategy loaded successfully!", type = "message", duration = 2)
    }, error = function(e) {
      showNotification(sprintf("Error loading strategy: %s", e$message), type = "error")
    })
  })

  # Generate edit UI for matrix gates
  generate_matrix_ui <- function(gate_matrix, gate_name) {
    if(is.null(gate_matrix)) return(NULL)

    n_vertices <- nrow(gate_matrix)
    col_names <- colnames(gate_matrix)

    rows <- lapply(1:n_vertices, function(i) {
      fluidRow(
        column(1, p(sprintf("V%d:", i), style = "margin-top: 5px;")),
        column(5,
               numericInput(paste0("creator_", gate_name, "_x", i),
                           label = col_names[1],
                           value = gate_matrix[i, 1],
                           width = "100%")
        ),
        column(5,
               numericInput(paste0("creator_", gate_name, "_y", i),
                           label = col_names[2],
                           value = gate_matrix[i, 2],
                           width = "100%")
        )
      )
    })

    div(rows)
  }

  # Render UI for each gate type
  output$creator_debris_ui <- renderUI({
    req(creator_rv$current_gates)
    generate_matrix_ui(creator_rv$current_gates$debris, "debris")
  })

  output$creator_singlet_ui <- renderUI({
    req(creator_rv$current_gates)
    generate_matrix_ui(creator_rv$current_gates$singlet, "singlet")
  })

  output$creator_live_ui <- renderUI({
    req(creator_rv$current_gates)
    generate_matrix_ui(creator_rv$current_gates$live_cells, "live")
  })

  output$creator_sphase_ui <- renderUI({
    req(creator_rv$current_gates)
    generate_matrix_ui(creator_rv$current_gates$s_phase_outliers, "sphase")
  })

  output$creator_fxcycle_ui <- renderUI({
    req(creator_rv$current_gates)
    gate <- creator_rv$current_gates$fxcycle_quantile

    div(
      numericInput("creator_fxcycle_prob_low", "Lower Percentile:",
                   value = gate$probs[1], min = 0, max = 1, step = 0.01),
      numericInput("creator_fxcycle_prob_high", "Upper Percentile:",
                   value = gate$probs[2], min = 0, max = 1, step = 0.01),
      textInput("creator_fxcycle_param", "Parameter:",
                value = gate$parameter),
      textInput("creator_fxcycle_desc", "Description:",
                value = gate$description)
    )
  })

  output$creator_edu_fxcycle_ui <- renderUI({
    req(creator_rv$current_gates)
    gate <- creator_rv$current_gates$edu_fxcycle_sphase

    div(
      numericInput("creator_edu_prob", "EdU Percentile (top %):",
                   value = gate$edu_prob, min = 0, max = 1, step = 0.01),
      textInput("creator_edu_param", "EdU Parameter:",
                value = gate$edu_parameter),
      numericInput("creator_edu_fxcycle_prob_low", "FxCycle Lower Percentile:",
                   value = gate$fxcycle_probs[1], min = 0, max = 1, step = 0.01),
      numericInput("creator_edu_fxcycle_prob_high", "FxCycle Upper Percentile:",
                   value = gate$fxcycle_probs[2], min = 0, max = 1, step = 0.01),
      textInput("creator_edu_fxcycle_param", "FxCycle Parameter:",
                value = gate$fxcycle_parameter),
      textInput("creator_edu_fxcycle_desc", "Description:",
                value = gate$description)
    )
  })

  output$creator_gate7_header <- renderUI({
    req(creator_rv$current_gates)
    is_quadrant <- !is.null(creator_rv$current_gates$quadrant)
    header_text <- if(is_quadrant) "Gate 7: Quadrant" else "Gate 7: HA Positive"

    tags$a(href = "#collapse_ha", `data-toggle` = "collapse",
           h5(HTML(paste0("<i class='glyphicon glyphicon-chevron-down'></i> ", header_text)),
              style = "color: #337ab7; cursor: pointer;"))
  })

  output$creator_ha_ui <- renderUI({
    req(creator_rv$current_gates)

    # Check if this is a quadrant gate or ha_positive gate
    is_quadrant <- !is.null(creator_rv$current_gates$quadrant)
    gate <- if(is_quadrant) creator_rv$current_gates$quadrant else creator_rv$current_gates$ha_positive

    if(is_quadrant) {
      # Determine quadrant type
      quad_type <- if(!is.null(gate$type)) gate$type else "paired_control_quadrant"
      is_control_based <- (quad_type == "control_based_quadrant")

      # Quadrant gate UI (dual thresholds for HA and EdU)
      div(
        h5("Quadrant Gating (HA vs EdU)", style = "margin-top: 0;"),
        radioButtons("creator_quad_type", "Quadrant Type:",
                     choices = c("Paired Control (per-sample thresholds)" = "paired_control_quadrant",
                                 "Control-Based (fixed thresholds)" = "control_based_quadrant"),
                     selected = quad_type),
        numericInput("creator_quad_ha_prob", "HA Control Percentile:",
                     value = if(is.null(gate$ha_prob)) 0.98 else gate$ha_prob,
                     min = 0, max = 1, step = 0.01),
        textInput("creator_quad_ha_param", "HA Parameter:",
                  value = if(is.null(gate$ha_parameter)) "FL3-A" else gate$ha_parameter),
        numericInput("creator_quad_edu_prob", "EdU Control Percentile:",
                     value = if(is.null(gate$edu_prob)) 0.02 else gate$edu_prob,
                     min = 0, max = 1, step = 0.01),
        textInput("creator_quad_edu_param", "EdU Parameter:",
                  value = if(is.null(gate$edu_parameter)) "FL1-A" else gate$edu_parameter),
        conditionalPanel(
          condition = "input.creator_quad_type == 'control_based_quadrant'",
          textInput("creator_quad_control_pattern", "Control Sample Pattern:",
                    value = if(is.null(gate$control_pattern)) "121_Empty_Vector_Dox-" else gate$control_pattern)
        ),
        conditionalPanel(
          condition = "input.creator_quad_type == 'paired_control_quadrant'",
          textInput("creator_quad_control_suffix", "Control Suffix:",
                    value = if(is.null(gate$control_suffix)) "Dox-" else gate$control_suffix),
          textInput("creator_quad_test_suffix", "Test Suffix:",
                    value = if(is.null(gate$test_suffix)) "Dox+" else gate$test_suffix)
        ),
        textInput("creator_quad_desc", "Description:",
                  value = if(is.null(gate$description)) {
                    if(is_control_based) "Quadrant gating using control-based thresholds" else "Quadrant gating using paired Dox- control"
                  } else gate$description)
      )
    } else {
      # Standard HA positive gate UI
      div(
        numericInput("creator_ha_prob", "Control Percentile:",
                     value = if(is.null(gate$prob)) 0.98 else gate$prob,
                     min = 0, max = 1, step = 0.01),
        textInput("creator_ha_param", "Parameter:",
                  value = if(is.null(gate$parameter)) "FL3-A" else gate$parameter),
        textInput("creator_ha_control", "Control Pattern:",
                  value = if(is.null(gate$control_pattern)) "Empty_Vector_Dox-" else gate$control_pattern),
        textInput("creator_ha_desc", "Description:",
                  value = if(is.null(gate$description)) "HA positive threshold" else gate$description)
      )
    }
  })

  # Get currently edited gates (read from UI inputs)
  get_edited_gates <- reactive({
    req(creator_rv$current_gates)

    gates <- list()

    # Helper function to safely get input value with default fallback
    safe_input <- function(input_id, default_val) {
      val <- input[[input_id]]
      if(is.null(val)) return(default_val)
      return(val)
    }

    # Debris gate
    if(!is.null(creator_rv$current_gates$debris)) {
      n <- nrow(creator_rv$current_gates$debris)
      coords <- matrix(nrow = n, ncol = 2)
      for(i in 1:n) {
        coords[i, 1] <- safe_input(paste0("creator_debris_x", i),
                                   creator_rv$current_gates$debris[i, 1])
        coords[i, 2] <- safe_input(paste0("creator_debris_y", i),
                                   creator_rv$current_gates$debris[i, 2])
      }
      colnames(coords) <- colnames(creator_rv$current_gates$debris)
      gates$debris <- coords
    }

    # Singlet gate
    if(!is.null(creator_rv$current_gates$singlet)) {
      n <- nrow(creator_rv$current_gates$singlet)
      coords <- matrix(nrow = n, ncol = 2)
      for(i in 1:n) {
        coords[i, 1] <- safe_input(paste0("creator_singlet_x", i),
                                   creator_rv$current_gates$singlet[i, 1])
        coords[i, 2] <- safe_input(paste0("creator_singlet_y", i),
                                   creator_rv$current_gates$singlet[i, 2])
      }
      colnames(coords) <- colnames(creator_rv$current_gates$singlet)
      gates$singlet <- coords
    }

    # Live cells gate
    if(!is.null(creator_rv$current_gates$live_cells)) {
      n <- nrow(creator_rv$current_gates$live_cells)
      coords <- matrix(nrow = n, ncol = 2)
      for(i in 1:n) {
        coords[i, 1] <- safe_input(paste0("creator_live_x", i),
                                   creator_rv$current_gates$live_cells[i, 1])
        coords[i, 2] <- safe_input(paste0("creator_live_y", i),
                                   creator_rv$current_gates$live_cells[i, 2])
      }
      colnames(coords) <- colnames(creator_rv$current_gates$live_cells)
      gates$live_cells <- coords
    }

    # S-phase outliers gate
    if(!is.null(creator_rv$current_gates$s_phase_outliers)) {
      n <- nrow(creator_rv$current_gates$s_phase_outliers)
      coords <- matrix(nrow = n, ncol = 2)
      for(i in 1:n) {
        coords[i, 1] <- safe_input(paste0("creator_sphase_x", i),
                                   creator_rv$current_gates$s_phase_outliers[i, 1])
        coords[i, 2] <- safe_input(paste0("creator_sphase_y", i),
                                   creator_rv$current_gates$s_phase_outliers[i, 2])
      }
      colnames(coords) <- colnames(creator_rv$current_gates$s_phase_outliers)
      gates$s_phase_outliers <- coords
    }

    # FxCycle quantile gate
    fxcycle_gate <- creator_rv$current_gates$fxcycle_quantile
    gates$fxcycle_quantile <- list(
      type = "quantile_range",
      parameter = safe_input("creator_fxcycle_param", fxcycle_gate$parameter),
      probs = c(safe_input("creator_fxcycle_prob_low", fxcycle_gate$probs[1]),
                safe_input("creator_fxcycle_prob_high", fxcycle_gate$probs[2])),
      description = safe_input("creator_fxcycle_desc", fxcycle_gate$description)
    )

    # EdU + FxCycle gate
    edu_gate <- creator_rv$current_gates$edu_fxcycle_sphase
    gates$edu_fxcycle_sphase <- list(
      type = "dual_quantile",
      edu_parameter = safe_input("creator_edu_param", edu_gate$edu_parameter),
      edu_prob = safe_input("creator_edu_prob", edu_gate$edu_prob),
      fxcycle_parameter = safe_input("creator_edu_fxcycle_param", edu_gate$fxcycle_parameter),
      fxcycle_probs = c(safe_input("creator_edu_fxcycle_prob_low", edu_gate$fxcycle_probs[1]),
                        safe_input("creator_edu_fxcycle_prob_high", edu_gate$fxcycle_probs[2])),
      description = safe_input("creator_edu_fxcycle_desc", edu_gate$description)
    )

    # Gate 7: HA positive or Quadrant gate
    if(!is.null(creator_rv$current_gates$quadrant)) {
      # Quadrant gate
      quad_gate <- creator_rv$current_gates$quadrant
      quad_type <- safe_input("creator_quad_type", quad_gate$type)

      # Build quadrant gate based on type
      quad_list <- list(
        type = quad_type,
        ha_parameter = safe_input("creator_quad_ha_param", quad_gate$ha_parameter),
        edu_parameter = safe_input("creator_quad_edu_param", quad_gate$edu_parameter),
        ha_prob = safe_input("creator_quad_ha_prob", quad_gate$ha_prob),
        edu_prob = safe_input("creator_quad_edu_prob", quad_gate$edu_prob),
        description = safe_input("creator_quad_desc", quad_gate$description)
      )

      # Add type-specific fields
      if(quad_type == "control_based_quadrant") {
        quad_list$control_pattern <- safe_input("creator_quad_control_pattern", quad_gate$control_pattern)
      } else {
        quad_list$control_suffix <- safe_input("creator_quad_control_suffix", quad_gate$control_suffix)
        quad_list$test_suffix <- safe_input("creator_quad_test_suffix", quad_gate$test_suffix)
      }

      gates$quadrant <- quad_list
    } else {
      # Standard HA positive gate
      ha_gate <- creator_rv$current_gates$ha_positive
      gates$ha_positive <- list(
        type = "control_based_threshold",
        parameter = safe_input("creator_ha_param", ha_gate$parameter),
        control_pattern = safe_input("creator_ha_control", ha_gate$control_pattern),
        prob = safe_input("creator_ha_prob", ha_gate$prob),
        description = safe_input("creator_ha_desc", ha_gate$description)
      )
    }

    gates
  })

  # Render sample name title for preview plot
  output$creator_preview_sample_title <- renderText({
    req(rv$experiments, input$creator_experiment, input$creator_sample)
    exp <- rv$experiments[[input$creator_experiment]]
    idx <- as.numeric(input$creator_sample)
    exp$metadata$sample_name[idx]
  })

  # Render preview plot with edited gates
  output$creator_preview_plot <- renderPlot({
    req(rv$experiments, input$creator_experiment, input$creator_sample)
    req(creator_rv$current_gates)

    exp <- rv$experiments[[input$creator_experiment]]
    idx <- as.numeric(input$creator_sample)
    sample_name <- exp$metadata$sample_name[idx]
    fcs <- exp$flowset[[idx]]

    # Get edited gates
    gates_to_use <- get_edited_gates()

    # Check if using quadrant strategy
    use_quadrant <- !is.null(creator_rv$current_strategy$analysis_type) &&
                    creator_rv$current_strategy$analysis_type == "quadrant_ratio" &&
                    !is.null(gates_to_use$quadrant)

    # Calculate thresholds based on strategy
    ha_threshold <- NULL
    edu_threshold <- NULL

    if(use_quadrant) {
      # Quadrant strategy: find paired control
      control_idx <- find_paired_control(sample_name, exp$metadata)
      if(!is.null(control_idx)) {
        control_fcs <- exp$flowset[[control_idx]]
        control_name <- exp$metadata$sample_name[control_idx]

        tryCatch({
          quadrant_result <- calculate_quadrant_from_paired_control(
            control_fcs, fcs, control_name, sample_name,
            gates_to_use, CHANNELS
          )
          ha_threshold <- quadrant_result$ha_threshold
          edu_threshold <- quadrant_result$edu_threshold
        }, error = function(e) {
          cat(sprintf("Error calculating quadrant thresholds: %s\n", e$message))
        })
      }
    } else {
      # Old strategy: use global Empty_Vector control
      control_idx <- find_control_sample(exp$metadata, "Empty_Vector_Dox-")
      if(!is.null(control_idx)) {
        control_fcs <- exp$flowset[[control_idx]]
        control_name <- exp$metadata$sample_name[control_idx]
        control_result <- calculate_ha_threshold_from_control(control_fcs, control_name,
                                                               gates = gates_to_use,
                                                               channels = CHANNELS)
        ha_threshold <- control_result$threshold
      }
    }

    # Set up multi-panel layout
    if(!is.null(ha_threshold)) {
      par(mfrow = c(4, 2), mar = c(5, 4, 3, 1))

      plot_debris_gate_single(fcs, sample_name, gates = gates_to_use, show_sample_name = FALSE)
      plot_singlet_gate_single(fcs, sample_name, gates = gates_to_use, show_sample_name = FALSE)
      plot_live_gate_single(fcs, sample_name, gates = gates_to_use, show_sample_name = FALSE)
      plot_sphase_outlier_gate_single(fcs, sample_name, gates = gates_to_use, show_sample_name = FALSE)
      plot_fxcycle_quantile_gate_single(fcs, sample_name, gates = gates_to_use, show_sample_name = FALSE)
      plot_edu_fxcycle_gate_single(fcs, sample_name, gates = gates_to_use, show_sample_name = FALSE)

      if(use_quadrant) {
        # Show quadrant plot for Gate 7
        plot_edu_ha_correlation_single(fcs, sample_name, ha_threshold,
                                       gates = gates_to_use, channels = CHANNELS,
                                       show_sample_name = FALSE,
                                       edu_threshold = edu_threshold)
        # Show quadrant plot again (same as Gate 7 for quadrant strategy)
        plot_edu_ha_correlation_single(fcs, sample_name, ha_threshold,
                                       gates = gates_to_use, channels = CHANNELS,
                                       show_sample_name = FALSE,
                                       edu_threshold = edu_threshold)
      } else {
        # Show traditional Gate 7 and correlation plot
        plot_ha_gate_single(fcs, sample_name, ha_threshold, gates = gates_to_use, show_sample_name = FALSE)
        plot_edu_ha_correlation_single(fcs, sample_name, ha_threshold,
                                       gates = gates_to_use, channels = CHANNELS, show_sample_name = FALSE)
      }
    } else {
      par(mfrow = c(3, 2), mar = c(5, 4, 3, 1))

      plot_debris_gate_single(fcs, sample_name, gates = gates_to_use, show_sample_name = FALSE)
      plot_singlet_gate_single(fcs, sample_name, gates = gates_to_use, show_sample_name = FALSE)
      plot_live_gate_single(fcs, sample_name, gates = gates_to_use, show_sample_name = FALSE)
      plot_sphase_outlier_gate_single(fcs, sample_name, gates = gates_to_use, show_sample_name = FALSE)
      plot_fxcycle_quantile_gate_single(fcs, sample_name, gates = gates_to_use, show_sample_name = FALSE)
      plot_edu_fxcycle_gate_single(fcs, sample_name, gates = gates_to_use, show_sample_name = FALSE)
    }
  })

  # Render individual plot based on user selection
  output$creator_single_plot_output <- renderPlot({
    req(rv$experiments, input$creator_experiment, input$creator_sample)
    req(creator_rv$current_gates, input$creator_single_plot)

    exp <- rv$experiments[[input$creator_experiment]]
    idx <- as.numeric(input$creator_sample)
    sample_name <- exp$metadata$sample_name[idx]
    fcs <- exp$flowset[[idx]]

    # Get edited gates
    gates_to_use <- get_edited_gates()

    # Check if using quadrant strategy
    use_quadrant <- !is.null(creator_rv$current_strategy$analysis_type) &&
                    creator_rv$current_strategy$analysis_type == "quadrant_ratio" &&
                    !is.null(gates_to_use$quadrant)

    # Calculate thresholds based on strategy
    ha_threshold <- NULL
    edu_threshold <- NULL

    if(use_quadrant) {
      # Quadrant strategy: find paired control
      control_idx <- find_paired_control(sample_name, exp$metadata)
      if(!is.null(control_idx)) {
        control_fcs <- exp$flowset[[control_idx]]
        control_name <- exp$metadata$sample_name[control_idx]

        tryCatch({
          quadrant_result <- calculate_quadrant_from_paired_control(
            control_fcs, fcs, control_name, sample_name,
            gates_to_use, CHANNELS
          )
          ha_threshold <- quadrant_result$ha_threshold
          edu_threshold <- quadrant_result$edu_threshold
        }, error = function(e) {
          cat(sprintf("Error calculating quadrant thresholds: %s\n", e$message))
        })
      }
    } else {
      # Old strategy: use global Empty_Vector control
      control_idx <- find_control_sample(exp$metadata, "Empty_Vector_Dox-")
      if(!is.null(control_idx)) {
        control_fcs <- exp$flowset[[control_idx]]
        control_name <- exp$metadata$sample_name[control_idx]
        control_result <- calculate_ha_threshold_from_control(control_fcs, control_name,
                                                               gates = gates_to_use,
                                                               channels = CHANNELS)
        ha_threshold <- control_result$threshold
      }
    }

    # Render the selected plot
    switch(input$creator_single_plot,
           "gate1" = plot_debris_gate_single(fcs, sample_name, gates = gates_to_use),
           "gate2" = plot_singlet_gate_single(fcs, sample_name, gates = gates_to_use),
           "gate3" = plot_live_gate_single(fcs, sample_name, gates = gates_to_use),
           "gate4" = plot_sphase_outlier_gate_single(fcs, sample_name, gates = gates_to_use),
           "gate5" = plot_fxcycle_quantile_gate_single(fcs, sample_name, gates = gates_to_use),
           "gate6" = plot_edu_fxcycle_gate_single(fcs, sample_name, gates = gates_to_use),
           "gate7" = {
             if(!is.null(ha_threshold)) {
               if(use_quadrant) {
                 # Show quadrant plot for Gate 7
                 plot_edu_ha_correlation_single(fcs, sample_name, ha_threshold,
                                                gates = gates_to_use, channels = CHANNELS,
                                                edu_threshold = edu_threshold)
               } else {
                 # Show traditional HA gate
                 plot_ha_gate_single(fcs, sample_name, ha_threshold, gates = gates_to_use)
               }
             } else {
               plot.new()
               text(0.5, 0.5, "No control sample found\n(needed to calculate thresholds)",
                    cex = 1.5, col = "red")
             }
           },
           "correlation" = {
             if(!is.null(ha_threshold)) {
               if(use_quadrant) {
                 # Show quadrant plot with thresholds
                 plot_edu_ha_correlation_single(fcs, sample_name, ha_threshold,
                                                gates = gates_to_use, channels = CHANNELS,
                                                edu_threshold = edu_threshold)
               } else {
                 # Show traditional correlation plot
                 plot_edu_ha_correlation_single(fcs, sample_name, ha_threshold,
                                                gates = gates_to_use, channels = CHANNELS)
               }
             } else {
               plot.new()
               text(0.5, 0.5, "No control sample found\n(needed to calculate thresholds)",
                    cex = 1.5, col = "red")
             }
           }
    )
  })

  # Save new strategy
  observeEvent(input$creator_save, {
    req(input$creator_new_id, input$creator_new_name)

    if(input$creator_new_id == "" || input$creator_new_name == "") {
      showNotification("Please provide both Strategy ID and Name", type = "warning")
      return()
    }

    # Get edited gates
    gates_to_save <- get_edited_gates()

    # Determine if this is a quadrant or ha_positive strategy
    is_quadrant <- !is.null(gates_to_save$quadrant)

    # Build Gate 7 section based on type
    gate7_section <- if(is_quadrant) {
      quad_type <- gates_to_save$quadrant$type
      if(quad_type == "control_based_quadrant") {
        # Control-based quadrant
        sprintf('# Gate 7: Quadrant gating with control-based thresholds
GATES$quadrant <- list(
  type = "control_based_quadrant",
  ha_parameter = "%s",
  edu_parameter = "%s",
  control_pattern = "%s",
  ha_prob = %g,
  edu_prob = %g,
  description = "%s"
)',
          gates_to_save$quadrant$ha_parameter,
          gates_to_save$quadrant$edu_parameter,
          gates_to_save$quadrant$control_pattern,
          gates_to_save$quadrant$ha_prob,
          gates_to_save$quadrant$edu_prob,
          gates_to_save$quadrant$description)
      } else {
        # Paired control quadrant
        sprintf('# Gate 7: Quadrant gating (HA vs EdU)
GATES$quadrant <- list(
  type = "paired_control_quadrant",
  ha_parameter = "%s",
  edu_parameter = "%s",
  control_suffix = "%s",
  test_suffix = "%s",
  ha_prob = %g,
  edu_prob = %g,
  description = "%s"
)',
          gates_to_save$quadrant$ha_parameter,
          gates_to_save$quadrant$edu_parameter,
          gates_to_save$quadrant$control_suffix,
          gates_to_save$quadrant$test_suffix,
          gates_to_save$quadrant$ha_prob,
          gates_to_save$quadrant$edu_prob,
          gates_to_save$quadrant$description)
      }
    } else {
      sprintf('# Gate 7: HA positive
GATES$ha_positive <- list(
  type = "control_based_threshold",
  parameter = "%s",
  control_pattern = "%s",
  prob = %g,
  description = "%s"
)',
        gates_to_save$ha_positive$parameter,
        gates_to_save$ha_positive$control_pattern,
        gates_to_save$ha_positive$prob,
        gates_to_save$ha_positive$description)
    }

    # Add analysis_type to metadata for quadrant strategies
    metadata_extra <- if(is_quadrant) {
      ',\n  analysis_type = "quadrant_ratio"  # Flag to indicate this uses ratio instead of correlation'
    } else {
      ''
    }

    # Create file content
    file_content <- sprintf('# ==============================================================================
# GATE DEFINITIONS - %s
# ==============================================================================

GATES <- list()

# Gate 1: Debris removal (FSC-A vs SSC-A)
GATES$debris <- matrix(c(
%s
), ncol = 2, byrow = TRUE)
colnames(GATES$debris) <- c("FSC-A", "SSC-A")

# Gate 2: Singlets (FSC-A vs FSC-H)
GATES$singlet <- matrix(c(
%s
), ncol = 2, byrow = TRUE)
colnames(GATES$singlet) <- c("FSC-A", "FSC-H")

# Gate 3: Live cells (DCM-A vs SSC-A)
GATES$live_cells <- matrix(c(
%s
), ncol = 2, byrow = TRUE)
colnames(GATES$live_cells) <- c("DCM-A", "SSC-A")

# Gate 4: S-phase outliers (FxCycle-A vs EdU-A)
GATES$s_phase_outliers <- matrix(c(
%s
), ncol = 2, byrow = TRUE)
colnames(GATES$s_phase_outliers) <- c("FxCycle-A", "EdU-A")

# Gate 5: FxCycle quantile
GATES$fxcycle_quantile <- list(
  type = "quantile_range",
  parameter = "%s",
  probs = c(%g, %g),
  description = "%s"
)

# Gate 6: Top %.0f%%%% EdU + FxCycle range
GATES$edu_fxcycle_sphase <- list(
  type = "dual_quantile",
  edu_parameter = "%s",
  edu_prob = %g,
  fxcycle_parameter = "%s",
  fxcycle_probs = c(%g, %g),
  description = "%s"
)

%s

# Gate strategy metadata
GATE_STRATEGY <- list(
  id = "%s",
  name = "%s",
  description = "%s",
  created = "%s"%s
)
',
      input$creator_new_id,
      # Debris coordinates
      paste(apply(gates_to_save$debris, 1, function(row) sprintf("  %g, %g", row[1], row[2])), collapse = ",\n"),
      # Singlet coordinates
      paste(apply(gates_to_save$singlet, 1, function(row) sprintf("  %g, %g", row[1], row[2])), collapse = ",\n"),
      # Live cells coordinates
      paste(apply(gates_to_save$live_cells, 1, function(row) sprintf("  %g, %g", row[1], row[2])), collapse = ",\n"),
      # S-phase outliers coordinates
      paste(apply(gates_to_save$s_phase_outliers, 1, function(row) sprintf("  %g, %g", row[1], row[2])), collapse = ",\n"),
      # FxCycle quantile
      gates_to_save$fxcycle_quantile$parameter,
      gates_to_save$fxcycle_quantile$probs[1],
      gates_to_save$fxcycle_quantile$probs[2],
      gates_to_save$fxcycle_quantile$description,
      # EdU + FxCycle
      gates_to_save$edu_fxcycle_sphase$edu_prob * 100,
      gates_to_save$edu_fxcycle_sphase$edu_parameter,
      gates_to_save$edu_fxcycle_sphase$edu_prob,
      gates_to_save$edu_fxcycle_sphase$fxcycle_parameter,
      gates_to_save$edu_fxcycle_sphase$fxcycle_probs[1],
      gates_to_save$edu_fxcycle_sphase$fxcycle_probs[2],
      gates_to_save$edu_fxcycle_sphase$description,
      # Gate 7
      gate7_section,
      # Metadata
      input$creator_new_id,
      input$creator_new_name,
      input$creator_new_desc,
      as.character(Sys.time()),
      metadata_extra
    )

    # Save to file
    file_name <- paste0("gates_", input$creator_new_id, ".r")
    file_path <- file.path("gate_definitions", file_name)

    tryCatch({
      writeLines(file_content, file_path)
      showNotification(sprintf("Strategy saved as %s", file_name), type = "message", duration = 5)

      # Rescan gate files
      rv$available_gate_files <- scan_gate_files()
      updateSelectInput(session, "creator_base_strategy",
                        choices = rv$available_gate_files,
                        selected = file_name)
    }, error = function(e) {
      showNotification(sprintf("Error saving file: %s", e$message), type = "error")
    })
  })

  # ============================================================================
  # END GATING STRATEGY CREATOR
  # ============================================================================

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

    # Get strategy metadata
    gate_strategy_key <- paste0("GATE_STRATEGY_", input$browse_gate_strategy)
    gate_strategy <- if(!is.null(rv$gate_strategies[[gate_strategy_key]])) {
      rv$gate_strategies[[gate_strategy_key]]
    } else {
      NULL
    }

    # Check if using quadrant strategy
    use_quadrant <- !is.null(gate_strategy$analysis_type) &&
                    gate_strategy$analysis_type == "quadrant_ratio" &&
                    !is.null(gates_to_use$quadrant)

    sample_name <- exp$metadata$sample_name[idx]

    # Debug output
    cat(sprintf("\n=== Gate 7 Plot Debug ===\n"))
    cat(sprintf("Sample: %s\n", sample_name))
    cat(sprintf("use_quadrant: %s\n", use_quadrant))
    cat(sprintf("gate_strategy: %s\n", if(!is.null(gate_strategy)) "PRESENT" else "NULL"))
    if(!is.null(gate_strategy)) {
      cat(sprintf("  analysis_type: %s\n", gate_strategy$analysis_type))
    }
    cat(sprintf("gates$quadrant: %s\n", if(!is.null(gates_to_use$quadrant)) "PRESENT" else "NULL"))

    if(use_quadrant) {
      # Quadrant strategy: find paired control and show quadrant plot
      control_idx <- find_paired_control(sample_name, exp$metadata)

      if(is.null(control_idx)) {
        plot.new()
        text(0.5, 0.5, sprintf("No paired Dox- control found for:\n%s", sample_name),
             cex = 1.2, col = "red")
        return()
      }

      # Calculate quadrant thresholds from paired control
      control_fcs <- exp$flowset[[control_idx]]
      test_fcs <- exp$flowset[[idx]]
      control_name <- exp$metadata$sample_name[control_idx]

      tryCatch({
        quadrant_result <- calculate_quadrant_from_paired_control(
          control_fcs, test_fcs, control_name, sample_name,
          gates_to_use, CHANNELS
        )

        # Show quadrant plot using the correlation plotting function with edu_threshold
        plot_edu_ha_correlation_single(test_fcs, sample_name,
                                       ha_threshold = quadrant_result$ha_threshold,
                                       gates = gates_to_use,
                                       channels = CHANNELS,
                                       show_sample_name = TRUE,
                                       edu_threshold = quadrant_result$edu_threshold)
      }, error = function(e) {
        plot.new()
        text(0.5, 0.5, sprintf("Error in quadrant analysis:\n%s", e$message),
             cex = 1.2, col = "red")
      })
    } else {
      # Old strategy: use global Empty_Vector control
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

      plot_ha_gate_single(exp$flowset[[idx]], sample_name, ha_threshold,
                          gates = gates_to_use)
    }
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

      # Create export data with summary statistics
      export_data <- rv$all_results
      export_data$Avg_All <- ""
      export_data$SD_All <- ""
      export_data$Avg_Week <- ""
      export_data$SD_Week <- ""

      # Calculate summary statistics per Cell_line
      if("Correlation" %in% names(rv$all_results) && "Cell_line" %in% names(rv$all_results)) {
        numeric_corr <- as.numeric(rv$all_results$Correlation)
        valid_data <- rv$all_results[!is.na(numeric_corr), ]
        valid_data$Correlation_num <- as.numeric(valid_data$Correlation)

        if(nrow(valid_data) > 0) {
          # Calculate stats per cell line (all experiments)
          all_stats <- valid_data %>%
            group_by(Cell_line) %>%
            summarize(
              Avg_All = mean(Correlation_num, na.rm = TRUE),
              SD_All = ifelse(n() > 1, sd(Correlation_num, na.rm = TRUE), NA),
              n_all = n(),
              .groups = 'drop'
            )

          # Calculate stats per cell line and week
          week_stats <- NULL
          if("Year_Week" %in% names(valid_data)) {
            week_stats <- valid_data %>%
              group_by(Cell_line, Year_Week) %>%
              summarize(
                Avg_Week = mean(Correlation_num, na.rm = TRUE),
                SD_Week = ifelse(n() > 1, sd(Correlation_num, na.rm = TRUE), NA),
                n_week = n(),
                .groups = 'drop'
              )
          }

          # Populate per-row stats
          for(i in seq_len(nrow(export_data))) {
            cell_line <- export_data$Cell_line[i]
            # All experiments stats
            cell_stats <- all_stats[all_stats$Cell_line == cell_line, ]
            if(nrow(cell_stats) > 0) {
              export_data$Avg_All[i] <- sprintf("%.4f", cell_stats$Avg_All[1])
              if(!is.na(cell_stats$SD_All[1])) {
                export_data$SD_All[i] <- sprintf("%.4f", cell_stats$SD_All[1])
              }
            }
            # Same week stats
            if(!is.null(week_stats) && "Year_Week" %in% names(export_data)) {
              year_week <- export_data$Year_Week[i]
              wk_stats <- week_stats[week_stats$Cell_line == cell_line & week_stats$Year_Week == year_week, ]
              if(nrow(wk_stats) > 0) {
                export_data$Avg_Week[i] <- sprintf("%.4f", wk_stats$Avg_Week[1])
                if(!is.na(wk_stats$SD_Week[1])) {
                  export_data$SD_Week[i] <- sprintf("%.4f", wk_stats$SD_Week[1])
                }
              }
            }
          }
        }
      }

      write.xlsx(export_data, file)
    }
  )
  
  # Multi-sample comparison
  
  # Create reactive value to store selected samples for comparison
  rv$comparison_samples <- reactiveVal(data.frame())
  
  # Display sample selector table (only analyzed samples)
  output$comparison_sample_selector <- renderDT({
    req(rv$all_results)

    # React to ui_refresh_trigger to update after analysis
    rv$ui_refresh_trigger

    # Debug: print what we have
    cat(sprintf("\n=== Multi-Sample Comparison Debug ===\n"))
    cat(sprintf("Total rows in rv$all_results: %d\n", nrow(rv$all_results)))
    cat(sprintf("Gate IDs present: %s\n", paste(unique(rv$all_results$Gate_ID), collapse = ", ")))
    cat(sprintf("Rows with non-NA Correlation: %d\n", sum(!is.na(rv$all_results$Correlation))))

    # Check quadrant specifically
    quadrant_rows <- rv$all_results[!is.na(rv$all_results$Gate_ID) & rv$all_results$Gate_ID == "quadrant", ]
    cat(sprintf("Quadrant rows total: %d\n", nrow(quadrant_rows)))
    cat(sprintf("Quadrant rows with Correlation: %d\n", sum(!is.na(quadrant_rows$Correlation))))
    if(nrow(quadrant_rows) > 0) {
      cat(sprintf("Sample quadrant row - Correlation: %s, HA_Pos_Pct: %s, Strength_Ratio: %s\n",
                  quadrant_rows$Correlation[1],
                  quadrant_rows$HA_Pos_Pct[1],
                  quadrant_rows$Strength_Ratio[1]))
    }

    # Start with all results
    analyzed <- rv$all_results

    # Filter out Dox- samples unless checkbox is checked
    if(!isTRUE(input$msc_show_dox_minus)) {
      analyzed <- analyzed[!grepl("Dox-", analyzed$Sample, ignore.case = TRUE), ]
    }

    # Filter out unanalyzed samples unless checkbox is checked
    if(!isTRUE(input$msc_show_unanalyzed)) {
      analyzed <- analyzed[!is.na(analyzed$Correlation) & analyzed$Correlation != "Not analyzed", ]
    }

    cat(sprintf("After filtering, rows: %d\n", nrow(analyzed)))
    quadrant_analyzed <- analyzed[!is.na(analyzed$Gate_ID) & analyzed$Gate_ID == "quadrant", ]
    cat(sprintf("Quadrant in filtered: %d\n", nrow(quadrant_analyzed)))

    if(nrow(analyzed) == 0) {
      return(datatable(data.frame(Message = "No samples match current filters")))
    }

    # Show relevant columns (include Date and Year_Week if available)
    cols_to_show <- c("Experiment")
    if("Date" %in% names(analyzed)) cols_to_show <- c(cols_to_show, "Date")
    if("Year_Week" %in% names(analyzed)) cols_to_show <- c(cols_to_show, "Year_Week")
    cols_to_show <- c(cols_to_show, "Sample", "Cell_line", "Gene",
                      "Mutation", "Correlation")

    # Add HA_Pos_Pct if it exists
    if("HA_Pos_Pct" %in% names(analyzed)) {
      cols_to_show <- c(cols_to_show, "HA_Pos_Pct")
      cat(sprintf("HA_Pos_Pct column exists - non-NA values: %d\n", sum(!is.na(analyzed$HA_Pos_Pct))))
    } else {
      cat("HA_Pos_Pct column does NOT exist in analyzed\n")
    }

    # Add Strength_Ratio if it exists
    if("Strength_Ratio" %in% names(analyzed)) {
      cols_to_show <- c(cols_to_show, "Strength_Ratio")
      cat(sprintf("Strength_Ratio column exists - non-NA values: %d\n", sum(!is.na(analyzed$Strength_Ratio))))
    } else {
      cat("Strength_Ratio column does NOT exist in analyzed\n")
    }

    cols_to_show <- c(cols_to_show, "N_cells", "Gate_ID")

    # Add Notes if it exists
    if("Notes" %in% names(analyzed)) {
      cols_to_show <- c(cols_to_show, "Notes")
    }
    display_data <- analyzed[, cols_to_show]

    # Calculate cell line statistics (across all experiments and same-week)
    if("Cell_line" %in% names(display_data) && "Year_Week" %in% names(display_data)) {
      analyzed_numeric <- analyzed
      analyzed_numeric$Correlation <- as.numeric(analyzed_numeric$Correlation)

      # Calculate stats across all experiments for each cell line
      all_stats <- analyzed_numeric %>%
        group_by(Cell_line) %>%
        summarise(
          Avg_Corr_All = mean(Correlation, na.rm = TRUE),
          Std_Corr_All = sd(Correlation, na.rm = TRUE),
          .groups = "drop"
        )

      # Calculate stats for same week for each cell line + week combination
      week_stats <- analyzed_numeric %>%
        group_by(Cell_line, Year_Week) %>%
        summarise(
          Avg_Corr_Week = mean(Correlation, na.rm = TRUE),
          Std_Corr_Week = sd(Correlation, na.rm = TRUE),
          .groups = "drop"
        )

      # Merge stats back to display_data
      display_data <- display_data %>%
        left_join(all_stats, by = "Cell_line") %>%
        left_join(week_stats, by = c("Cell_line", "Year_Week"))
    }

    # Ensure columns are proper types for filtering
    display_data$Correlation <- as.numeric(display_data$Correlation)
    display_data$Gate_ID <- as.character(display_data$Gate_ID)
    display_data$Experiment <- as.character(display_data$Experiment)
    if("Notes" %in% names(display_data)) {
      display_data$Notes <- as.character(display_data$Notes)
    }

    # Ensure numeric columns are numeric (not character)
    if("HA_Pos_Pct" %in% names(display_data)) {
      display_data$HA_Pos_Pct <- as.numeric(display_data$HA_Pos_Pct)
      # Convert from 0-100 range to 0-1 range for formatPercentage
      display_data$HA_Pos_Pct <- display_data$HA_Pos_Pct / 100
    }
    if("Strength_Ratio" %in% names(display_data)) {
      display_data$Strength_Ratio <- as.numeric(display_data$Strength_Ratio)
    }

    # Create datatable with numeric columns (formatting applied separately)
    dt <- datatable(display_data,
              selection = 'multiple',
              options = list(
                pageLength = 10,
                scrollX = TRUE,
                search = list(regex = FALSE, caseInsensitive = TRUE)
              ),
              filter = list(position = 'top', clear = FALSE)) %>%
      formatRound('Correlation', digits = 4)

    # Add formatting for optional columns
    if("HA_Pos_Pct" %in% names(display_data)) {
      dt <- dt %>% formatPercentage('HA_Pos_Pct', digits = 2)
    }
    if("Strength_Ratio" %in% names(display_data)) {
      dt <- dt %>% formatRound('Strength_Ratio', digits = 4)
    }
    # Format statistics columns
    if("Avg_Corr_All" %in% names(display_data)) {
      dt <- dt %>% formatRound('Avg_Corr_All', digits = 4)
    }
    if("Std_Corr_All" %in% names(display_data)) {
      dt <- dt %>% formatRound('Std_Corr_All', digits = 4)
    }
    if("Avg_Corr_Week" %in% names(display_data)) {
      dt <- dt %>% formatRound('Avg_Corr_Week', digits = 4)
    }
    if("Std_Corr_Week" %in% names(display_data)) {
      dt <- dt %>% formatRound('Std_Corr_Week', digits = 4)
    }

    dt
  })
  
  # Store previous selection to detect additions vs removals
  rv$prev_selected_rows <- reactiveVal(integer(0))
  
  # Auto-select samples from same cell line (but allow deselection)
  observeEvent(input$comparison_sample_selector_rows_selected, {
    req(rv$all_results)

    analyzed <- rv$all_results[!is.na(rv$all_results$Correlation), ]
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

    analyzed <- rv$all_results[!is.na(rv$all_results$Correlation), ]
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

    analyzed <- rv$all_results[!is.na(rv$all_results$Correlation), ]
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
                           max(c(0, plot_data$Correlation, plot_summary$mean_corr), na.rm = TRUE) + 0.3),
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
          select(Experiment, Correlation, HA_Pos_Pct, Strength_Ratio, N_cells)

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

      # Create separate dataframes for HA_Pos_Pct and Strength_Ratio
      ha_pos_df <- data.frame(Experiment = unique_experiments, stringsAsFactors = FALSE)
      strength_df <- data.frame(Experiment = unique_experiments, stringsAsFactors = FALSE)

      for(i in 1:nrow(unique_groups)) {
        cell_line <- unique_groups$Cell_line[i]
        mutation <- unique_groups$Mutation[i]
        col_name <- col_names[i]

        cell_data <- plot_data %>%
          filter(Cell_line == cell_line) %>%
          select(Experiment, HA_Pos_Pct, Strength_Ratio)

        # HA_Pos_Pct values
        ha_values <- sapply(unique_experiments, function(exp) {
          match_idx <- which(cell_data$Experiment == exp)
          if(length(match_idx) > 0) {
            ha_pct <- as.numeric(cell_data$HA_Pos_Pct[match_idx[1]])
            if(!is.na(ha_pct)) {
              return(round(ha_pct, 2))
            }
          }
          return(NA_real_)
        })

        # Strength_Ratio values
        strength_values <- sapply(unique_experiments, function(exp) {
          match_idx <- which(cell_data$Experiment == exp)
          if(length(match_idx) > 0) {
            strength <- as.numeric(cell_data$Strength_Ratio[match_idx[1]])
            if(!is.na(strength)) {
              return(round(strength, 4))
            }
          }
          return(NA_real_)
        })

        ha_pos_df[[col_name]] <- ha_values
        strength_df[[col_name]] <- strength_values
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
      
      # Add warnings to other sheets
      ha_pos_df$Low_Cell_Warning <- warnings_col
      strength_df$Low_Cell_Warning <- warnings_col

      notes_df <- data.frame(
        Note = c("Three sheets: Correlation, HA_Pos_Pct, and Strength_Ratio",
                 "Correlation values: 4 decimal places",
                 "HA_Pos_Pct values: percentage (2 decimal places)",
                 "Strength_Ratio values: 4 decimal places",
                 "Each row = one experiment, each column = one cell line",
                 "Low_Cell_Warning lists cell lines with N_cells < 500",
                 paste0("Data exported: ", Sys.time())),
        stringsAsFactors = FALSE
      )

      write_xlsx(list(
        "Correlation" = output_df,
        "HA_Pos_Pct" = ha_pos_df,
        "Strength_Ratio" = strength_df,
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