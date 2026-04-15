# Flow_EdU: EdU-HA Correlation Analysis Tool

An R Shiny application for analyzing the relationship between EdU incorporation and HA-tagged protein expression in flow cytometry experiments.

## Overview

This tool analyzes multi-color flow cytometry data to quantify the correlation between EdU incorporation (S-phase marker) and HA-tagged protein expression. It implements a standardized 8-gate sequential gating strategy and calculates correlation metrics across samples and experiments.

**Target users:** Flow cytometry researchers studying cell cycle and protein expression dynamics.

## What It Does

The application:
- Loads and processes FCS files from multiple experiments
- Applies sequential gating (debris → singlets → live → S-phase → quantiles → EdU/HA)
- Calculates EdU vs HA correlation using linear regression on log-transformed data
- Generates publication-quality plots with customizable formatting
- Exports data to Excel format compatible with GraphPad Prism
- Supports multiple gating strategies for comparing different analysis parameters

## Key Features

### Analysis Capabilities
- **Automated gating**: Sequential 8-gate strategy removes debris, doublets, dead cells, and outliers
- **Correlation analysis**: Linear regression (OLS) on log10-transformed EdU and HA signals
- **Batch processing**: Analyze multiple experiments and samples simultaneously
- **Gating strategy comparison**: Test different gate parameters and compare results
- **Statistical testing**: Paired/unpaired t-tests, ANOVA, non-parametric alternatives

### Visualization
- **Sample Overview**: All 8 gates displayed in a 2×4 grid for quality control
- **Cross-Experiment View**: Compare one sample across multiple experiments
- **Multi-Sample Comparison**: Bar plots with statistical comparisons
- **Publication-quality exports**: SVG, PDF, and PNG formats with customizable dimensions
- **Interactive plots**: Adjustable font sizes, manual axis limits, and plot dimensions

### Data Export
- **Excel export**: Multi-sheet workbooks with correlation, slope, and gating strategy metadata
- **Prism-ready format**: Data structured for direct import into GraphPad Prism
- **Vector graphics**: High-resolution SVG/PDF plots for publications

## Analysis Workflow

### Sequential Gating Strategy

1. **Gate 1: Debris Removal** - Exclude low FSC-A/SSC-A events
2. **Gate 2: Singlets** - FSC-A vs FSC-H doublet discrimination
3. **Gate 3: Live Cells** - DCM-negative (live cell gate)
4. **Gate 4: S-phase Outliers** - Remove outliers from EdU vs FxCycle
5. **Gate 5: FxCycle Quantile** - Select middle 90% of DNA content (1%-90% percentile)
6. **Gate 6: EdU + FxCycle Range** - Top 45% EdU+ cells within DNA content range
7. **Gate 7: HA-Positive** - HA expression above threshold (calculated from Dox- control)
8. **Final: Correlation** - EdU vs HA correlation in HA+ population

### Correlation Calculation

The tool calculates:
- **Pearson correlation (r)**: Linear correlation coefficient between log10(EdU) and log10(HA)
- **R²**: Coefficient of determination (fraction of variance explained)
- **Slope**: Linear regression slope (OLS) - how much EdU changes per unit HA
- **Cell count**: Number of cells in the final gated population

**Why log transformation?** Flow cytometry fluorescence spans orders of magnitude. Log transformation:
- Normalizes the dynamic range
- Makes relationships more linear
- Reduces the influence of extreme outliers
- Matches how flow cytometry data is typically visualized

**Why OLS regression?** HA expression is the independent variable (predictor) and EdU incorporation is the dependent variable (biological response). Ordinary least squares regression is the appropriate method for this directional relationship.

## Installation & Setup

### Requirements
- R (version 4.0 or higher)
- RStudio (recommended)

### R Packages
```r
install.packages(c("shiny", "flowCore", "DT", "openxlsx", "sp"))
```

### Directory Structure
```
Flow_EdU/
├── Flow_GUI.r              # Main Shiny application
├── plotting_functions.r    # Plotting functions
├── Experiments/            # FCS data files (organized by experiment)
├── gate_definitions/       # Gate parameter files (gates_*.r)
└── analysis_cache/         # Cached analysis results (auto-generated)
```

## Usage

### 1. Launch the App
```r
shiny::runApp("Flow_GUI.r")
```

### 2. Load Data
- **Data Selection tab**: Set master folder path (e.g., "Experiments/")
- Select experiments and gating strategy
- Click "Analyze Selected Experiments"

### 3. View Results
- **Results Table**: Overview of all samples with correlation metrics
- **Sample Overview**: Quality control - view all gates for individual samples
- **Multi-Sample Comparison**: Compare correlations across samples with statistics
- **Final: Correlation New**: Publication-quality EdU vs HA scatter plots

### 4. Export Data
- **Download for Prism**: Excel file with correlation, slope, and metadata
- **Download Plot**: Export plots as SVG, PDF, or PNG

## Gating Strategies

Gate parameters are defined in `gate_definitions/gates_*.r` files. Each file specifies:
- Gate coordinates (polygons or thresholds)
- Quantile cutoffs
- Channel assignments

**Default strategy**: `gates_gdef.r`

To create a new strategy:
1. Copy an existing `gates_*.r` file
2. Modify gate parameters as needed
3. Save with a unique identifier (e.g., `gates_custom_20250101.r`)
4. Restart the Shiny app
5. Select your strategy in the Data Selection tab

## Tips for Best Results

### Sample Preparation
- Include Dox- controls (for HA threshold calculation)
- Maintain consistent staining protocols across experiments
- Acquire sufficient events (>10,000 recommended)

### Data Organization
- Organize FCS files by experiment in subdirectories
- Use consistent naming conventions
- Include experimental metadata in filenames or sample names

### Quality Control
- Always check Sample Overview plots to verify gating
- Look for unusual distributions or low cell counts
- Compare results across different gating strategies if uncertain

### Publication Plots
- Use "Final: Correlation New" tab for publication figures
- Export as PDF or SVG for vector graphics
- Adjust dimensions to match journal requirements (typically 600-1200 px)
- Use manual axis limits for consistency across samples

## Troubleshooting

**"No control sample found" error**
- Ensure you have a Dox- Empty Vector control sample
- Check sample naming includes "Dox-" and "Empty_Vector"

**Gates not appearing in dropdown**
- Check that corresponding `gates_*.r` file exists in `gate_definitions/`
- Restart the Shiny app after adding new gate files

**Large SVG file sizes**
- SVG files can be 20-30 MB for multi-panel plots with many cells
- Use PDF format instead - much smaller and equally publication-ready

**Plots look different in exported files**
- This is normal due to graphics device differences
- Exported files are publication-quality and should be used for final figures

## Citation

If you use this tool in your research, please cite: "Pending"

## Technical Details

**Language**: R  
**Framework**: Shiny  
**License**: [To be added]  
**Author**: Petar Mitev

**Contact**: pmitev93@gmail.com

## Version History

- **v1.0** (2025): Initial release with core functionality
  - Sequential gating strategy
  - EdU-HA correlation analysis
  - Publication-quality plots
  - Excel export for Prism

## Acknowledgments

Analysis pipeline developed for studying EdU incorporation dynamics in HA-tagged protein expression systems.
