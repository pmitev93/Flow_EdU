# ==============================================================================
# GATE DEFINITIONS - Simplified (Gates 3 & 4 use vertical thresholds)
# ==============================================================================

GATES <- list()

# Gate 1: Debris removal (FSC-A vs SSC-A)
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
  1.4e6, 4.1e6
), ncol = 2, byrow = TRUE)
colnames(GATES$debris) <- c("FSC-A", "SSC-A")

# Gate 2: Singlets (FSC-A vs FSC-H)
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
  5.8e6, 2.5e6
), ncol = 2, byrow = TRUE)
colnames(GATES$singlet) <- c("FSC-A", "FSC-H")

# Gate 3: Live cells (DCM-A threshold - keep everything LEFT of the line)
GATES$live_cells <- list(
  type = "vertical_threshold",
  parameter = "DCM-A",
  threshold = 18000,  # Keep cells with DCM-A < 18000
  direction = "less_than",
  description = "Live cells (DCM-A below threshold)"
)

# Gate 4: S-phase selection (FxCycle-A range - keep everything BETWEEN two lines)
GATES$s_phase_outliers <- list(
  type = "vertical_range",
  parameter = "FxCycle-A",
  lower_threshold = 1.5e6,  # Left vertical line
  upper_threshold = 9e6,     # Right vertical line
  description = "S-phase cells (FxCycle-A between thresholds)"
)

# Gate 5: FxCycle quantile
GATES$fxcycle_quantile <- list(
  type = "quantile_range",
  parameter = "FL6-A",
  probs = c(0.01, 0.90),
  description = "Remove FxCycle outliers (keep 1st-90th percentile)"
)

# Gate 6: Top 50% EdU + FxCycle range
GATES$edu_fxcycle_sphase <- list(
  type = "dual_quantile",
  edu_parameter = "FL1-A",
  edu_prob = 0.50,
  fxcycle_parameter = "FL6-A",
  fxcycle_probs = c(0.01, 0.90),
  description = "Select top 50% EdU AND FxCycle 1st-90th percentile"
)

# Gate 7: HA positive
GATES$ha_positive <- list(
  type = "control_based_threshold",
  parameter = "FL3-A",
  control_pattern = "Empty_Vector",
  prob = 0.98,
  description = "HA-positive cells (above 98th percentile of empty vector control)"
)

# Gate strategy metadata
GATE_STRATEGY <- list(
  id = "simplified",
  name = "Simplified Gate Strategy",
  description = "Gates 3 & 4 use vertical thresholds instead of polygons",
  created = Sys.time()
)
