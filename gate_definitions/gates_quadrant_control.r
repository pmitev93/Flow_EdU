# ==============================================================================
# GATE DEFINITIONS - quadrant_control (Quadrant with Control-Based Thresholds)
# ==============================================================================
# This strategy uses quadrant gating with EdU and HA thresholds calculated
# from the 121_Empty_Vector_Dox- control sample and applied to all samples

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
  1.9e6, 0.972e6,
  3.4e6, 1.5e6,
  5.8e6, 2.5e6,
  6.8e6, 2.3e6,
  4.9e6, 1.6e6,
  3.7e6, 1.1e6,
  2.8e6, 0.748e6,
  1.9e6, 0.486e6,
  1.7e6, 0.491e6,
  1.5e6, 0.508e6,
  1.6e6, 0.660e6,
  1.9e6, 0.972e6
), ncol = 2, byrow = TRUE)
colnames(GATES$singlet) <- c("FSC-A", "FSC-H")

# Gate 3: Live cells (DCM-A vs SSC-A)
GATES$live_cells <- matrix(c(
  122, 181000,
  18000, 181000,
  18000, 13000000,
  122, 13000000,
  122, 181000
), ncol = 2, byrow = TRUE)
colnames(GATES$live_cells) <- c("DCM-A", "SSC-A")

# Gate 4: S-phase outliers (FxCycle-A vs EdU-A)
GATES$s_phase_outliers <- matrix(c(
  1.5e6, 100,
  9e6, 100,
  9e6, 5e6,
  1.5e6, 5e6,
  1.5e6, 100
), ncol = 2, byrow = TRUE)
colnames(GATES$s_phase_outliers) <- c("FxCycle-A", "EdU-A")

# Gate 5: FxCycle quantile
GATES$fxcycle_quantile <- list(
  type = "quantile_range",
  parameter = "FL6-A",
  probs = c(0.01, 0.90),
  description = "Remove FxCycle outliers (keep 1st-90th percentile)"
)

# Gate 6: EdU + FxCycle S-phase (45-55th percentile EdU)
GATES$edu_fxcycle_sphase <- list(
  type = "dual_quantile",
  edu_parameter = "FL1-A",
  edu_prob = 0.45,  # Threshold at 45th percentile (45-55% range)
  fxcycle_parameter = "FL6-A",
  fxcycle_probs = c(0.30, 0.70),
  description = "Select EdU 45-55th percentile AND FxCycle 30th-70th percentile"
)

# Gate 7: Quadrant gating with control-based thresholds
GATES$quadrant <- list(
  type = "control_based_quadrant",
  ha_parameter = "FL3-A",
  edu_parameter = "FL1-A",
  control_pattern = "121_Empty_Vector_Dox-",
  ha_prob = 0.98,
  edu_prob = 0.02,
  description = "Quadrant gating using 121_Empty_Vector_Dox- control for both HA and EdU thresholds"
)

# Gate strategy metadata
GATE_STRATEGY <- list(
  id = "quadrant_control",
  name = "Quadrant with Control-Based Thresholds",
  description = "Quadrant gating strategy with EdU and HA thresholds calculated from 121_Empty_Vector_Dox- control and applied to all samples",
  created = Sys.time(),
  analysis_type = "quadrant_ratio"  # Flag to indicate this uses ratio instead of correlation
)
