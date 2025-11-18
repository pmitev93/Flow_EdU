# ==============================================================================
# GATE DEFINITIONS - quadrant (Quadrant Gating Strategy)
# ==============================================================================

GATES <- list()

# Gate 1: Debris removal (FSC-A vs SSC-A) - from gid17
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

# Gate 2: Singlets (FSC-A vs FSC-H) - from gid17
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

# Gate 3: Live cells (DCM-A vs SSC-A) - from gid17
GATES$live_cells <- matrix(c(
  122, 181000,
  18000, 181000,
  18000, 13000000,
  122, 13000000,
  122, 181000
), ncol = 2, byrow = TRUE)
colnames(GATES$live_cells) <- c("DCM-A", "SSC-A")

# Gate 4: S-phase outliers (FxCycle-A vs EdU-A) - from gid17
GATES$s_phase_outliers <- matrix(c(
  1.5e6, 100,
  9e6, 100,
  9e6, 5e6,
  1.5e6, 5e6,
  1.5e6, 100
), ncol = 2, byrow = TRUE)
colnames(GATES$s_phase_outliers) <- c("FxCycle-A", "EdU-A")

# Gate 5: FxCycle quantile - 5% to 95%
GATES$fxcycle_quantile <- list(
  type = "quantile_range",
  parameter = "FL6-A",
  probs = c(0.05, 0.95),
  description = "Remove FxCycle outliers (keep 5th-95th percentile)"
)

# Gate 6: Top 45% EdU + FxCycle 50-85%
GATES$edu_fxcycle_sphase <- list(
  type = "dual_quantile",
  edu_parameter = "FL1-A",
  edu_prob = 0.55,  # Top 45% means threshold at 55th percentile
  fxcycle_parameter = "FL6-A",
  fxcycle_probs = c(0.50, 0.85),
  description = "Select top 45% EdU AND FxCycle 50th-85th percentile"
)

# Gate 7: Quadrant gating (HA vs EdU)
GATES$quadrant <- list(
  type = "paired_control_quadrant",
  ha_parameter = "FL3-A",
  edu_parameter = "FL1-A",
  control_suffix = "Dox-",
  test_suffix = "Dox+",
  ha_prob = 0.98,
  edu_prob = 0.02,
  description = "Quadrant gating using paired Dox- control for both HA and EdU thresholds"
)

# Gate strategy metadata
GATE_STRATEGY <- list(
  id = "quadrant",
  name = "Quadrant Gating Strategy",
  description = "Quadrant gating with paired Dox- controls, calculates Strength_Ratio",
  created = Sys.time(),
  analysis_type = "quadrant_ratio"  # Flag to indicate this uses ratio instead of correlation
)
