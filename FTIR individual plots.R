# =============================================================================
# FTIR Analysis - EPS from Anaerobic Digester Sludge
# Aalborg West (AW) - 5 g/L
# Processing: Savitzky-Golay smoothing, rubberband baseline correction,
#             zero-one normalisation (ir package, Teickner 2025)
# =============================================================================

# --- 1. Load packages ---------------------------------------------------------
library(ir)
library(tidyverse)

# --- 2. Define file paths -----------------------------------------------------
data_dir <- r"(C:\Users\huong\OneDrive - Aalborg Universitet\8. semester\Project\Lab\FTIR)"

files_aw5 <- c(
  "AW 5 gL Alkaline 1.1.dpt",
  "AW 5 gL Alkaline 1.2.dpt",
  "AW 5 gL Alkaline 1.3.dpt",
  "AW 5 gL CER 1.1.dpt",
  "AW 5 gL CER 1.2.dpt",
  "AW 5 gL CER 2.1.dpt"
)

# --- 3. Colour palettes -------------------------------------------------------
# AW concentrations: sequential blues (all visible on white)
aw_colours <- c(
  "5"  = "#92C5DE",
  "8"  = "#4393C3",
  "10" = "#2166AC",
  "14" = "#1A6099",
  "18" = "#08519C",
  "22" = "#08306B"
)

# Køng samples: distinct colours, visible on white
kong_colours <- c(
  "Køng 1" = "#E69F00",
  "Køng 2" = "#009E73",
  "Køng 3" = "#CC79A7",
  "Køng 4" = "#D55E00",
  "Køng 5" = "#000000"
)

# Extraction method colours for AW 5 g/L plots
method_colours <- c("Alkaline" = "#2166AC", "CER" = "#2166AC")

# --- 4. Helper: parse filename into metadata ----------------------------------
parse_filename <- function(fname) {
  base   <- tools::file_path_sans_ext(fname)
  pellet <- as.integer(str_extract(base, "(?<=\\s)(\\d+)(?=\\.\\d+$)"))
  side   <- as.integer(str_extract(base, "(?<=\\.)(\\d+)$"))
  method <- case_when(
    str_detect(base, "Alkaline") ~ "Alkaline",
    str_detect(base, "CER")      ~ "CER",
    TRUE                         ~ NA_character_
  )
  conc <- str_extract(base, "\\d+(?= gL)")
  tibble(
    filename  = fname,
    method    = method,
    conc_gL   = conc,
    pellet    = pellet,
    side      = side,
    sample_id = paste0("P", pellet, "S", side)
  )
}

# --- 5. Load spectra ----------------------------------------------------------
# ir package requires spectra as data frames with columns x (wavenumber) and y (absorbance)
load_dpt <- function(fname, dir) {
  path <- file.path(dir, fname)
  df   <- read_csv(path, col_names = c("x", "y"), show_col_types = FALSE)
  # Sort ascending so ir processing works correctly
  # Trim to biologically relevant window (adjust lower limit to 400 if clean below 600)
  df |>
    filter(x >= 600 & x <= 4000) |>
    arrange(x)
}

# Build metadata table
meta <- map_dfr(files_aw5, parse_filename)

# Load spectra as list of data frames (x = wavenumber, y = absorbance)
spectra_list <- map(files_aw5, ~ load_dpt(.x, data_dir))

# Construct ir object
ir_obj <- ir_new_ir(
  spectra  = spectra_list,
  metadata = meta
)

# --- 6. Pre-processing pipeline -----------------------------------------------
# 6a. Savitzky-Golay smoothing
ir_smooth_obj <- ir_smooth(ir_obj, method = "sg")

# 6b. Rubberband baseline correction
ir_bc_obj <- ir_bc(ir_smooth_obj, method = "rubberband")

# 6c. Zero-one normalisation
ir_norm_obj <- ir_normalise(ir_bc_obj, method = "zeroone")

# --- 7. Extract to long-format data frame ------------------------------------
ir_to_long <- function(ir_obj) {
  ir_obj |>
    as_tibble() |>
    mutate(spectra = map(spectra, as_tibble)) |>
    unnest(spectra) |>
    rename(wavenumber = x, absorbance = y)
}

df_processed <- ir_to_long(ir_norm_obj)

# Add facet label and sort ascending for correct plotting with scale_x_reverse
df_processed <- df_processed |>
  mutate(method_label = case_when(
    method == "Alkaline" ~ "Alkaline extraction",
    method == "CER"      ~ "CER extraction"
  )) |>
  arrange(method, sample_id, wavenumber)

# --- 8. Compute mean +/- SD ---------------------------------------------------
df_mean <- df_processed |>
  group_by(method, method_label, wavenumber) |>
  summarise(
    mean_abs = mean(absorbance, na.rm = TRUE),
    sd_abs   = sd(absorbance,   na.rm = TRUE),
    .groups  = "drop"
  ) |>
  arrange(method, wavenumber)

# --- 9. Plot: mean +/- SD ribbon, faceted by extraction method ---------------
p_combined <- ggplot() +
  # SD ribbon
  geom_ribbon(
    data = df_mean,
    aes(x    = wavenumber,
        ymin  = mean_abs - sd_abs,
        ymax  = mean_abs + sd_abs,
        fill  = method),
    alpha  = 0.25,
    colour = NA
  ) +
  # Mean line
  geom_line(
    data = df_mean,
    aes(x      = wavenumber,
        y      = mean_abs,
        colour = method),
    linewidth = 0.7
  ) +
  scale_colour_manual(values = method_colours, name = "Extraction method") +
  scale_fill_manual(values   = method_colours, name = "Extraction method") +
  # Reverse x-axis: data is ascending, axis display is high -> low (FTIR convention)
  scale_x_reverse(
    breaks = seq(4000, 600, by = -500),
    limits = c(4000, 600),
    expand = c(0, 0)
  ) +
  facet_wrap(~ method_label, ncol = 1) +
  labs(
    title    = "ATR-FTIR spectra \u2014 Aalborg West EPS, 5 g/L",
    subtitle = "Mean \u00b1 SD ribbon across all replicates",
    x        = expression("Wavenumber (cm"^{-1}*")"),
    y        = "Normalised absorbance (a.u.)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "grey92"),
    strip.text       = element_text(face = "bold"),
    legend.position  = "none",   # method shown in facet label
    panel.grid.minor = element_blank()
  )

print(p_combined)
