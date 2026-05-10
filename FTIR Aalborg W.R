# =============================================================================
# FTIR Analysis - Aalborg West EPS, all concentrations
# Individual plots (one per concentration) + waterfall comparison
# Processing: Savitzky-Golay smoothing, rubberband baseline correction,
#             zero-one normalisation (ir package, Teickner 2025)
# Output: saved to plots\ subfolder
# =============================================================================

# --- 1. Load packages ---------------------------------------------------------
library(ir)
library(tidyverse)

# --- 2. Define file paths -----------------------------------------------------
data_dir <- r"(C:\Users\huong\OneDrive - Aalborg Universitet\8. semester\Project\Lab\FTIR)"
out_dir  <- r"(C:\Users\huong\OneDrive - Aalborg Universitet\8. semester\Project\Lab\FTIR\plots)"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --- 3. Colour palette --------------------------------------------------------
# Sequential blues, lightest = lowest concentration
aw_colours <- c(
  "5"  = "#6BAED6",
  "8"  = "#4393C3",
  "10" = "#2171B5",
  "14" = "#1561A9",
  "18" = "#08519C",
  "22" = "#08306B"
)

offset_step  <- 1
conc_levels  <- c("5", "8", "10", "14", "18", "22")
conc_labels  <- paste0(conc_levels, " g/L")
names(conc_labels) <- conc_levels

# --- 4. Define all AW files ---------------------------------------------------
# Note: file names corrected to match actual files on disk
aw_files <- tribble(
  ~conc, ~method,    ~filename,
  # AW 5 g/L
  "5",  "Alkaline", "AW 5 gL Alkaline 1.1.dpt",
  "5",  "Alkaline", "AW 5 gL Alkaline 1.2.dpt",
  "5",  "Alkaline", "AW 5 gL Alkaline 1.3.dpt",
  "5",  "Alkaline", "AW 5 gL Alkaline 1.4.dpt",
  "5",  "CER",      "AW 5 gL CER 1.1.dpt",
  "5",  "CER",      "AW 5 gL CER 1.2.dpt",
  "5",  "CER",      "AW 5 gL CER 2.1.dpt",
  # AW 8 g/L
  "8",  "Alkaline", "AW 8 gL Alkaline 1.1.dpt",
  "8",  "Alkaline", "AW 8 gL Alkaline 1.2.dpt",
  "8",  "Alkaline", "AW 8 gL Alkaline 2.1.dpt",
  "8",  "Alkaline", "AW 8 gL Alkaline 2.2.dpt",
  "8",  "CER",      "AW 8 gL CER 1.1.dpt",
  "8",  "CER",      "AW 8 gL CER 2.1.dpt",
  "8",  "CER",      "AW 8 gL CER 2.2.dpt",
  "8",  "CER",      "AW 8 gL CER 2.3.dpt",
  # AW 10 g/L
  "10", "Alkaline", "AW 10 gL Alkaline 1.1.dpt",
  "10", "Alkaline", "AW 10 gL Alkaline 1.2.dpt",
  "10", "Alkaline", "AW 10 gL Alkaline 2.1.dpt",
  "10", "Alkaline", "AW 10 gL Alkaline 2.2.dpt",
  "10", "CER",      "AW 10 gL CER 1.1.dpt",
  "10", "CER",      "AW 10 gL CER 1.2.dpt",
  "10", "CER",      "AW 10 gL CER 2.1.dpt",
  "10", "CER",      "AW 10 gL CER 2.2.dpt",
  # AW 14 g/L
  "14", "Alkaline", "AW 14 gL Alkaline 1.1.dpt",
  "14", "Alkaline", "AW 14 gL Alkaline 1.2.dpt",
  "14", "Alkaline", "AW 14 gL Alkaline 2.1.dpt",
  "14", "Alkaline", "AW 14 gL Alkaline 2.2.dpt",
  "14", "CER",      "AW 14 gL CER 1.1.dpt",
  "14", "CER",      "AW 14 gL CER 1.2.dpt",
  "14", "CER",      "AW 14 gL CER 1.3.dpt",
  "14", "CER",      "AW 14 gL CER 2.1.dpt",
  # AW 18 g/L
  "18", "Alkaline", "AW 18 gL Alkaline 1.1.dpt",
  "18", "Alkaline", "AW 18 gL Alkaline 1.2.dpt",
  "18", "Alkaline", "AW 18 gL Alkaline 2.1.dpt",
  "18", "Alkaline", "AW 18 gL Alkaline 2.2.dpt",
  "18", "CER",      "AW 18 gL CER 1.1.dpt",
  "18", "CER",      "AW 18 gL CER 1.2.dpt",
  "18", "CER",      "AW 18 gL CER 2.1.dpt",
  "18", "CER",      "AW 18 gL CER 2.2.dpt",
  # AW 22 g/L
  "22", "Alkaline", "AW 22 gL Alkaline 1.1.dpt",
  "22", "Alkaline", "AW 22 gL Alkaline 1.2.dpt",
  "22", "Alkaline", "AW 22 gL Alkaline 1.3.dpt",
  "22", "Alkaline", "AW 22 gL Alkaline 1.4.dpt",
  "22", "CER",      "AW 22 gL CER 1.1.dpt",
  "22", "CER",      "AW 22 gL CER 1.2.dpt",
  "22", "CER",      "AW 22 gL CER 2.1.dpt",
  "22", "CER",      "AW 22 gL CER 2.2.dpt"
)

# --- 5. Helper functions ------------------------------------------------------
load_dpt <- function(fname, dir) {
  path <- file.path(dir, fname)
  df   <- read_delim(path, col_names = c("x", "y"),
                     show_col_types = FALSE, delim = NULL)
  df |> filter(x >= 600 & x <= 4000) |> arrange(x)
}
ir_to_long <- function(ir_obj) {
  ir_obj |>
    as_tibble() |>
    mutate(spectra = map(spectra, as_tibble)) |>
    unnest(spectra) |>
    rename(wavenumber = x, absorbance = y)
}

# --- 6. Load and process all spectra ------------------------------------------
aw_files <- aw_files |>
  filter(file.exists(file.path(data_dir, filename))) |>
  mutate(
    pellet    = as.integer(str_extract(
      tools::file_path_sans_ext(filename), "(?<=\\s)(\\d+)(?=\\.\\d+$)")),
    side      = as.integer(str_extract(
      tools::file_path_sans_ext(filename), "(?<=\\.)(\\d+)$")),
    sample_id = paste0("P", pellet, "S", side)
  )

spectra_list <- map(aw_files$filename, ~ load_dpt(.x, data_dir))

ir_obj <- ir_new_ir(spectra = spectra_list, metadata = aw_files)
ir_obj <- ir_smooth(ir_obj,    method = "sg")
ir_obj <- ir_bc(ir_obj,        method = "rubberband")
ir_obj <- ir_normalise(ir_obj, method = "zeroone")

df <- ir_to_long(ir_obj) |>
  arrange(conc, method, pellet, side, wavenumber)

# --- 7. Compute mean +/- SD ---------------------------------------------------
df_mean <- df |>
  group_by(conc, method, wavenumber) |>
  summarise(
    mean_abs = mean(absorbance, na.rm = TRUE),
    sd_abs   = sd(absorbance,   na.rm = TRUE),
    .groups  = "drop"
  ) |>
  arrange(conc, method, wavenumber)

# --- 8. Shared theme ----------------------------------------------------------
ftir_theme <- theme_bw(base_size = 11) +
  theme(
    strip.background   = element_rect(fill = "grey92"),
    strip.text         = element_text(face = "bold"),
    legend.position    = "none",
    panel.grid         = element_blank()
  )

ftir_x <- scale_x_reverse(
  breaks = seq(4000, 600, by = -500),
  limits = c(4000, 600),
  expand = c(0, 0)
)

# --- 9. Individual concentration plots ----------------------------------------
# One plot per concentration, faceted by extraction method
for (cc in conc_levels) {
  
  df_cc <- df_mean |>
    filter(conc == cc) |>
    mutate(
      method_label = factor(
        case_when(
          method == "Alkaline" ~ "Alkaline extraction",
          method == "CER"      ~ "CER extraction"
        ),
        levels = c("Alkaline extraction", "CER extraction")
      )
    )
  
  colour <- aw_colours[[cc]]   # double brackets to get the value not a named vector
  
  p <- ggplot(df_cc,
              aes(x = wavenumber, group = method_label)) +
    geom_ribbon(aes(ymin = mean_abs - sd_abs,
                    ymax = mean_abs + sd_abs),
                alpha = 0.25, colour = NA, fill = colour) +
    geom_line(aes(y = mean_abs), linewidth = 0.7, colour = colour) +
    facet_wrap(~ method_label, ncol = 1) +
    ftir_x +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.05))) +
    labs(
      title    = paste0("ATR-FTIR spectra \u2014 Aalborg West EPS, ", cc, " g/L"),
      subtitle = "Mean \u00b1 SD ribbon across all replicates",
      x        = expression("Wavenumber (cm"^{-1}*")"),
      y        = "Normalised absorbance (a.u.)"
    ) +
    ftir_theme
  
  ggsave(file.path(out_dir, paste0("FTIR_AW_", cc, "gL.pdf")),
         p, width = 18, height = 14, units = "cm")
  ggsave(file.path(out_dir, paste0("FTIR_AW_", cc, "gL.png")),
         p, width = 18, height = 14, units = "cm", dpi = 300)
  message("Saved: FTIR_AW_", cc, "gL")
}

# --- 10. Waterfall: all concentrations ----------------------------------------
# Add offsets
offset_vals        <- (seq_along(conc_levels) - 1) * offset_step
names(offset_vals) <- conc_labels

df_wf <- df_mean |>
  mutate(
    conc     = factor(conc, levels = conc_levels),
    offset   = (as.integer(conc) - 1) * offset_step,
    mean_off = mean_abs + offset,
    ymin_off = (mean_abs - sd_abs) + offset,
    ymax_off = (mean_abs + sd_abs) + offset,
    method_label = factor(
      case_when(
        method == "Alkaline" ~ "Alkaline extraction",
        method == "CER"      ~ "CER extraction"
      ),
      levels = c("Alkaline extraction", "CER extraction")
    )
  )

wf_scales <- list(
  scale_colour_manual(values = aw_colours, labels = conc_labels,
                      name = "Concentration"),
  scale_fill_manual(values   = aw_colours, labels = conc_labels,
                    name = "Concentration"),
  scale_x_reverse(
    breaks = seq(4000, 600, by = -500),
    limits = c(4000, 600),
    expand = c(0, 0)
  ),
  scale_y_continuous(
    breaks = offset_vals,
    labels = names(offset_vals),
    expand = expansion(mult = c(0.08, 0.05))
  )
)

wf_theme <- theme_bw(base_size = 11) +
  theme(
    legend.position    = "none",
    panel.grid         = element_blank(),
    strip.background   = element_rect(fill = "grey92"),
    strip.text         = element_text(face = "bold")
  )

# Combined waterfall (both methods, two rows)
p_wf <- ggplot(df_wf,
               aes(x = wavenumber, colour = conc, fill = conc, group = conc)) +
  geom_ribbon(aes(ymin = ymin_off, ymax = ymax_off), alpha = 0.25, colour = NA) +
  geom_line(aes(y = mean_off), linewidth = 0.7) +
  facet_wrap(~ method_label, ncol = 1) +
  wf_scales +
  labs(
    title    = "ATR-FTIR spectra \u2014 Aalborg West EPS, all concentrations",
    subtitle = "Mean \u00b1 SD ribbon; spectra offset vertically for clarity",
    x        = expression("Wavenumber (cm"^{-1}*")"),
    y        = "Absorbance (offset for clarity)"
  ) +
  wf_theme

ggsave(file.path(out_dir, "FTIR_AW_waterfall_combined.pdf"),
       p_wf, width = 18, height = 24, units = "cm")
ggsave(file.path(out_dir, "FTIR_AW_waterfall_combined.png"),
       p_wf, width = 18, height = 24, units = "cm", dpi = 300)
message("Saved: FTIR_AW_waterfall_combined")

# Separate waterfall per method
for (mm in c("Alkaline", "CER")) {
  df_wf_m <- df_wf |> filter(method == mm)

  p_wf_m <- ggplot(df_wf_m,
                   aes(x = wavenumber, colour = conc, fill = conc, group = conc)) +
    geom_ribbon(aes(ymin = ymin_off, ymax = ymax_off), alpha = 0.25, colour = NA) +
    geom_line(aes(y = mean_off), linewidth = 0.7) +
    wf_scales +
    labs(
      title    = paste0("ATR-FTIR spectra \u2014 Aalborg West EPS, ", mm, " extraction"),
      subtitle = "Mean \u00b1 SD ribbon; spectra offset vertically for clarity",
      x        = expression("Wavenumber (cm"^{-1}*")"),
      y        = "Absorbance (offset for clarity)"
    ) +
    wf_theme

  ggsave(file.path(out_dir, paste0("FTIR_AW_waterfall_", mm, ".pdf")),
         p_wf_m, width = 18, height = 14, units = "cm")
  ggsave(file.path(out_dir, paste0("FTIR_AW_waterfall_", mm, ".png")),
         p_wf_m, width = 18, height = 14, units = "cm", dpi = 300)
  message("Saved: FTIR_AW_waterfall_", mm)
}

message("\nAll AW plots saved to: ", out_dir)
