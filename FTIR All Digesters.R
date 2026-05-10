# =============================================================================
# FTIR Analysis - Koeng 1-5 + Aalborg West 10 g/L
# Individual plots (one per sample) + waterfall comparison
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
sample_colours <- c(
  "AW 10"   = "#2171B5",
  "Koeng 1" = "#E69F00",
  "Koeng 2" = "#009E73",
  "Koeng 3" = "#CC79A7",
  "Koeng 4" = "#D55E00",
  "Koeng 5" = "#000000"
)

sample_labels <- c(
  "AW 10"   = "AW 10 g/L",
  "Koeng 1" = "K\u00f8ng 1",
  "Koeng 2" = "K\u00f8ng 2",
  "Koeng 3" = "K\u00f8ng 3",
  "Koeng 4" = "K\u00f8ng 4",
  "Koeng 5" = "K\u00f8ng 5"
)

# Waterfall order: AW on top, Koeng 1-5 going downward
# Bottom to top: Koeng 5, Koeng 4, Koeng 3, Koeng 2, Koeng 1, AW 10
sample_levels <- c("Koeng 5", "Koeng 4", "Koeng 3", "Koeng 2", "Koeng 1", "AW 10")
offset_step   <- 1.0

# --- 4. Define all files ------------------------------------------------------
all_files <- tribble(
  ~sample,   ~method,    ~filename,
  # AW 10 g/L
  "AW 10", "Alkaline", "AW 10 gL Alkaline 1.1.dpt",
  "AW 10", "Alkaline", "AW 10 gL Alkaline 1.2.dpt",
  "AW 10", "Alkaline", "AW 10 gL Alkaline 2.1.dpt",
  "AW 10", "Alkaline", "AW 10 gL Alkaline 2.2.dpt",
  "AW 10", "CER",      "AW 10 gL CER 1.1.dpt",
  "AW 10", "CER",      "AW 10 gL CER 1.2.dpt",
  "AW 10", "CER",      "AW 10 gL CER 2.1.dpt",
  "AW 10", "CER",      "AW 10 gL CER 2.2.dpt",
  # Koeng 1
  "Koeng 1", "Alkaline", "K\u00f8ng 1 Alkaline 1.1.dpt",
  "Koeng 1", "Alkaline", "K\u00f8ng 1 Alkaline 1.2.dpt",
  "Koeng 1", "Alkaline", "K\u00f8ng 1 Alkaline 2.1.dpt",
  "Koeng 1", "Alkaline", "K\u00f8ng 1 Alkaline 2.2.dpt",
  "Koeng 1", "CER",      "K\u00f8ng 1 CER 1.1.dpt",
  "Koeng 1", "CER",      "K\u00f8ng 1 CER 1.2.dpt",
  "Koeng 1", "CER",      "K\u00f8ng 1 CER 2.1.dpt",
  "Koeng 1", "CER",      "K\u00f8ng 1 CER 2.2.dpt",
  # Koeng 2
  "Koeng 2", "Alkaline", "K\u00f8ng 2 Alkaline 1.1.dpt",
  "Koeng 2", "Alkaline", "K\u00f8ng 2 Alkaline 1.2.dpt",
  "Koeng 2", "Alkaline", "K\u00f8ng 2 Alkaline 2.1.dpt",
  "Koeng 2", "Alkaline", "K\u00f8ng 2 Alkaline 2.2.dpt",
  "Koeng 2", "CER",      "K\u00f8ng 2 CER 1.1.dpt",
  "Koeng 2", "CER",      "K\u00f8ng 2 CER 1.2.dpt",
  "Koeng 2", "CER",      "K\u00f8ng 2 CER 2.1.dpt",
  "Koeng 2", "CER",      "K\u00f8ng 2 CER 2.2.dpt",
  # Koeng 3
  "Koeng 3", "Alkaline", "K\u00f8ng 3 Alkaline 1.1.dpt",
  "Koeng 3", "Alkaline", "K\u00f8ng 3 Alkaline 1.2.dpt",
  "Koeng 3", "Alkaline", "K\u00f8ng 3 Alkaline 2.1.dpt",
  "Koeng 3", "Alkaline", "K\u00f8ng 3 Alkaline 2.2.dpt",
  "Koeng 3", "CER",      "K\u00f8ng 3 CER 1.1.dpt",
  "Koeng 3", "CER",      "K\u00f8ng 3 CER 1.2.dpt",
  "Koeng 3", "CER",      "K\u00f8ng 3 CER 2.1.dpt",
  "Koeng 3", "CER",      "K\u00f8ng 3 CER 2.2.dpt",
  # Koeng 4
  "Koeng 4", "Alkaline", "K\u00f8ng 4 Alkaline 1.1.dpt",
  "Koeng 4", "Alkaline", "K\u00f8ng 4 Alkaline 1.2.dpt",
  "Koeng 4", "Alkaline", "K\u00f8ng 4 Alkaline 1.3.dpt",
  "Koeng 4", "Alkaline", "K\u00f8ng 4 Alkaline 1.4.dpt",
  "Koeng 4", "CER",      "K\u00f8ng 4 CER 1.1.dpt",
  "Koeng 4", "CER",      "K\u00f8ng 4 CER 1.2.dpt",
  "Koeng 4", "CER",      "K\u00f8ng 4 CER 2.1.dpt",
  "Koeng 4", "CER",      "K\u00f8ng 4 CER 2.2.dpt",
  # Koeng 5
  "Koeng 5", "Alkaline", "K\u00f8ng 5 Alkaline 1.1.dpt",
  "Koeng 5", "Alkaline", "K\u00f8ng 5 Alkaline 1.2.dpt",
  "Koeng 5", "Alkaline", "K\u00f8ng 5 Alkaline 2.1.dpt",
  "Koeng 5", "Alkaline", "K\u00f8ng 5 Alkaline 2.2.dpt",
  "Koeng 5", "CER",      "K\u00f8ng 5 CER 1.1.dpt",
  "Koeng 5", "CER",      "K\u00f8ng 5 CER 1.2.dpt",
  "Koeng 5", "CER",      "K\u00f8ng 5 CER 2.1.dpt",
  "Koeng 5", "CER",      "K\u00f8ng 5 CER 2.2.dpt"
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
all_files <- all_files |>
  filter(file.exists(file.path(data_dir, filename))) |>
  mutate(
    pellet    = as.integer(str_extract(
      tools::file_path_sans_ext(filename), "(?<=\\s)(\\d+)(?=\\.\\d+$)")),
    side      = as.integer(str_extract(
      tools::file_path_sans_ext(filename), "(?<=\\.)(\\d+)$")),
    sample_id = paste0("P", pellet, "S", side)
  )

spectra_list <- map(all_files$filename, ~ load_dpt(.x, data_dir))

ir_obj <- ir_new_ir(spectra = spectra_list, metadata = all_files)
ir_obj <- ir_smooth(ir_obj,    method = "sg")
ir_obj <- ir_bc(ir_obj,        method = "rubberband")
ir_obj <- ir_normalise(ir_obj, method = "zeroone")

df <- ir_to_long(ir_obj) |>
  arrange(sample, method, pellet, side, wavenumber)

# --- 7. Compute mean +/- SD ---------------------------------------------------
df_mean <- df |>
  group_by(sample, method, wavenumber) |>
  summarise(
    mean_abs = mean(absorbance, na.rm = TRUE),
    sd_abs   = sd(absorbance,   na.rm = TRUE),
    .groups  = "drop"
  ) |>
  arrange(sample, method, wavenumber)

# --- 8. Shared theme and x scale ----------------------------------------------
ftir_theme <- theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "grey92"),
    strip.text       = element_text(face = "bold"),
    legend.position  = "none",
    panel.grid       = element_blank()
  )

ftir_x <- scale_x_reverse(
  breaks = seq(4000, 600, by = -500),
  limits = c(4000, 600),
  expand = c(0, 0)
)

# --- 9. Individual plots per sample -------------------------------------------
for (ss in names(sample_colours)) {

  df_ss <- df_mean |>
    filter(sample == ss) |>
    mutate(
      method_label = factor(
        case_when(
          method == "Alkaline" ~ "Alkaline extraction",
          method == "CER"      ~ "CER extraction"
        ),
        levels = c("Alkaline extraction", "CER extraction")
      )
    )

  colour <- sample_colours[[ss]]
  label  <- sample_labels[[ss]]

  p <- ggplot(df_ss, aes(x = wavenumber, group = method_label)) +
    geom_ribbon(aes(ymin = mean_abs - sd_abs,
                    ymax = mean_abs + sd_abs),
                alpha = 0.25, colour = NA, fill = colour) +
    geom_line(aes(y = mean_abs), linewidth = 0.7, colour = colour) +
    facet_wrap(~ method_label, ncol = 1) +
    ftir_x +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.05))) +
    labs(
      title    = paste0("ATR-FTIR spectra \u2014 ", label, " EPS, 10 g/L"),
      subtitle = "Mean \u00b1 SD ribbon across all replicates",
      x        = expression("Wavenumber (cm"^{-1}*")"),
      y        = "Normalised absorbance (a.u.)"
    ) +
    ftir_theme

  # Clean filename (remove spaces and slashes)
  fname <- paste0("FTIR_", gsub(" ", "_", ss))
  ggsave(file.path(out_dir, paste0(fname, ".pdf")),
         p, width = 18, height = 14, units = "cm")
  ggsave(file.path(out_dir, paste0(fname, ".png")),
         p, width = 18, height = 14, units = "cm", dpi = 300)
  message("Saved: ", fname)
}

# --- 10. Waterfall: all samples -----------------------------------------------
offset_vals        <- (seq_along(sample_levels) - 1) * offset_step
names(offset_vals) <- sample_labels[sample_levels]

df_wf <- df_mean |>
  mutate(
    sample   = factor(sample, levels = sample_levels),
    offset   = (as.integer(sample) - 1) * offset_step,
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
  scale_colour_manual(values = sample_colours, labels = sample_labels,
                      name = "Sample"),
  scale_fill_manual(values   = sample_colours, labels = sample_labels,
                    name = "Sample"),
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
               aes(x = wavenumber, colour = sample, fill = sample, group = sample)) +
  geom_ribbon(aes(ymin = ymin_off, ymax = ymax_off), alpha = 0.25, colour = NA) +
  geom_line(aes(y = mean_off), linewidth = 0.7) +
  facet_wrap(~ method_label, ncol = 1) +
  wf_scales +
  labs(
    title    = "ATR-FTIR spectra \u2014 K\u00f8ng 1\u20135 & AW 10 g/L",
    subtitle = "Mean \u00b1 SD ribbon; spectra offset vertically for clarity",
    x        = expression("Wavenumber (cm"^{-1}*")"),
    y        = "Absorbance (offset for clarity)"
  ) +
  wf_theme

ggsave(file.path(out_dir, "FTIR_Kong_waterfall_combined.pdf"),
       p_wf, width = 18, height = 24, units = "cm")
ggsave(file.path(out_dir, "FTIR_Kong_waterfall_combined.png"),
       p_wf, width = 18, height = 24, units = "cm", dpi = 300)
message("Saved: FTIR_Kong_waterfall_combined")

# Separate waterfall per method
for (mm in c("Alkaline", "CER")) {
  df_wf_m <- df_wf |> filter(method == mm)

  p_wf_m <- ggplot(df_wf_m,
                   aes(x = wavenumber, colour = sample, fill = sample,
                       group = sample)) +
    geom_ribbon(aes(ymin = ymin_off, ymax = ymax_off), alpha = 0.25, colour = NA) +
    geom_line(aes(y = mean_off), linewidth = 0.7) +
    wf_scales +
    labs(
      title    = paste0("ATR-FTIR spectra \u2014 K\u00f8ng 1\u20135 & AW 10 g/L, ",
                        mm, " extraction"),
      subtitle = "Mean \u00b1 SD ribbon; spectra offset vertically for clarity",
      x        = expression("Wavenumber (cm"^{-1}*")"),
      y        = "Absorbance (offset for clarity)"
    ) +
    wf_theme

  ggsave(file.path(out_dir, paste0("FTIR_Kong_waterfall_", mm, ".pdf")),
         p_wf_m, width = 18, height = 14, units = "cm")
  ggsave(file.path(out_dir, paste0("FTIR_Kong_waterfall_", mm, ".png")),
         p_wf_m, width = 18, height = 14, units = "cm", dpi = 300)
  message("Saved: FTIR_Kong_waterfall_", mm)
}

message("\nAll Koeng + AW 10 g/L plots saved to: ", out_dir)
