############################################################
## 00_plot_style_and_palettes.r
## Shared plotting colors and themes
############################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(paletteer)
})


############################################################
## 1. SAMPLE TYPES
############################################################

sample_cols <- c(
  Tumor  = "#BF816B",
  Plasma = "#467378",
  Buffy  = "#C29961"
)


############################################################
## 2. SOLID CLINICAL ANNOTATIONS
############################################################

grade_cols <- c(
  `2` = "#DFA398",
  `3` = "#68855C",
  `4` = "#526A83"
)


############################################################
## 3. DMR DIRECTION COLORS
############################################################

solid_direction_cols <- c(
  Tumor_higher  = sample_cols["Tumor"],
  Plasma_higher = sample_cols["Plasma"]
)

octane_direction_cols <- c(
  Buffy_higher  = sample_cols["Buffy"],
  Plasma_higher = sample_cols["Plasma"]
)


############################################################
## 4. CORRELATION HEATMAP PALETTE
##
## Lower correlation = muted blue
## Mid correlation   = off-white
## Higher correlation = muted clay/red
############################################################

correlation_palette <- grDevices::colorRampPalette(
  c(
    "#526A83",
    "#F4F1EA",
    "#BF816B"
  )
)(101)


############################################################
## 5. DELTA-BETA / Z-SCORE HEATMAP PALETTE
##
## Negative = blue
## Zero     = off-white
## Positive = clay/red
############################################################

diverging_methylation_palette <- grDevices::colorRampPalette(
  c(
    "#526A83",
    "#AADCE0",
    "#F4F1EA",
    "#ECC9A0",
    "#BF816B"
  )
)(101)


############################################################
## 6. NEUTRAL COLORS
############################################################

neutral_cols <- c(
  background = "#D9D9D9",
  light      = "#F4F1EA",
  medium     = "#A5A596",
  dark       = "#5A4B3C"
)


############################################################
## 7. INTEGRATION COLORS
############################################################

integration_cols <- c(
  Background = "#D9D9D9",
  P1 = "#68855C",
  P2 = "#526A83",
  P3 = "#C29961"
)


############################################################
## 8. SHARED GGPLOT THEME
############################################################

theme_project <- function(base_size = 13) {

  theme_minimal(
    base_size = base_size
  ) +

    theme(

      plot.title = element_text(
        face = "plain",
        size = base_size + 4
      ),

      plot.subtitle = element_text(
        face = "plain",
        size = base_size + 1
      ),

      axis.title = element_text(
        face = "plain"
      ),

      legend.title = element_text(
        face = "plain"
      ),

      strip.text = element_text(
        face = "plain"
      ),

      panel.grid.minor = element_blank(),

      legend.position = "right"
    )
}