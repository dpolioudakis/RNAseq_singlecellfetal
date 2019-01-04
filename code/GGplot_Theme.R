# Damon Polioudakis
# 2018-09-09
# Set ggplot2 theme

theme_set(theme_bw())
theme_set(theme_get() + theme(text = element_text(size = 14)))
theme_update(plot.title = element_text(size = 14))
theme_update(
  axis.line = element_line(colour = "black")
  , axis.text = element_text(colour = "black")
  , plot.background = element_blank()
  , panel.border = element_blank()
)

# Publication theme
ggplot_set_theme_publication <-
  theme_bw() +
  theme(
    , text = element_text(size = 10, colour = "black")
    , axis.text = element_text(colour = "black", size = 10)
    , strip.text = element_text(size = 10)
    , plot.title = element_text(size = 10)
    , axis.line = element_line(colour = "black")
    , panel.border = element_blank()
    , panel.grid.major = element_blank()
    , panel.grid.minor = element_blank()
  )
ggplot_set_theme_publication_nolabels <-
  theme_bw() +
  theme(
    , text = element_text(size = 10, colour = "black")
    , strip.text = element_text(size = 10)
    , plot.title = element_text(size = 10)
    , panel.border = element_blank()
    , panel.grid.major = element_blank()
    , panel.grid.minor = element_blank()
    , axis.title = element_blank()
    , axis.text = element_blank()
    , axis.ticks = element_blank()
  )
