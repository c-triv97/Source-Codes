if (!require("pacman")) install.packages("pacman")
cran.packages <- c("tidyverse",
                   "ggplot2",
                   "ggrepel",
                   "RColorBrewer",
                   "stringr",
                   "viridis",
                   "data.table",
                   "tidyfast",
                   "rstatix",
                   "ggpubr",
                   "ggprism"
)

pacman::p_load(cran.packages, character.only = TRUE)
              
#lapply(cran.packages, packages.load, character.only = TRUE)

theme_scientific <-  theme_light(base_size = 18) + theme(
    panel.background = element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA),
    prism.ticks.length.y = unit(4, "pt"), 
    strip.background = element_rect(fill = "transparent"), 
    strip.text = element_text(color = "black", size = 18, face = "bold.italic"), 
    panel.spacing = unit(2, "lines")
    )

theme_paper <-  theme_light(base_size = 18) + theme(
    panel.background = element_rect(fill='transparent'),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    plot.background = element_rect(fill='transparent', color=NA),
    prism.ticks.length.y = unit(4, "pt"), 
    strip.background = element_rect(fill = "transparent"), 
    strip.text = element_text(color = "black", size = 18, face = "bold.italic"),
    panel.border = element_blank(),
    axis.line = element_line(colour = "grey"),
    panel.spacing = unit(2, "lines"))
