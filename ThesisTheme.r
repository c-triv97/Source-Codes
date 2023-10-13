require(installr)
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
lapply(cran.packages, require2, character.only = TRUE)

theme_thesis <-  theme_light(base_size = 18) + theme(
    panel.background = element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA),
    prism.ticks.length.y = unit(4, "pt")
    )
