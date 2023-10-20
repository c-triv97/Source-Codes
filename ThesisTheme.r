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

packages.load <- function(x) {
  if(require(x)){
    print(paste(x, " is loaded correctly"))
} else {
    print(paste(x, " is trying to install"))
    install.packages(x)
    if(require(x)){
        print(paste(x, " is installed and loaded"))
    } else {
        stop("could not install")
    }
}
              
lapply(cran.packages, packages.load, character.only = TRUE)

theme_thesis <-  theme_light(base_size = 18) + theme(
    panel.background = element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA),
    prism.ticks.length.y = unit(4, "pt")
    )
