## RQ value:
#generate delta Ct: target - housekeeper
#generate mean delta Ct per sample
#generate mean delta ct per group/condition
#set `control` value as reference group or condition mean
#generate ddCt per sample: dCt - control
#log ddCt
#calc. RQ error bars `2^-(dCt -/+ se.dCt)`

if (!require("pacman")) install.packages("pacman")

cran.packages <- c("tidyverse",
                   "rstatix"
)

pacman::p_load(cran.packages, character.only = TRUE)
source("https://raw.githubusercontent.com/c-triv97/Source-Codes/main/GDPR.R")
source("https://raw.githubusercontent.com/c-triv97/Source-Codes/main/GraphThemes.r")

dummy_dat <- data.frame(
    sample = c(1:6),
    group = c(
        rep(
            "control",
            3), 
        rep(
            "exp",
            3)
        ),
    housekeeping.ct = c(
        runif(
            3, 
            min=18,
            max=22
            ),
        runif(
            3, 
            min=18,
            max=22
            )
        ),
    gene.of.interest.ct = c(
        runif(
            3, 
            min=22,
            max=26
            ),
        runif(
            3, 
            min=26,
            max=28
            )
        )
)       

#### functions for qPCR
# data.table needs to have columns for:
# grouping variable, ct of housekeeper and ct for gene of interst - can be multiple genes of interest but functions needs to be performed separately in this case
# look at "dummy_dat" for minimum reproducible data

calc.dCt <- function(data,
                     housekeeper,
                     gene.of.interest){

    dat1 <- data %>% 
            mutate_if(is.character, as.factor) %>% 
            mutate(dCt = {{gene.of.interest}}-{{housekeeper}})
    
    return(dat1)
}

calc.RQ <- function(data,
                    grouping = group, 
                    control = "control"){

    cont.dCt <- subset(
        data, {{group}} == control
    ) %>%
    pull("dCt")

    dat1 <- data %>% 
            group_by({{grouping}}) %>%
            summarise(mean.dCt = mean(dCt), 
                      sd.dCt = sd(dCt),
                      sem.dCt = standard_error(dCt)) %>%
            mutate(ddCt = mean.dCt - mean(cont.dCt)) %>%
            mutate(RQ = 2^-(ddCt), 
                   RQ.max = 2^-(ddCt - sem.dCt),
                   RQ.min = 2^-(ddCt + sem.dCt))
    
    return(dat1)
            
}

dummy_dat1 <- calc.dCt(data = dummy_dat, 
                       housekeeper = housekeeping.ct, 
                       gene.of.interest = gene.of.interest.ct)

dummy_dat2 <- calc.RQ(dummy_dat1)

# example plot 
plt <- ggbarplot(dummy_dat2,
          x= "group",
          y= "RQ",
          color="white",
          fill= "grey") +
    geom_errorbar(aes(group = group,
                  ymax = RQ.max,
                  ymin = RQ.min, 
                  color = group), width = 0.1,
                  linewidth = 1) +
    xlab("") +
    guides(color = "none") +
    scale_y_continuous(guide = "prism_minor", 
                       expand = expansion(mult = c(0, .1))) +
    theme_scientific
    
