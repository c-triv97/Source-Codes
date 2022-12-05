# functions for performing calculations from echo measurements in imageJ 
library(tidyverse)

calc_hr <- function(dat, msec = 1286.67){  #dat is table containing ID and number of peaks per image
  dat %>%
    group_by(animal_id)%>%
    dplyr::summarise(mean_peaks = mean(peaks)) %>%
    mutate(bpm = 60/(msec/mean_peaks)*1000)
}

average_echo <- function(data, measures = c(AWTd, EDD, PWTd, AWTs, ESD, PWTs)){ # animal id in one column, AWTd, EDD, PWTd, AWTs, ESD and PWTs in columns 
    data %>%
        group_by(animal_id) %>%
        summarise(across(measures, ~ mean(.x, na.rm = TRUE)))
}

calculations <- function(data){ # ASE cube = 1.04*((IVSTd+LVIDd+PWTd)3–LVIDd3) Devereux correction = 0.8(ASEcube)+0.6
  data %>%
    mutate(ASE = 1.04 * ((((awtd + edd + pwtd) * 10)^3) - ((edd * 10)^3))) %>%
    mutate(LVM = (0.8 * (ASE) + 0.6) / 1000) %>%
    mutate(FS = (edd - esd) / edd * 100) %>%
    mutate(RWT = ((awtd + pwtd) * 10) / (edd * 10)) %>%
    mutate(EDDV = (edd^3) * 1.047,
           ESDV = (esd^3) * 1.047)
}

functionalcalcs <- function(data){
  data %>%
    mutate(SV = (EDDV - ESDV)) %>%
    mutate(CO = (SV * bpm)) %>%
    mutate(EF = (SV / EDDV) * 100)
}

corrections <- function(x) {
  x / 1000
}

remove_whitespace <- function(x) {
  gsub("\\s+", "", x) 
}
