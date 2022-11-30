# functions for performing calculations from echo measurements in imageJ 

calc_hr <- function(dat, animal_id, msec = 1286.67){  #dat is table containing ID and number of peaks per image
  dat %>%
    group_by(animal_id)%>%
    summarise(mean_peaks = mean(peaks)) %>%
    mutate(bpm = 60/(msec/mean_peaks)*1000)
}

average_echo <- function(data, measures = c(AWTd, EDD, PWTd, AWTs, ESD, PWTs)){ # animal id in one column, AWTd, EDD, PWTd, AWTs, ESD and PWTs in columns 
    data %>%
        group_by(animal_id) %>%
        summarise(across(measures, ~ mean(.x, na.rm = TRUE)))
}

calculations <- function(data){ # ASE cube = 1.04*((IVSTd+LVIDd+PWTd)3–LVIDd3) Devereux correction = 0.8(ASEcube)+0.6
  data %>%
    mutate(ASE = 1.04*(((AWTd*10)+(EDD*10)+(PWTd*10))^3 - (EDD*10)^3)) %>%
    mutate(LVM = (0.8*(ASE)+0.6)/1000) %>%
    mutate(FS = (EDD - ESD)/EDD*100) %>%
    mutate(RWT = ((AWTd*10)+(PWTd*10))/(EDD*10)) %>% 
    mutate(EDDV = (EDD^3)*1.047, 
           ESDV = (ESD^3)*1.047) %>%
    mutate(SV = (EDDV - ESDV)/1000) %>%
    mutate(CO = (SV * bpm)/1000) %>%
    mutate(EF = SV/EDDV)
}



