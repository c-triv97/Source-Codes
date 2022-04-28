find_files <- function(pat){
  require(stringr)
  require(dplyr)
  
  files = list.files(pattern = pat, 
                     recursive = T)
  
  data = list()
  for (i in files){
    d = read.table(i, 
                   sep = ",", 
                   header = T)
    d = d %>%
      dplyr::select(c(Length))%>%
      mutate(measure = c(rep(c("PWTd",
                             "EDD",
                             "AWTd",
                             "PWTs",
                             "ESD",
                             "AWTs"), 3)))
    data[[i]] = d
  }
  
  names(data) = c(str_replace(files, ".csv", ""))
  
  return(data)
}

calculations = function(data){
  d = spread(data, key = measure, value = length)
  
  d1 = d %>%
    mutate(calc = (AWTd*10)+(PWTd*10)+(EDD*10),
      ASE.cube = ((calc^3)-((EDD*10)^3))*1.04,
      LVM = (0.8*ASE.cube+0.6)/1000,
      FS = (EDD-ESD)/EDD*100,
      RWT = ((AWTd*10)+(PWTd*10))/(EDD*10),
      SV = (1.047*(EDD^3))-(1.047*(ESD^3)), 
      EF.pct.1 = (SV/(1.047*(EDD^3)))*100,
      area.d = 7/(2.4 + EDD)*(EDD^3),
      EF.pct.2 = SV/area.d)
  
  return(d1)
}

corrections <- function(data, msec = 1286.67){
  data %>%
    mutate(hr = 60/(msec/peaks), 
           bpm = hr * 1000, 
           CO = SV*bpm,
           SV.corr = SV/1000, 
           LVM.corr =LVM/1000, 
           CO.corr =  CO/1000)
}
