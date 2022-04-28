## source code for BCA assays 

sample_gen <- function(strain, 
                       start,
                       n){
  n=n-1
  
  sequence <- seq(start, (start+n), by=1)
  
  paste(strain, 
        sequence, 
        sep = "-")
}

standard_gen <- function(type, 
                         replicates, 
                         max_col){
  matrix(rep(type, replicates), ncol = max_col, 
         byrow = T)
}

row_gen <- function(sample,
                    replicates){
  matrix(c(rep(sample[1], replicates),
           NA, 
           rep(sample[2], replicates), 
           NA, 
           rep(sample[3], replicates)), 
              nrow =1)
}

plate2graph <- function(x, 
                        row){
  x %>%
    mutate(row = 1:row) %>%
    pivot_longer(-row, names_to = "col", values_to = "value")%>%
    mutate(col = as.integer(str_remove(col, "V")))%>%
    mutate(well = paste0(LETTERS[row],0, col))
}

plot_96well <- function(x, y){
  ggplot(data = x)+
    geom_circle(aes(x0 = col, 
                    y0 = row, 
                    r= 0.5, 
                    fill = value))+
    coord_equal()+
    scale_x_continuous(breaks = 1:12, 
                       expand = expansion(mult = c(0.01, 0.01)))+
    scale_y_continuous(breaks = 1:8, labels = LETTERS[1:8], 
                       expand = expansion(mult = c(0.01, 0.01)), 
                       trans = reverse_trans())+
    labs(title = "96 well plate", 
         subtitle = paste0(y),
         x="Col", y="Row")+
    theme_bw(base_size =20)+
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          legend.position = "none")+
    geom_text(aes(x = col, 
                  y = row, 
                  label = paste0(value)), size =3)
}

plot_96wellABS <- function(x, measure, y){
  ggplot(data = x)+
    geom_circle(aes(x0 = col, 
                    y0 = row, 
                    r= 0.5, 
                    fill = measure))+
    coord_equal()+
    scale_x_continuous(breaks = 1:12, 
                       expand = expansion(mult = c(0.01, 0.01)))+
    scale_y_continuous(breaks = 1:8, labels = LETTERS[1:8], 
                       expand = expansion(mult = c(0.01, 0.01)), 
                       trans = reverse_trans())+
    labs(title = "96 well plate", 
         subtitle = paste0(y),
         x="Col", y="Row")+
    theme_bw(base_size =20)+
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          legend.position = "right")+
    geom_text(aes(x = col, 
                  y = row, 
                  label = paste0(value)), size =3)+
    scale_fill_gradientn(colours=rainbow(4))
}

curve_fit <- function(x, 
                      y){
  require(ggplot2)
  
  ggplot(x[-1,], 
         aes(x = as.numeric(value), 
             y = mean_abs))+
    geom_point(size = 3)+
    geom_smooth(method = drm, 
                method.args = list(fct = LL.4()), 
                se = FALSE, 
                colour = "grey")+
    theme_pubr(legend = "right")+
    geom_point(data = y, aes(x = estimate, 
                             y = mean_abs,
                             colour = value), 
               shape = 18, 
               size =4)+
    labs(title = "Standard Curve and Unknown Samples", 
         x = "concentration (ug/ml)", 
         y = "absorbance (A)")+
    scale_x_log10()
}


final_conc <- function(data, 
                       df, 
                       tube){
  data %>%
    mutate(`final conc (ug/ml)` = estimate * df)%>%
    mutate(tissue = tube)
}

western_vols <- function(data, 
                         ug, 
                         vol){
  data %>%
    mutate(`vol req (uL)` = (ug / (`final conc (ug/ml)` / 1) * 1000))%>%
    mutate(`diluent (uL)` = vol - `vol req (uL)`)
}


data_import <- function(pat, plate.dat, replicates, path=getwd()){
  
  data.file = list.files(pattern = pat, 
                         recursive = T, 
                         path = getwd())
  data = read_excel(data.file)
  
  colnames(data)[grepl('Absorbance',colnames(data))] = "absorbance"
  
  wells = plate.dat[!(is.na(plate.dat$value)),] %>%
    pull(well)
  
  data = data %>%
    mutate(well = wells)
  
  data1 = merge(data, 
                 plate.dat, 
                 all =T)
  
  blank = subset(data1, value == 0)%>%
    pull("absorbance")  
  
  blank = blank/replicates
  
  data2 <- data1 %>%
    mutate(corrected_absorbance = absorbance - blank) 
  
  data3 <- data2 %>%
    group_by(value)%>%
    summarise_at(vars(corrected_absorbance), 
                 list(mean_abs = mean, 
                      sd = sd))
  
  data.list = list(raw_data = data1, summary_data = data3)
  
  return(data.list)
}
 
fit2curve <- function(summary.dat, standards, samples, DF, tissue, ug, vol){

  standard_curve = subset(summary.dat, 
                           value %in% standards) #subsetting a table to use for standard curve generation
  sample_data = subset(summary.dat,
                        value %in% samples) #subetting table of unknowns for estimation of protein conc using model created by standards
  
  model = drc::drm(mean_abs ~ as.numeric(value), 
                    fct = LL.4(), 
                    data = standard_curve) #4-parameter quadratic fit of data 
  
  concentration = ED(model,
                     sample_data$mean_abs, 
                      type = "absolute", 
                      display = F)%>%
    as.data.frame() #generating estimate protein concs of the samples based on model generated by known standards 
  
  sample_data = sample_data %>%
    mutate(estimate = concentration$Estimate, 
           std_error = concentration$`Std. Error`) #adding estimate and SD to the sample_data table 
  
  fit = curve_fit(standard_curve,
                  sample_data)
  
  sample_data = final_conc(sample_data, 
                           df = DF, 
                           paste(tissue))  # applying dilution factor and annotating with some more sample information 
  
  #generating required ul for given ug and diluent for given uL total volume excluding SDS sample buffer  
  sample_data = western_vols(sample_data, 
                             ug = ug, 
                             vol = vol)
  
  list.model = list(drc.fit = model, 
                    fitted.model = fit, 
                    protein.conc = sample_data, 
                    standards = standard_curve)
  
  return(list.model)
}

  
  
  
  
  
  