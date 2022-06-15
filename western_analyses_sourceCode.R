desktop <- "C:/Users/2454683t/"
pc <- "C:/Users/carat/"
onedrive <- "OneDrive - University of Glasgow/PhD/"
# source code for analysing westerns 
find_files <- function(pat, replace, samples){
  files = list.files(pattern = pat, 
                     recursive = T)
  
  data = list()
  for (i in files){
    d = read_excel(path = i, 
                   col_names = T)
    data[[i]] = d
  }
  
  names(data) = c(str_replace(files, replace, ""))
  
  data = lapply(data, function(x){
    
    lanes = paste("lane", seq(1, length(samples)))
    
    x %>%
      mutate(Name = samples, 
             lane = lanes)%>%
      dplyr::select(c(Name, Signal, lane))
  })
}

normalisation <- function(dat){
  data1 = bind_rows(dat, 
                    .id = "Result")
  
  max = max(data1$Signal)
  
  LNF = subset(data1, 
               Result == "TPS")$Signal / max
  
  data = lapply(dat, function(x){
    x %>%
      mutate(LNF = LNF, 
             NS = Signal / LNF)%>%
      mutate_if(is.character, as.factor)
})
  
  return(data)
}


stat_results <- function(dat, groups, signal){
  
  dat = dat %>%
    mutate(condition = groups,
           NS = signal)
  
  stat.test = dat %>% 
    anova_test(NS ~ condition) %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance("p.adj")
  
  print(stat.test)
  
  pwc = dat %>%
    tukey_hsd(NS ~ condition) %>%
    add_xy_position(x = "condition")
  
  return(pwc)
} 

plot_gen <- function(dat, groups, pwc){
  condition = groups 
  
  data = lapply(dat, function(x){
    x %>%
      mutate(condition = condition)%>%
      mutate(lane = as.character(lane))
  })
  
  names(data) = c("protein", "raw")
  
  plotProtein = ggbarplot(data$protein, 
                          x = "condition", 
                          y = "NS", 
                          fill = "condition",
                          color = "condition",
                          add = c("mean_se", "jitter"),
                          add.params = list(color = "black", fill = "black"),
                          xlab = "", 
                          ylab = "Normalised Signal", 
                          ggtheme = theme_pubr(base_size = 20), 
                          legend = "right", 
                          palette = "Dark2")+
    stat_pvalue_manual(pwc, hide.ns = T)
  
  plotRawProtein = ggbarplot(data$protein,
                             x = "lane",
                             y = "Signal",
                             fill = "condition",
                             color = "condition",
                             xlab = "",
                             ylab = "Raw Signal",
                             palette = "Dark2", 
                             ggtheme = theme_pubr(), 
                             legend ="none")%>%
    ggpar(x.text.angle = 45)
  
  plotRaw = ggbarplot(data$raw, 
                      x = "lane", 
                      y = "Signal",
                      fill = "condition",
                      color = "condition",
                      palette = "Dark2",
                      xlab = "", 
                      ylab = "Total Protein Signal",
                      ggtheme = theme_pubr(), 
                      legend = "none")%>%
    ggpar(x.text.angle = 45)
  
  plot = ggarrange(plotProtein,
                   ggarrange(plotRaw, plotRawProtein, 
                             ncol = 2, labels = c("B", "C")), 
                   nrow=2, 
                   labels = "A", 
                   heights = c(2,1))
  
  return(plot)
}

z_score <- function(data.tab, value){
  dat = data.tab %>%
    mutate_at(vars(value), funs(z = (. - min(.)) / (max(.) - min(.))))
  
  return(dat)
}


