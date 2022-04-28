date_col  <- function(x, 
                  n){
  rep(x, n)
}
age_col <- function(x, 
                    n){
  paste(c(rep(x, 
              n)))
}
sex_col <- function(n_m, 
                    n_f){
  c(rep("male", n_m), rep("female", n_f))
}

id_mother <- function(mother, 
                      n){
  paste(c(rep(mother, 
              n)))
}

id_neo <- function(mother,
                   seq_start,
                   n, 
                   age){
  sequence <- seq(seq_start, (n+seq_start-1), by=1) 
  
  paste(mother, 
        age, 
        sequence, 
        sep = "_")
}

id_5wk <- function(strain, 
                   start, 
                   n){
  n=n-1
  
  sequence <- seq(start, (start+n), by=1)
  
  paste(strain, 
        sequence, 
        sep = "")
}

read.data <- function(files, 
                      sep = ","){
 
  tables = lapply(files, function(x){
    file = x
    tables = read.table(file = file, 
                        header = T, sep = sep, 
                        stringsAsFactors = T)})
  
  names(tables) = files
  
  return(tables)
}

wt2ratio <- function(x, y){
  x/y
}

not_all_na <- function(x) any(!is.na(x))

stats <- function(data, 
                  dv, 
                  iv){
  data %>%
    anova_test(dv ~ iv)%>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance("p.adj")
}

plot_data <- function(dat, 
                      measure, 
                      comparisons){
  ggbarplot(dat, 
            x = "genotype", 
            y = paste(measure), 
            add = c("mean_se", "jitter"), 
            fill = "genotype")
}
  