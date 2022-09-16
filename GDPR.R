find_files <- function(pat, replace, path, sep, skip){
  require(stringr)
  setwd(path)
  
  files = list.files(pattern = pat, 
                     recursive = T)
  
  data = list()
  for (i in files){
    d = read.table(file = i,
                   sep = sep, 
                   skip = skip, 
                   header = TRUE)
    data[[i]] = d
  }
  
  names(data) = c(str_replace(files, replace, ""))
  
  return(data)
}

standard_error <- function(x) sd(x, na.rm = TRUE) / sqrt(length(na.omit(x))) 

wd <- function(file_path = "PhD/", sub="PhD.*"){
  wd = gsub(sub, "", getwd())
  dir =  paste0(wd, file_path)
  setwd(dir)
}

`%!in%` <- negate(`%in%`)
