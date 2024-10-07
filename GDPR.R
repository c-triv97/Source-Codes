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

wd <- function(file_path, sub){
  wd = gsub(sub, "", getwd())
  dir =  paste0(wd, file_path, sep = "")
  setwd(dir)
}

`%!in%` <- negate(`%in%`)

# reading multiple excel sheets
require("readxl") 
multiplesheets <- function(fname) { 
   
  # getting info about all excel sheets 
  sheets <- readxl::excel_sheets(fname) 
  tibble <- lapply(sheets, function(x) readxl::read_excel(fname, sheet = x)) 
  data_frame <- lapply(tibble, as.data.frame) 
    
  # assigning names to data frames 
  names(data_frame) <- sheets 

  return(data_frame) 
} 
