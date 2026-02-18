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

# mixed order ordering function so that SPL1 to SPL20 is logically ordered by alphabet and then number
require("gtools") 
mixedrank = function(x) order(gtools::mixedorder(x))

# group colours 
require("colorspace")
require("RColorBrewer") 
# Function: generate replicate colors with variable replicates per group
make_replicate_palette <- function(n_reps, base_colors = group_colors) {
  n_groups <- length(n_reps)
  pal_list <- vector("list", n_groups)
  for (i in seq_len(n_groups)) {
    pal_list[[i]] <- rev(colorRampPalette(c("white", base_colors[i]))(n_reps[i] + 1)[-1])
  }
  unlist(pal_list)
}
                   
