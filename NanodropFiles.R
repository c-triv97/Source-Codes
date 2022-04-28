load.files = function(pattern, sep = "\t", 
                      write.table = TRUE, 
                      filename = NA){
  require(data.table)
  require(dplyr)
  
  files = list.files(pattern = pattern, 
                     recursive = T)
  tables = lapply(files, function(x){
    file = x
    tables = read.table(file = file, 
                        header = T, sep = sep)})
  
  names(tables) = files
  
  tables2 = data.table::rbindlist(tables, idcol = "filename")%>%
    dplyr::select(c("filename", 
                    "Sample.ID", 
                    "ng.ul", 
                    "X260.280", 
                    "X260.230", 
                    "Date"))
  
  if (write.table == TRUE){
    write.csv(tables2, file = paste("Nanodrop", filename, "csv", 
                                    sep = "."),
              row.names = FALSE)
  }
  
  return(tables2)
}
