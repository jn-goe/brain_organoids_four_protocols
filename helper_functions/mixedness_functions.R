sample_stratified <- function(obj, group, size, df_bin = F) {
  
  if(!df_bin) { 
    df <- obj@meta.data
  }
  
  else {
    df <- obj
  }

  group.table <- unique(df[[group]])
  size.per.group <- round(size/length(group.table))
  group.table <- names(sort(table(df[[group]])))
  
  cells <- c()
  count <- 1
  
  for(g in group.table) {
    cells <- c(cells, 
               sample(rownames(df)[which(df[[group]] == g)], 
                      min(sum(df[[group]] == g), size.per.group)))
    if(sum(df[[group]] == g) < size.per.group) {
      size <- size - sum(df[[group]] == g)
      size.per.group <- round(size/(length(group.table) - count))
      count <- count + 1
    }
  }
  return(cells)
}
