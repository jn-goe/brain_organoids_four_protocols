map_symbols <- function(ensembl, map) {
  symbols = vector()
  i = 1
  for(id in ensembl) {
    symbols[i] = map$symbol[map$ensembl == id]
    i = i+1
  }
  symbols
}

map_symbols_rev <- function(gene, map) {
  ensembl <- gene_symbols[which(gene_symbols[,2] == gene),1]
  return(ensembl)
}