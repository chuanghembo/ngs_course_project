library(tidyverse)
library(biomaRt)

WT_vs_WTLeukemia <- read.csv('../WT_vs_WTLeukemia.csv', row.names = 'X') 
WT_vs_PPM1DTLeukemia <- read.csv('../WT_vs_PPM1DTLeukemia.csv', row.names = 'X')
PPM1DTLeukemia_vs_WTLeukemia <- read.csv('../PPM1DTLeukemia_vs_WTLeukemia_clean.csv')

# Get ensembl id
ensembl_id <- WT_vs_WTLeukemia |> 
        pull(row) |> 
        as.character()

# Specify Mart
ensembl = useMart("ensembl", dataset = 'mmusculus_gene_ensembl')

# Submit Ensembl id in chunks and merge the result in the end
# Ensembl can take around 500-1000 ensembl ids at a time
chunk_size <- 100

id_chunks <- ensembl_id[1:100] |>  
        enframe() |>  
        mutate(chunk = ceiling(row_number() / chunk_size)) |> 
        nest(data = c(value))

# Function for batch querying
batch_query <- function(ids, mart) {
        getBM(attributes = c('ensembl_gene_id_version', 
                             'external_gene_name',
                             'go_id',
                             'name_1006',
                             'definition_1006'),
              filters = 'ensembl_gene_id_version', 
              values = ids, 
              mart = mart)
}

# Perform the batch querying
gene_info <- id_chunks |> 
        mutate(results = map(data, ~batch_query(.x$value, mart = ensembl))) |> 
        dplyr::select(results) |> 
        unnest(cols = c(results))

# Write the results to a CSV file
write.csv(gene_info, "../biomart_query_results.csv", row.names = FALSE)