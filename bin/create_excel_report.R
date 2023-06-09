#!/usr/bin/env Rscript
# rm(list = ls())
library(tidyverse)
library(data.table)
library(glue)
library(openxlsx)
library(jsonlite)

args <- commandArgs(trailingOnly = TRUE)
root <- args[1]
xls_name <- args[2]

# Functions

# Define a function called "rr" that takes a file path as an argument
rr <- function(x) {
  # Read the 7th line of the file at the given path, and then create a tibble with two columns: "sample" and "line"
  read_lines(x)[7] %>%
    tibble(sample = basename(x) %>% str_remove(., ".flagstat"), line = .) %>%
    # Extract the number of human reads and human proportion from the "line" column, and add them as new columns
    mutate(human_reads = str_extract(line,"^[0-9]{1,10}")) %>%
    mutate(human_proportion = str_extract(line,"[0-9]{1,2}\\.[0-9]{1,3}")) %>%
    # Remove the "line" column
    select(-line)
}

# Define a function called "read_amr" that takes a file path as an argument
read_amr <- function(x) {
  # Read the file at the given path using the "fread" function from the "data.table" package
  tmp <- fread(x)
  # Check if the resulting data frame has any rows
  if (nrow(tmp) > 0) {
    # If the data frame is not empty, add several new columns to it using the "mutate" function from the "dplyr" package
    tmp %>%
      # Add a "sample" column that contains the file name without the ".res" extension
      mutate(sample = basename(x) %>% str_remove(., ".res")) %>%
      # Split the "#Template" column into eight parts using the "str_split_fixed" function from the "stringr" package
      mutate(split = str_split_fixed(`#Template`,"\\|",8)) %>%
      # Extract the sixth, seventh, and eighth parts of the "#Template" column, and use them to create new columns called "AMR_gene", "AMR_gene_var", and "AMR_desc", respectively
      mutate(AMR_gene = split[,6]) %>%
      mutate(AMR_gene_var = split[,7]) %>%
      mutate(AMR_desc = split[,8]) %>%
      # Remove the "split" column
      select(-split) %>%
      # Move the "sample" column to the front of the data frame
      select(sample, everything())
  }
}

# Define a function called "read_taxa" that takes a file path as an argument
read_taxa <- function(x) {
  # Read the file at the given path using the "fread" function from the "data.table" package, selecting the first three columns and skipping the first row
  tmp <- fread(x, select = c(1:3), skip = 1)
  # Check if the resulting data frame has more than one row
  if (nrow(tmp) > 1) {
    # If the data frame has more than one row, add a "sample" column to it using the "mutate" function from the "dplyr" package
    tmp %>%
      # Add a "sample" column that contains the file name without the "_metaphlan_report.csv" extension
      mutate(sample = basename(x) %>% str_remove(., "_metaphlan_report.csv")) %>%
      # Move the "sample" column to the front of the data frame
      select(sample, everything())
  } else {
    # If the data frame has only one row or is empty, return NULL
    NULL
  }
}

read_json <- function(json_file){
  json_data <- fromJSON(json_file)
  # Create a data frame with the required columns
  df <- data.frame(reads = json_data$reads,
                   bases = json_data$bases,
                   n50 = json_data$n50,
                   longest = json_data$longest,
                   shortest = json_data$shortest,
                   mean_length = json_data$mean_length,
                   median_length = json_data$median_length,
                   mean_quality = json_data$mean_quality,
                   median_quality = json_data$median_quality,
                   top_lengths = I(list(json_data$top_lengths)[[1]] %>% paste(., collapse=",")),
                   top_qualities = I(list(json_data$top_qualities)[[1]] %>% paste(., collapse=","))
  ) %>%
    mutate(sample = basename(json_file)) %>%
    mutate(filter = ifelse(grepl("pre", sample), "pre", "post")) %>%
    mutate(sample = str_remove(sample, "\\.stats\\.(pre|post)\\_filtered\\.json")) %>%
    select(sample, filter, everything())
  return(df)
}


flagstats <- list.files(root, "*.flagstat", full.names = T)

if (length(flagstats) > 0) {
  human_reads <- map_df(flagstats,rr)
} else {
  human_reads <- tibble()
}


# AMR

amr_files <- list.files(root, "*.res", full.names = T)
if (length(amr_files) > 0) {
  amr <- map_df(amr_files, read_amr)
} else {
  amr <- tibble()
}




read_stat_files <- list.files(root,"*.json", full.names=T)

if (length(read_stat_files) > 0) {
  read_stats <- map_df(read_stat_files, read_json)
} else {
  read_stats <- tibble()
}



taxa_files <- list.files(root, "*_metaphlan_report.csv", full.names = T)

if (length(taxa_files) > 0) {
  taxa <- map_df(taxa_files, read_taxa) %>%
    pivot_wider(names_from = sample, values_from = relative_abundance) %>%
    mutate(rank_level0 = str_count(`#clade_name`,"\\|")) %>%
    mutate(rank_level = case_when(
      rank_level0 == "0" ~ "kingdom",
      rank_level0 == "1" ~ "phylum",
      rank_level0 == "2" ~ "class",
      rank_level0 == "3" ~ "order",
      rank_level0 == "4" ~ "family",
      rank_level0 == "5" ~ "genus",
      rank_level0 == "6" ~ "species",
      rank_level0 == "7" ~ "strain",
      TRUE ~ ""
    )) %>%
    select(-rank_level0)
} else {
  taxa <- tibble()
}


xls <- list()

for (df in c("read_stats", "human_reads", "taxa", "amr")) {
  if (nrow(get(df)) > 0 ) {
    xls[[toupper(df)]] <- get(df)
  }
}

write.xlsx(xls,file = glue("{root}/{xls_name}.xlsx"),asTable = T,firstRow=T)


