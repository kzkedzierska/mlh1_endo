#!/usr/bin/env Rscript
library("tidyverse")
library("optparse")

# Test for other dependecies
not_loaded_packages <- c("biomaRt", "data.table") 
test_not_loaded_packages <- not_loaded_packages %in% installed.packages()

in (any(!test_not_loaded_packages)) {
  error_message <- paste0("Missing the following package(s): ",
                          paste(not_loaded_packages[!test_not_loaded_packages], 
                                collapse = " "),
                          ". Please install it/them and rerun the script.")
  stop(error_message)
}

# Arguments
arguments_list = list(
  make_option(c("-m", "--manifest"), 
              type = "character", 
              default = NULL, 
              help = "Path to the manifest with expression&methylation files.", 
              metavar = "character"),
  make_option(c("-e", "--expression"), 
              type = "character", 
              default = NULL, 
              help = "Directry with expression files", 
              metavar = "character"),
  make_option(c("-t", "--methylation"), 
              type = "character", 
              default = NULL, 
              help = "Directory with methylation files.", 
              metavar = "character"),
  make_option(c("-g", "--genotype_manifest"), 
              type = "character", 
              default = NULL, 
              help = "Path to the manifest with raw WES sequencing.", 
              metavar = "character"),
  make_option(c("-s", "--status"), 
              type = "character", 
              default = NULL, 
              help = "Path to file with Microsatelite instability status.", 
              metavar = "character"),
  make_option(c("-p", "--platypus"), 
              type = "character", 
              default = NULL, 
              help = "Path to the directory with platypus calls", 
              metavar = "character")
) 

args_parser <- OptionParser(option_list = arguments_list)
args <- parse_args(args_parser);

# Test arguments
args_test <- lapply(args, function(path) {
  check_path <- file.exists(path)
  if (!check_path) {
    print(paste0("Path: ", path, " is incorrect!"))
  }
  check_path
})

stopifnot(args_test)

# assign args to variables
manifest_path <- args$manifest
manifest_vcf_path <- args$genotype_manifest
mss_status_path <- args$status
expression_path <- args$expression
methylation_path <- args$methylation

# Hard coded constants
snp_position <- 36993455

# Read in data
manifest <- data.table::fread(manifest_path)
manifest_vcf <- data.table::fread(manifest_vcf_path)
mss_status <- data.table::fread(mss_status_path)

# Prepare text file to download sample table
format <- "TSV"
fields <- c("file_name", 
            "cases.samples.portions.analytes.aliquots.submitter_id",
            "data_category")


Part1 <- '{"filters":{"op":"in","content":{"field":"files.file_id","value":[ '

Part2 <- paste0('] }},', 
                '"format":"', format, '",',
                '"fields":"', paste(fields, collapse = ","), '",',
                '"size":')

id <- toString(sprintf('"%s"', manifest$id))
Part3 <- paste0('"', nrow(manifest), '"}')
Sentence <- paste(Part1, id, Part2, Part3, collapse = " ")
write.table(Sentence, "Payload.txt",
            quote = FALSE, col.names = FALSE, row.names = FALSE)

# Download sample table
cmd <- paste0('curl --request POST --header "Content-Type: application/json" --data @', 
              'Payload.txt',
              '"https://api.gdc.cancer.gov/files"')
sample_table <- data.table::fread(cmd) 

colnames(sample_table) <- c("file_name", "file_id", "data_category", "id")
sample_table <- sample_table %>%
  mutate(sample_id = gsub(pattern = "(TCGA-[^-]*-[^-]*)-.*", 
                          replacement = "\\1", 
                          file_id),
         sample_type = ifelse(grepl(pattern = "0[1-9][A-Z]", 
                                    substr(file_id, 14, 16)),
                              "Tumor",
                              "Control"))
  

# Read in expression files
ensembl <- biomaRt::useMart("ensembl", 
                           dataset = "hsapiens_gene_ensembl")
gene_id_name <- biomaRt::getBM(attributes = c("ensembl_gene_id", 
                                              "external_gene_name"),
                               mart = ensembl)
colnames(gene_id_name) <- c("gene_id", "gene_name")

read_expr <- function(fname, 
                      dir = expression_path, gene = "MLH1") {
  bname <- basename(fname)
  if (grepl(".gz$", fname)) {
    fname <- paste0("gunzip -c ", fname)
    } 
  data.table::fread(file.path(dir, fname)) %>%
    mutate(file_name = bname) %>%
    select(gene_id = V1, expression = V2, file_name) %>% 
    mutate(gene_id = gsub(x = gene_id, pattern = "\\..*", # not interested with gene_id versions
                          replacement = "")) %>%
    left_join(gene_id_name) %>%
    filter(gene_name == gene) %>%
    select(gene_name, expression, file_name)
} 

# Extract expression file names and read in expression
expression_files <- sample_table %>%
  filter(data_category == "Transcriptome Profiling") %>%
  mutate(ptf = file.path(id, file_name)) %>% 
  pull(ptf)

expression_df <- lapply(expression_files, read_expr) %>%
  bind_rows() %>%
  left_join(sample_table, by = "file_name") 


# Download methylation data
methylation_files <-  sample_table %>%
  filter(data_category == "DNA Methylation") %>%
  mutate(ptf = file.path(id, file_name)) %>% 
  pull(ptf)

read_meth <- function(fname, 
                      dir = methylation_path, gene_name = "MLH1") {
  bname <- basename(fname)
  fname <- file.path(dir, fname)
  if (grepl(".gz$", process_cmd)) {
    # adding tail -n+2 so column names will be the same in all and gene name
    process_cmd <- paste0("gunzip -c ", fname, "| tail -n+2")
  } else {
    process_cmd <- paste0("tail -n+2 ", fname) 
    } 
  if (gene != "all") {
    process_cmd <- paste0(process_cmd, "| awk '$6 ~ \"", gene_name, "\"'")
  }
  data.table::fread(process_cmd) %>%
    select(probe_id = V1,
           beta = V2,
           gene = V6) %>%
    separate_rows(gene, sep = ";") %>%
    distinct(probe_id, beta, gene) %>%
    filter(gene == gene_name) %>%
    mutate(file_name = bname)
}

methylation_df <- lapply(methylation_files, read_meth) %>%
  bind_rows() %>%
  left_join(sample_table, by = "file_name") 

# Processing the genotype calls
vcf_process_cmd <- paste0("bash process_vcf.sh ", vcf_path)
data.table::fread(vcf_process_cmd)
genotype <- data.table::fread(vcf_process_cmd) %>% 
  rename(id = V1, 
         position = V2, 
         ref_allele = V3, alt_allele = V4, 
         total_coverage = V5, alt_reads = V6)



# Prepare text file to download sample table
format <- "TSV"
fields <- c("file_name", 
            "cases.samples.portions.analytes.aliquots.submitter_id",
            "data_category")


Part1 <- '{"filters":{"op":"in","content":{"field":"files.file_id","value":[ '

Part2 <- paste0('] }},', 
                '"format":"', format, '",',
                '"fields":"', paste(fields, collapse = ","), '",',
                '"size":')

id_vcf = toString(sprintf('"%s"', manifest_vcf$id))
Part3_vcf = paste0('"', nrow(manifest_vcf), '"}')
Sentence <- paste(Part1, id_vcf, Part2, Part3_vcf, collapse = " ")
write.table(Sentence, "Payload_vcf.txt",
            quote = F, col.names = F, row.names = F)

cmd_vcf <- paste0('curl --request POST --header ',
                  '"Content-Type: application/json" --data @', 
                  'Payload_vcf.txt',
                  ' "https://api.gdc.cancer.gov/files"')
sample_table_vcf <- data.table::fread(cmd_vcf) 

colnames(sample_table_vcf) <- c("file_name", "file_id", "data_category", "id2")
sample_table_vcf <- sample_table_vcf %>%
  mutate(sample_id = gsub(pattern = "(TCGA-[^-]*-[^-]*)-.*", 
                          replacement = "\\1", 
                          file_id),
         sample_type = ifelse(grepl(pattern = "0[1-9][A-Z]", 
                                    substr(file_id, 14, 16)),
                              "Somatic",
                              "Germline"),
         id = gsub("_.*", "", file_name))


genotype2 <- genotype %>% 
  left_join(sample_table_vcf) %>% 
  filter(position == snp_position) %>% 
  mutate(fract = alt_reads / total_coverage, 
         genotype = ifelse(is.na(fract), 
                           NA,
                           ifelse(fract > .2 & fract < .8, 
                                  paste0(ref_allele, alt_allele), 
                                  ifelse(fract >= .8, 
                                         paste0(alt_allele, alt_allele),
                                         paste(ref_allele, ref_allele))))) %>%
  group_by(sample_id, sample_type) %>% 
  summarise(genotype = paste(unique(genotype[!is.na(genotype)]), collapse=",")) %>%
  mutate(genotype = ifelse(genotype == "", "GG", genotype)) %>%
  spread(sample_type, genotype)


# Final df created by joining all subsequent dataframes
full_df <- expression_df %>% 
  select(sample_id, expression, sample_type) %>% 
  full_join((methylation_df %>% 
               select(probe_id, beta, sample_id, sample_type) %>%
               full_join(genotype2))) %>% 
  full_join(mss_status) 

