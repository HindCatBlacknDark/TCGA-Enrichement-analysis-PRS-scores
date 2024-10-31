######## Set up######

# %% BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("TCGAbiolinks")
BiocManager::install("SummarizedExperiment")
BiocManager::install("GSVA")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("GSEABase")
BiocManager::install("AnnotationDbi")

# %%
if (!requireNamespace("pheatmap", quietly = TRUE)) {
    install.packages("pheatmap")
}

# %% packages
install.packages("dplyr")
install.packages("tidyverse")
install.packages("purrr")
install.packages("readxl")
install.packages("writexl")
install.packages("parallel")


# %% libraries
library(TCGAbiolinks)
library(survival)
library(tidyverse)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(GSVA)
library(parallel)
library(clusterProfiler)
library(org.Hs.eg.db)
library(readxl)
library(writexl)
library(dplyr)
library(biomaRt)
library(purrr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GSEABase)
library(AnnotationDbi)
library(pheatmap)
library(data.table)


# %% Set up parallel processing
n_cores <- detectCores() - 1
cl <- makeCluster(n_cores)
clusterEvalQ(cl, {
  library(DESeq2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

# %% Data Query for Gene Expression Data
query_expr <- GDCquery(
  project = "TCGA-LUAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)
GDCdownload(query_expr)
se <- GDCprepare(query_expr)

# Filter for protein-coding genes and save
se_mrna <- se[rowData(se)$gene_type == "protein_coding", ]
saveRDS(se_mrna, file = "TCGA_LUAD_mRNA_protein_coding.rds")
saveRDS(se, file = "TCGA_LUAD_complete.rds")

# Save clinical data
clinical_data <- colData(se)
saveRDS(clinical_data, file = "TCGA_LUAD_clinical_data.rds")

# Output dimensions and sample data for verification
cat("Expression data dimensions:", dim(se_mrna), "\n")
cat("Clinical data dimensions:", dim(clinical_data), "\n")
cat("Head of clinical data:\n")
print(head(clinical_data[, 1:66]))

# %% Data Query for SNP Data
query_snp <- GDCquery(
  project = "TCGA-LUAD",
  data.category = "Simple Nucleotide Variation",
  data.type = "Masked Somatic Mutation"
)
GDCdownload(query_snp)
snp_data <- GDCprepare(query_snp)
saveRDS(snp_data, file = "TCGA_LUAD_SNP_data.rds")

# Output dimensions and a sample for verification
cat("SNP data dimensions:", dim(snp_data), "\n")
print(head(snp_data))

#################### %% clinical data #########################
# Load the clinical data
#clinical_data <- readRDS("TCGA_LUAD_clinical_data.rds")

# Convert clinical_data to a data frame
clinical_data_df <- as.data.frame(clinical_data)

# Extract relevant columns
selected_clinical_data <- clinical_data_df %>%
  dplyr::select(sample_submitter_id, shortLetterCode, definition, sample_type, 
                race, gender, age_at_diagnosis, days_to_birth, tissue_or_organ_of_origin, 
                vital_status, disease_type, primary_site, year_of_diagnosis)

# Display dimensions of the selected clinical data and preview the data
cat("Selected Clinical Data Dimensions:", dim(selected_clinical_data), "\n")
print(head(selected_clinical_data, 10))

# Filter to keep only primary tumor samples (where shortLetterCode == "TP")
primary_tumors <- selected_clinical_data %>% filter(shortLetterCode == "TP")

# Summary statistics for primary tumors
cat("Primary Tumor Data Dimensions:", dim(primary_tumors), "\n")
cat("Gender Distribution:\n")
print(table(primary_tumors$gender))
cat("Race Distribution:\n")
print(table(primary_tumors$race))
cat("Vital Status:\n")
print(table(primary_tumors$vital_status))

# Save the primary tumor data
saveRDS(primary_tumors, file = "TCGA_LUAD_primary_tumor_clinical_data.rds")




############# %% Preprocessing Expression Data #############
# Load gene expression data
se_mrna <- readRDS("TCGA_LUAD_mRNA_protein_coding.rds")

# Extract TPM, FPKM, and unstranded counts and convert Ensembl IDs to Hugo gene symbols
expr_tpm <- assay(se_mrna, "tpm_unstrand")
expr_fpkm <- assay(se_mrna, "fpkm_unstrand")
expr_counts <- assay(se_mrna, "unstranded")
symbol_mrna <- rowData(se_mrna)$gene_name
rownames(expr_tpm) <- symbol_mrna
rownames(expr_fpkm) <- symbol_mrna
rownames(expr_counts) <- symbol_mrna

# Remove rows with NA gene symbols
expr_tpm <- expr_tpm[!is.na(rownames(expr_tpm)), ]
expr_fpkm <- expr_fpkm[!is.na(rownames(expr_fpkm)), ]
expr_counts <- expr_counts[!is.na(rownames(expr_counts)), ]

# Aggregate duplicates by taking the mean for each dataset
aggregate_by_gene_symbol <- function(data) {
  aggregated_data <- aggregate(data, by = list(gene_symbol = rownames(data)), FUN = mean)
  rownames(aggregated_data) <- aggregated_data$gene_symbol
  aggregated_data$gene_symbol <- NULL
  return(aggregated_data)
}
expr_tpm <- aggregate_by_gene_symbol(expr_tpm)
expr_fpkm <- aggregate_by_gene_symbol(expr_fpkm)
expr_counts <- aggregate_by_gene_symbol(expr_counts)

# Filter low-expressed genes in TPM data
LUAD_tpm <- as_tibble(expr_tpm, rownames = "symbol_mrna") %>%
  mutate(meanrow = rowMeans(.[, -1], na.rm = TRUE), .before = 2) %>%
  filter(meanrow >= 1) %>%
  arrange(desc(meanrow)) %>%
  distinct(symbol_mrna, .keep_all = TRUE) %>%
  dplyr::select(-meanrow) %>%
  column_to_rownames(var = "symbol_mrna") %>%
  as.data.frame()

# Verify that gene symbols are still consistent after filtering
cat("First few row names of LUAD_tpm after filtering low-expressed genes:\n")
print(head(rownames(LUAD_tpm)))

# Define function to keep only primary tumor samples and standardize sample IDs
filter_primary_tumor_samples <- function(data) {
  colnames(data) <- substr(colnames(data), 1, 15)
  data <- data[, !duplicated(colnames(data))]
  data <- data[, substr(colnames(data), 14, 15) == "01"]
  return(data)
}

# Apply the primary tumor filtering to each dataset
LUAD_tpm <- filter_primary_tumor_samples(LUAD_tpm)
expr_fpkm <- filter_primary_tumor_samples(expr_fpkm)
expr_counts <- filter_primary_tumor_samples(expr_counts)

# Verify the filtered column names for primary tumor samples
cat("First few column names of LUAD_tpm after filtering for primary tumor samples:\n")
print(head(colnames(LUAD_tpm)))

# Remove "MT-" prefix from gene symbols in each dataset
rownames(LUAD_tpm) <- gsub("^MT-", "", rownames(LUAD_tpm))
rownames(expr_fpkm) <- gsub("^MT-", "", rownames(expr_fpkm))
rownames(expr_counts) <- gsub("^MT-", "", rownames(expr_counts))

# Output the first few rows and columns to confirm the structure
cat("First few rows and columns of cleaned TPM data:\n")
print(LUAD_tpm[1:5, 1:5])
cat("First few rows and columns of cleaned FPKM data:\n")
print(expr_fpkm[1:5, 1:5])
cat("First few rows and columns of cleaned Counts data:\n")
print(expr_counts[1:5, 1:5])

# Save each processed dataset for future use
saveRDS(LUAD_tpm, file = "TCGA_LUAD_TPM_cleaned.rds")
saveRDS(expr_fpkm, file = "TCGA_LUAD_FPKM_cleaned.rds")
saveRDS(expr_counts, file = "TCGA_LUAD_Counts_cleaned.rds")



############# CNV Data Preparation Phase #############
# %% Load and preprocess CNV data from GISTIC file
gistic_file <- "C:/Users/Hind.RAKI/hindcodes/TCGA_lung/datafiles/TCGA-LUAD.gistic/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes"
LUAD_cnv <- read.delim(gistic_file, header = TRUE, check.names = FALSE)

rownames(LUAD_cnv) <- LUAD_cnv$`Gene Symbol`
LUAD_cnv$`Gene Symbol` <- NULL

# Truncate sample IDs to the first 15 characters for standardization
colnames(LUAD_cnv) <- substr(colnames(LUAD_cnv), 1, 15)

# Filter to keep only primary tumor samples (sample type "01")
LUAD_cnv <- LUAD_cnv[, substr(colnames(LUAD_cnv), 14, 15) == "01"]

# Remove "MT-" prefix from gene symbols, NA, and duplicate rows
rownames(LUAD_cnv) <- gsub("^MT-", "", rownames(LUAD_cnv))
LUAD_cnv <- LUAD_cnv[!is.na(rownames(LUAD_cnv)) & rownames(LUAD_cnv) != "", ]
LUAD_cnv <- LUAD_cnv[!duplicated(rownames(LUAD_cnv)), ]

# Output dimensions and sample data for verification
cat("CNV data dimensions after cleaning:", dim(LUAD_cnv), "\n")
print(head(LUAD_cnv))

############# Align CNV and TPM Sample IDs #############
# %% Align sample IDs between TPM and CNV datasets
# Truncate sample IDs in the TPM dataset for consistency
colnames(expr_tpm) <- substr(colnames(expr_tpm), 1, 15)

# Find common sample IDs between TPM and CNV data
common_samples <- intersect(colnames(expr_tpm), colnames(LUAD_cnv))

# Subset TPM and CNV data to include only common sample IDs
expr_tpm_common <- expr_tpm[, common_samples]
cnv_common <- LUAD_cnv[, common_samples]

# Output dimensions of the subset data
cat("Aligned TPM data dimensions:", dim(expr_tpm_common), "\n")
cat("Aligned CNV data dimensions:", dim(cnv_common), "\n")

print(head(expr_tpm_common))
print(head(cnv_common))



############# Keep Only Common Genes between Expression and CNV Data #############
# %% Extract common genes and subset data accordingly
common_genes <- intersect(rownames(expr_tpm_common), rownames(cnv_common))
expr_tpm_common <- expr_tpm_common[common_genes, ]
cnv_common <- cnv_common[common_genes, ]

# Output dimensions of the subset data
cat("Final TPM data dimensions after gene alignment:", dim(expr_tpm_common), "\n")
cat("Final CNV data dimensions after gene alignment:", dim(cnv_common), "\n")
print(head(expr_tpm_common))
print(head(cnv_common))

############# Save Aligned TPM and CNV Data #############
# %% Save the aligned TPM data
saveRDS(expr_tpm_common, file = "aligned_TPM_data.rds")
cat("Aligned TPM data saved as 'aligned_TPM_data.rds'\n")

# Save the aligned CNV data
saveRDS(cnv_common, file = "aligned_CNV_data.rds")
cat("Aligned CNV data saved as 'aligned_CNV_data.rds'\n")

# If we want to save as CSV for inspection or other software
write.csv(expr_tpm_common, "aligned_TPM_data.csv", row.names = TRUE)
write.csv(cnv_common, "aligned_CNV_data.csv", row.names = TRUE)
cat("Aligned TPM and CNV data saved as CSV files for inspection.\n")


############# Gene Set Analysis (ssGSEA) #############
# %% Prepare the expression matrix and gene sets
expr_matrix <- as.matrix(expr_tpm_common)

# Load gene sets from GMT file
gmt_file <- "./datafiles/c4.all.v2024.1.Hs.symbols.gmt"
gene_sets <- getGmt(gmt_file)

# Verify the structure of the data
cat("First few row names of expression matrix:\n")
print(head(rownames(expr_matrix)))
cat("First few gene symbols in gene sets:\n")
print(head(geneIds(gene_sets)[[1]]))

# Filter out constant genes if necessary
constant_genes <- apply(expr_matrix, 1, function(x) all(x == x[1]))
expr_matrix <- expr_matrix[!constant_genes, ]
cat("Number of genes after removing constants:", nrow(expr_matrix), "\n")

# %% Set up parameters and run ssGSEA analysis with minSize adjustment
gsvaPar <- gsvaParam(expr_matrix, gene_sets, minSize = 2)
gsva.es <- gsva(gsvaPar, verbose = FALSE)

# %% Save GSVA enrichment scores to a CSV file
write.csv(gsva.es, "gsva_enrichment_scoresC4.csv")
cat("GSVA enrichment scores saved to 'gsva_enrichment_scoresC4.csv'\n")


############# Visualization of Top 5 Enriched Gene Sets #############
# %% Load GSVA enrichment scores
gsva_results <- read.csv("gsva_enrichment_scoresC4.csv", row.names = 1)
cat("GSVA enrichment scores loaded from 'gsva_enrichment_scoresC4.csv'\n")

# Calculate mean enrichment score for each gene set and identify the top 5
gsva_summary <- apply(gsva_results, 1, mean)
top_gene_sets <- sort(gsva_summary, decreasing = TRUE)[1:5]

# Save top 5 enriched gene sets to a CSV file
write.csv(top_gene_sets, "top_5_enriched_gene_sets.csv")
cat("Top 5 enriched gene sets saved to 'top_5_enriched_gene_sets.csv'\n")

# View the top 5 enriched gene sets
top_gene_sets <- head(sort(gsva_summaryC4, decreasing = TRUE), 5)
print(top_gene_sets)


# %% Generate bar plot for top 5 enriched gene sets
png("top_5_gene_sets_barplot.png", width = 800, height = 600)
barplot(top_gene_sets, las = 2, col = "#adaee6",
        main = "Top 5 Enriched Gene Sets",
        ylab = "Mean Enrichment Score", xlab = "Gene Sets",
        names.arg = names(top_gene_sets), cex.names = 0.8,  # Adjust label size
        ylim = c(0, max(top_gene_sets) + 0.1 * max(top_gene_sets)))  # Dynamic y-axis
dev.off()
cat("Bar plot of top 5 enriched gene sets saved to 'top_5_gene_sets_barplot.png'\n")

# Generate heatmap for the top 5 enriched gene sets
top_gene_set_names <- names(top_gene_sets)
png("top_5_gene_sets_heatmap.png", width = 800, height = 600)
pheatmap::pheatmap(gsva_results[top_gene_set_names, ], 
                   show_rownames = TRUE, show_colnames = FALSE,
                   cluster_rows = TRUE, cluster_cols = TRUE,
                   main = "Top 5 Enriched Gene Sets Across Samples",
                   color = colorRampPalette(c("navy", "white", "firebrick3"))(50))  # Color gradient
dev.off()
cat("Heatmap of top 5 enriched gene sets saved to 'top_5_gene_sets_heatmap.png'\n")



############# Filter and Prepare SNP Data for PLINK Analysis #############

# %% Load SNP data
snp_data <- readRDS("TCGA_LUAD_SNP_data.rds")
cat("Original SNP data dimensions:", dim(snp_data), "\n")

# Check row and column names of the SNP data
cat("Row names of SNP data:\n")
print(head(rownames(snp_data)))
cat("Column names of SNP data:\n")
print(head(colnames(snp_data)))

# Print a sample of the data for verification
print(head(snp_data))



# %% Define top 5 gene symbols from enriched gene sets
top_genes <- c("SPINK1", "SERPINI2", "IGFBP1", "EGFR")  # MODULE_215 is a cluster, not an individual gene

# Filter SNP data using the correct column name "Hugo_Symbol"
snp_filtered <- snp_data[snp_data$Hugo_Symbol %in% top_genes, ]

# Check dimensions and sample of the filtered SNP data
cat("Filtered SNP data dimensions:", dim(snp_filtered), "\n")
print(head(snp_filtered))

# %% Save the filtered SNP data for the top genes to a CSV file
write.csv(snp_filtered, "filtered_SNP_data_top_genes.csv", row.names = FALSE)

cat("Filtered SNP data saved to 'filtered_SNP_data_top_genes.csv'\n")



############# Convert to VCF Format for PLINK #############
# %%Define a function to convert MAF data to VCF if necessary
maf_to_vcf <- function(maf_data, output_vcf) {
  # Reformat columns to fit VCF format requirements
  vcf_data <- data.frame(
    CHROM = maf_data$Chromosome, 
    POS = maf_data$Start_Position,
    ID = paste0(maf_data$Hugo_Symbol, "_", maf_data$Variant_Classification),
    REF = maf_data$Reference_Allele, 
    ALT = maf_data$Tumor_Seq_Allele2,
    QUAL = ".", 
    FILTER = "PASS", 
    INFO = paste0("GENE=", maf_data$Hugo_Symbol)
  )
  
  # Write the VCF file with required header for PLINK compatibility
  cat("##fileformat=VCFv4.2\n", file = output_vcf)
  cat("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n", file = output_vcf, append = TRUE)
  write.table(vcf_data, file = output_vcf, append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  cat("Converted MAF data to VCF format and saved to", output_vcf, "\n")
}

# %% Run the MAF to VCF conversion for your SNP data
maf_to_vcf(snp_filtered, "filtered_snp_data.vcf")
cat("Filtered SNP data saved as 'filtered_snp_data.vcf' for PLINK analysis\n")


############# Prepare VCF File for PLINK #############
# %%Confirm VCF format aligns with PLINK requirements
# Sample mapping may be required for genotype-level analyses with PLINK
# For SNP association analysis, this VCF file is ready for PLINK

# Check the first few lines of the VCF file to confirm header and formatting
vcf_file <- "filtered_snp_data.vcf"
vcf_preview <- readLines(vcf_file, n = 10)
cat(vcf_preview, sep = "\n")

# %% Some adjustements (remove the "chr" prefix0
# Load the VCF data without skipping lines
vcf_data <- read.delim("filtered_snp_data.vcf", header = FALSE, sep = "\t", stringsAsFactors = FALSE)

# Identify the line with column names and set it as the header row
header_line <- grep("^#CHROM", vcf_data[,1])
colnames(vcf_data) <- as.character(unlist(vcf_data[header_line, ]))
vcf_data <- vcf_data[(header_line + 1):nrow(vcf_data), ]

# Remove "chr" prefix from the CHROM column if it exists
vcf_data$CHROM <- gsub("^chr", "", vcf_data$CHROM)

# %% Verify successful prefix removal
cat("Unique values in CHROM column after prefix removal:\n")
print(unique(vcf_data$CHROM))

# %% Save the corrected VCF file
output_vcf <- "filtered_snp_data_final.vcf"
cat("##fileformat=VCFv4.2\n", file = output_vcf)
cat("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n", file = output_vcf, append = TRUE)
write.table(vcf_data, file = output_vcf, append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
cat("Final VCF file saved as 'filtered_snp_data_final.vcf'\n")




# %% Stop parallel processing
stopCluster(cl)
