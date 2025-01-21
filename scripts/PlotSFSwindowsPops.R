# Author: Mostafa Rahnama
# Plotting copyprobsperlocus.out data from CP based on chromosomes and strains of interest 

library(reshape2)
library(scales)
library(ggplot2)
library(grid)
library(dplyr)
library(RColorBrewer)
library(gridExtra)
library(stringr)
library(readxl)
library(readr)
library(tidyr)

'%notin%' <- Negate('%in%')

# Specify the file path
haplotypesFile <- "~/FgWardPlusHaplotypes.complete.txt"

# Chromosome to analyze
Chr <- "sequence3"

# Read haplotype data
dt <- as.data.frame(read.table(haplotypesFile, header=FALSE, skip=8, as.is=TRUE, colClasses="character"))
colnames(dt) <- c("id", "pop", "sequence", "haplotype")


# Now extract the haplotype sequences from the filtered dataset for the given chromosome
dt_filtered <- dt[dt$sequence == Chr, ]

# Filter the dataset to keep only NA1, NA2, and NA3 populations

population = c("NA1")#, "NA2", "NA3")
titlePop <- paste(population, collapse = ", ")

dt_filtered <- dt_filtered %>%
  filter(pop %in% c(population))

# Ensure that we only keep the relevant columns ('haplotype' column contains the SNP data)
haplotype_data <- dt_filtered$haplotype

print(paste("Length of the first haplotype:", nchar(haplotype_data[1])))  # Length of the first haplotype
print(paste("Number of haplotypes:", length(haplotype_data)))  # Number of haplotypes

# The number of sites is determined by the length of the 'haplotype' string (each SNP corresponds to a position)
filtered_length <- nchar(haplotype_data[1])  # Assuming all haplotypes have the same length for the filtered sequence

# Print to ensure filtered_length is correct
print(paste("Filtered length:", filtered_length))

# Ensure filtered_length is valid
if (is.na(filtered_length) || filtered_length <= 0) {
  stop("Invalid filtered length. Please check your input data.")
}

# Define 20 equal-sized windows for the filtered SNP sites
window_size <- floor(filtered_length / 20)  # Ensure windows are of equal size

# Print to check window_size calculation
print(paste("Window size:", window_size))

# Ensure window_size is valid
if (window_size <= 0) {
  stop("Window size is too small. Please check your input data.")
}

num_blocks <- 20  # Define 20 blocks

# Create an empty list to store the SFS for each block
block_sfs <- list()

# Loop through the dataset, creating 20 equal-sized blocks of SNPs
for (block in 1:num_blocks) {
  
  # Define the start and end SNP positions for this block (based on filtered data)
  block_start <- (block - 1) * window_size + 1
  block_end <- block * window_size
  
  # Adjust the block_end for the last block to capture any remaining SNPs
  if (block == num_blocks) {
    block_end <- filtered_length  # Ensure the last block includes the remainder of the SNPs
  }
  
  # Check if block_end is within range, prevent NA or negative values
  if (block_start > filtered_length | block_end > filtered_length) {
    next  # Skip this block if the range is out of bounds
  }
  
  # Extract the SNP positions for this block from the haplotypes (columns within the window)
  block_haplotypes <- substr(haplotype_data, block_start, block_end)
  
  # Split the haplotypes into columns (one per individual)
  haplotype_split <- str_split(block_haplotypes, "")
  hdf_block <- as.data.frame(do.call(cbind, haplotype_split))
  colnames(hdf_block) <- dt_filtered$id
  rownames(hdf_block) <- block_start:block_end  # SNP positions for this block
  
  # Convert to numeric genotypes (0/1 format)
  df_block <- as.data.frame(t(hdf_block))
  df_block <- df_block %>%
    mutate(across(everything(), ~ as.numeric(as.character(.))))
  
  # Convert to long format for this block
  long_data_block <- df_block %>%
    pivot_longer(cols = everything(), names_to = "snp", values_to = "genotype") %>%
    mutate(snp = as.numeric(as.character(snp)))  # Ensure SNP positions are treated as numbers
  
  # Calculate allele frequencies for this block
  allele_freqs_block <- long_data_block %>%
    group_by(snp) %>%
    summarize(allele_freq = mean(genotype), .groups = "drop")
  
  # Calculate the folded SFS: min(allele_freq, 1 - allele_freq)
  allele_freqs_block <- allele_freqs_block %>%
    mutate(folded_sfs = pmin(allele_freq, 1 - allele_freq))
  
  # Add a column for block number to facilitate faceting
  allele_freqs_block$block <- factor(block, levels = 1:num_blocks)
  
  # Store the folded SFS for this block
  block_sfs[[block]] <- allele_freqs_block
}

# Combine the results from all blocks into a single data frame
sfs_all_blocks <- do.call(rbind, block_sfs)

# Filter out rows where folded_sfs is 0 before plotting
sfs_all_blocks <- sfs_all_blocks %>%
  filter(folded_sfs > 0)  # Remove rows where folded_sfs is 0

# Plot the folded SFS for all blocks, faceted by block number
p <- ggplot(sfs_all_blocks, aes(x = folded_sfs)) +
  geom_histogram(binwidth = 0.025, fill = "blue", color = "black", alpha = 0.7) +
  facet_wrap(~ block, ncol = 5, scales = "free") +  # Facet by block
  labs(title = paste0("Folded Site Frequency Spectrum Across ", num_blocks, " Equal-Sized Blocks of ", window_size, " SNPs (After Filtering): Chromosome ", Chr, "; population, ", titlePop),
       x = "Folded Allele Frequency",
       y = "Count") +
  theme_minimal() +
  theme(strip.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1)) 

print(p)

pdf(paste0("SFS_", Chr, "_", titlePop, ".pdf"), 11, 8.5)
print(p)
dev.off()
