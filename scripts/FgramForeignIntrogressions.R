# Load necessary libraries
library(tidyverse)

# Set the directory where the files are located
dir_path <- "~/FgramCopyProbs"

# Get the list of files in the directory
file_list <- list.files(dir_path, full.names = TRUE)

# Define the donor populations (including NA1, NA2, NA3) with specific order
donor_order <- c("FMER", "NA2", "FCOR", "FASI", "Fl", "Fb", "FAUS", "Fger", "NA3", "NA1")

# Define donor colors
donor_colors <- c("FMER" = "red", "NA2" = "gray", "FCOR" = "blue", "FASI" = "green", 
                  "Fl" = "purple", "Fb" = "orange", "FAUS" = "cyan", "Fger" = "magenta", 
                  "NA3" = "gray", "NA1" = "gray")

# Initialize an empty list to store data frames
all_data <- list()

# Read and process each file
for (file in file_list) {
  # Extract strain ID from filename (assuming the strain ID is before the first period)
  strain_id <- str_extract(basename(file), "^[^.]+")
  
  # Read the data
  df <- read.table(file, header = TRUE, stringsAsFactors = FALSE)
  
  # Add strain_id to the data
  df$strain_id <- strain_id
  
  # Append to the list of data frames
  all_data[[strain_id]] <- df
}

# Combine all the data frames into one
combined_data <- bind_rows(all_data)

# Function to assign donor names (not colors) based on conditions
assign_color <- function(row, donor_order, donor_colors) {
  # Extract donor probabilities for each row
  donor_probs <- row[donor_order]
  
  # 1. Check if any donor has probability > 0.9, and assign donor name
  donor_name <- NA
  for (donor in donor_order) {
    if (donor_probs[donor] > 0.9) {
      donor_name <- donor
      break  # Exit the loop once we find a donor with probability > 0.9
    }
  }
  
  # 2. If no donor with probability > 0.9, check if NA1, NA2, NA3 have probabilities between 0.1 and 0.4
  if (is.na(donor_name)) {
    na_probs <- row[c("NA1", "NA2", "NA3")]
    if (all(na_probs > 0.1 & na_probs < 0.4)) {
      donor_name <- "gray"  # Assign "gray" if NA1, NA2, NA3 meet the condition
    }
  }
  
  # Return the assigned donor name (or "gray" or NA if no color is assigned)
  return(donor_name)
}

# Apply the function row-wise to assign donor names
combined_data$color_assigned <- pmap_chr(combined_data, function(...) {
  row <- list(...)
  assign_color(row, donor_order, donor_colors)
})

# Check the first few rows of the data with color_assigned
print(combined_data %>% select(pos, strain_id, color_assigned) %>% head(20))


# Apply the function row-wise to assign colors using pmap
combined_data$color_assigned <- pmap_chr(combined_data, function(...) {
  row <- list(...)
  assign_color(row, donor_order, donor_colors)
})

# Ensure the color_assigned column is a character vector
combined_data$color_assigned <- as.character(combined_data$color_assigned)

# Check the color_assigned values before filtering
print(combined_data %>% select(pos, strain_id, color_assigned) %>% head(20))

# Keep only rows with valid colors (non-NA, non-gray)
combined_data_filtered <- combined_data %>%
  filter(!is.na(color_assigned)) # & color_assigned != "gray")

# Print the filtered data for verification (including strain_id)
print(combined_data_filtered)

combined_data_filtered$donor <- combined_data_filtered$color_assigned

# Create a new column 'point_size' based on donor values
combined_data_filtered <- combined_data_filtered %>%
  mutate(point_size = case_when(
    donor %in% c("NA1", "NA2", "NA3") ~ 2,  # Set size 2 for NA1, NA2, NA3
    TRUE ~ 10  # Set size 10 for all other donors
  ))

# Check how many rows are left after filtering
cat("Number of rows after filtering: ", nrow(combined_data_filtered), "\n")



# Remove gray circles if needed
#combined_data_filtered <- combined_data_filtered[!combined_data_filtered$color_assigned %in% c("NA1", "NA2", "NA3"),]

# Create the plot with faceting by strain_id
# Create the plot with faceting by strain_id
pdf("~/FgramCPplots.pdf", 60, 24)
ggplot(combined_data_filtered, aes(x = pos, y = 0, color = color_assigned)) +
  geom_point(aes(size = point_size), shape = 16, stroke = 0.5) +  # Circle plot
  xlim(6300000, 7500000) +
  coord_cartesian(clip = "off") +
  facet_grid(strain_id ~ .) +  # Facet by strain
  scale_color_manual(values = donor_colors) +  # Custom color mapping based on donor names
  theme_minimal() +  # Minimal theme
  labs(x = "Chromosome Position", y = NULL, color = "Donor Population", size = "Probability") +
  theme(legend.position = "bottom",
        axis.text.y = element_blank(),  # Hide y-axis labels (strain names are in facets)
        axis.ticks.y = element_blank(),
        strip.placement = "right",
        strip.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
        panel.spacing = unit(0.5, "lines"),
        plot.margin = unit(c(5, 5, 5, 15), "mm"))
dev.off()
