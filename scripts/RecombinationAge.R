# Load necessary libraries
library(ggplot2)
library(dplyr)

# Path to your data file (update this path to where your file is located)
file_path <- "~/Downloads/FgramARGchr3_1-1000000.txt_145118-999584.143.smc.arg"  # Replace with your actual file path

# Read the data file
data <- read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE, skip = 1)


# Filter the data to include only rows where the 'event' column contains "recomb"
recomb_data <- data[grep("recomb", data$event), ]

# Create a new column for the 100 KB bin based on 'pos' column
recomb_data <- recomb_data %>%
  mutate(pos_bin = floor(pos / 100000) * 100000)  # Group by 100 KB intervals

# Create a new column to bin the 'age' values into specific ranges (e.g., every 1000 units)
recomb_data <- recomb_data %>%
  mutate(age_bin = cut(age, breaks = seq(0, max(age), by = 10000), 
                       right = FALSE, include.lowest = TRUE))  # Adjusted bin width to 1000

# Now we count the number of recombination events in each 'age_bin' for each 'pos_bin'
recomb_counts <- recomb_data %>%
  group_by(pos_bin, age_bin) %>%
  summarise(count = n(), .groups = "drop")  # Count the number of events in each bin

# Custom function to format facet labels as "100 KB", "200 KB", etc.
recomb_counts <- recomb_counts %>%
  mutate(pos_bin_label = paste0(format(pos_bin / 1000, big.mark = ","), " KB"))

# Plotting the distribution of recombination events by age (timing) within 100 KB bins
ggplot(recomb_counts, aes(x = age_bin, y = count)) +
  geom_bar(stat = "identity", fill = "skyblue", color = "black", alpha = 0.7) +  # Bar plot for count of events
  facet_wrap(~ pos_bin_label, scales = "free_y") +  # Facet by 100 KB bins with formatted labels
  labs(title = "Recombination Event Count by Timing (Age) and Chromosomal Position (100 KB bins)",
       x = "Recombination Event Timing (Age bins)",
       y = "Count of Recombination Events") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
        axis.text.y = element_text(size = 10),               # Adjust y-axis label size
        strip.text = element_text(size = 10)) +              # Adjust facet labels size
  scale_x_discrete(drop = FALSE) +  # Don't drop any age bins
  scale_y_continuous(labels = scales::comma)  # Use comma for y-axis labels
