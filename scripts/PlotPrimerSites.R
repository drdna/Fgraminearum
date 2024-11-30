#R

library(ggplot2)
library(scales)

# Read the file into R
data <- read.table("Picked_primers.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)

# Assign column names
colnames(data) <- c("Chromosome", "Position", "Seq1", "Seq2")

# Create a new data frame with just the Chromosome and Position columns
plot_data <- data.frame(
  Chromosome = data$Chromosome,
  Position = data$Position
)

# You can add a column to represent the sequence number or ID if desired
plot_data$Sequence_ID <- seq_along(plot_data$Position)


# Plot the positions using facet_grid
ggplot(plot_data, aes(x = Position, y = 0)) +
  geom_point(alpha = 0.5, color = "blue") +  # Plot the positions as dots
  facet_grid(Chromosome ~ ., scales = "free_x") +  # Facet by Chromosome (rows)
  labs(title = "Positions of MonsterPlex Target Loci Along Chromosomes",
       subtitle= "(high frequency variant sites f(A)<0.9, flanked by invariant primer sites)",
       x = "Position on Chromosome (Mb)") +
  scale_x_continuous(
    breaks = seq(0, max(plot_data$Position), by = 1e6),   # Set breaks every 1 million
    labels = scales::comma  # Format labels with commas (optional)
  ) +
  theme(axis.text.y = element_blank(),  # Remove y-axis labels (Sequence IDs)
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank())  # Remove y-axis ticks
