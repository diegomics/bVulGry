# Chromosome-Scale Genome Assembly Provides Insights into Condor Evolution and Conservation. De Panis et al. 2025.
# Fig. S7

library(ggplot2)
library(dplyr)
library(tidyr)

##### A: AnCo
pi_data <- read.delim("AnCo_diversity_analysis.windowed.pi", header=TRUE)

# Define chromosome list to show (NULL to include all)
chrom_list <- c("SUPER_1", "SUPER_2", "SUPER_3", "SUPER_4", "SUPER_5",
                "SUPER_6", "SUPER_7", "SUPER_8", "SUPER_9", "SUPER_10",
                "SUPER_11", "SUPER_12", "SUPER_13", "SUPER_14")
#chrom_list <- NULL

show_sd <- FALSE

# Filter data with chrom_list
if (!is.null(chrom_list)) {
  pi_data_filtered <- pi_data %>% filter(CHROM %in% chrom_list)
}

# Calculate average π
genome_wide_pi <- mean(pi_data$PI)
cat("Genome-wide average π:", format(genome_wide_pi, scientific=TRUE), "\n")
genome_wide_pi_filtered <- mean(pi_data_filtered$PI)
cat("Genome-wide average π (filtered):", format(genome_wide_pi_filtered, scientific=TRUE), "\n")

# Calculate chromosome-specific averages
chrom_averages <- pi_data_filtered %>%
  group_by(CHROM) %>%
  summarise(
    mean_pi = mean(PI),
    sd_pi = sd(PI),
    windows = n()
  ) %>%
  mutate(CHROM_numeric = as.numeric(gsub("SUPER_", "", CHROM))) %>%
  arrange(CHROM_numeric) %>%
  mutate(CHROM = factor(CHROM, levels = CHROM)) %>%  # Ensure ggplot respects order
  select(-CHROM_numeric)  # Drop helper column

# Print chromosome averages
print("Chromosome averages:")
print(chrom_averages, n=50)

# Create visualization
ggplot(chrom_averages, aes(x=CHROM, y=mean_pi)) +
  geom_bar(stat="identity", fill="#7B9EE8") +
  geom_hline(yintercept = genome_wide_pi, linetype="dashed", color="#960018", linewidth=1) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 11, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 11),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold")
  ) +
  labs(
    title = "",
    x = "Chromosome",
    y = "Mean π"
  )


##### B: CaCo 
pi_data <- read.delim("CaCo_diversity_analysis.windowed.pi", header=TRUE)

# Define chromosome list to show (NULL to include all)
chrom_list <- c("NC_059471.1", "NC_059472.1", "NC_059473.1", "NC_059474.1", "NC_059475.1",
                "NC_059476.1", "NC_059477.1", "NC_059478.1", "NC_059479.1", "NC_059480.1",
                "NC_059481.1", "NC_059482.1", "NC_059483.1", "NC_059484.1")
#caco_chrom_list <- NULL

show_sd <- FALSE

# Filter data with chrom_list
if (!is.null(chrom_list)) {
  pi_data_filtered <- pi_data %>% filter(CHROM %in% chrom_list)
}

# Calculate average π
genome_wide_pi <- mean(pi_data$PI)
cat("Genome-wide average π:", format(genome_wide_pi, scientific=TRUE), "\n")
genome_wide_pi_filtered <- mean(pi_data_filtered$PI)
cat("Genome-wide average π (filtered):", format(genome_wide_pi_filtered, scientific=TRUE), "\n")

# Calculate chromosome-specific averages
chrom_averages <- pi_data_filtered %>%
  group_by(CHROM) %>%
  summarise(
    mean_pi = mean(PI),
    sd_pi = sd(PI),
    windows = n()
  ) %>%
  # Extract the numeric part from the CHROM name
  # Remove the "NC_" prefix and any leading zeros, captures the integer part
  mutate(CHROM_numeric = as.numeric(sub("NC_0*([0-9]+).*", "\\1", CHROM))) %>%
  arrange(CHROM_numeric) %>%
  # Convert CHROM to a factor with levels in the order of appearance (now sorted)
  mutate(CHROM = factor(CHROM, levels = CHROM)) %>%
  select(-CHROM_numeric)  # Remove the helper column


# Print chromosome averages
print("Chromosome averages:")
print(chrom_averages, n=50)

# Create visualization
ggplot(chrom_averages, aes(x=CHROM, y=mean_pi)) +
  geom_bar(stat="identity", fill="#FF8C69") +
  geom_hline(yintercept = genome_wide_pi, linetype="dashed", color="#960018", linewidth=1) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 11, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 11),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold")
  ) +
  labs(
    title = "",
    x = "Chromosome",
    y = "Mean π"
  )
