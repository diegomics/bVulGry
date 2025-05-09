# Chromosome-Scale Genome Assembly Provides Insights into Condor Evolution and Conservation. De Panis et al. 2025.
# Fig. S8

library(ggplot2)
library(dplyr)

# Function to read and process bootstrap data
read_psmc_data <- function(file_paths, individual_name) {
  # Read all bootstrap files
  data_list <- lapply(file_paths, function(file) {
    # Read the .txt file into a data frame
    df <- read.table(file, header = FALSE)
    
    # Assign column names to match the data structure
    colnames(df) <- c("Years", "EffectivePopulationSize", "NA1", "NA2", "NA3")
    
    # Filter out invalid rows (e.g., rows with Years <= 0 or EffectivePopulationSize <= 0)
    df <- df %>% filter(Years > 0, EffectivePopulationSize > 0)
    
    # Add metadata: individual name and bootstrap file ID
    df$Individual <- individual_name
    df$BootstrapID <- basename(file)  # Use the filename as a unique ID
    
    return(df)  # Return processed data frame
  })
  
  # Combine all bootstrap data into a single data frame
  return(do.call(rbind, data_list))
}

# Function to read and process consensus data
read_consensus_data <- function(file_path, individual_name) {
  # Read the consensus file into a data frame
  df <- read.table(file_path, header = FALSE)
  
  # Assign column names to match the data structure
  colnames(df) <- c("Years", "EffectivePopulationSize", "NA1", "NA2", "NA3")
  
  # Filter out invalid rows (e.g., rows with Years <= 0 or EffectivePopulationSize <= 0)
  df <- df %>% filter(Years > 0, EffectivePopulationSize > 0)
  
  # Add metadata: individual name
  df$Individual <- individual_name
  
  return(df)  # Return processed data frame
}

# Define individuals and their corresponding data paths
# Each individual has a name, a directory for bootstrap files, and a consensus file
individuals <- list(
  list(name = "AnCo1", path = "../psmc", pattern = "Pehuel_plot\\.\\d+\\.txt$", consensus = "../PEHUEL/psmc/Pehuel_plot.0.txt"),
  list(name = "AnCo2", path = "../psmc", pattern = "SAMN18477436_plot\\.\\d+\\.txt$", consensus = "../psmc/SAMN18477436_plot.0.txt"),
  list(name = "CaCo1", path = "../psmc", pattern = "SAMN18477434_plot\\.\\d+\\.txt$", consensus = "../psmc/SAMN18477434_plot.0.txt"),
  list(name = "CaCo2", path = "../psmc", pattern = "SAMN18477435_plot\\.\\d+\\.txt$", consensus = "../psmc/SAMN18477435_plot.0.txt")
)

# Specify which individuals to include in the plot
# Modify this vector to plot different individuals
selected_individuals <- c("AnCo1","AnCo2","CaCo1","CaCo2")

# Initialize lists to store bootstrap and consensus data for selected individuals
all_bootstrap_data <- list()
all_consensus_data <- list()

# Loop through all individuals and read their data if they are selected
for (ind in individuals) {
  if (ind$name %in% selected_individuals) {
    # Get the list of bootstrap files for the individual
    bootstrap_files <- list.files(ind$path, full.names = TRUE, pattern = ind$pattern)
    
    # Get the consensus file for the individual
    consensus_file <- ind$consensus
    
    # Read bootstrap and consensus data
    bootstrap_data <- read_psmc_data(bootstrap_files, ind$name)
    consensus_data <- read_consensus_data(consensus_file, ind$name)
    
    # Store the data in the lists
    all_bootstrap_data[[ind$name]] <- bootstrap_data
    all_consensus_data[[ind$name]] <- consensus_data
  }
}

# Combine all bootstrap and consensus data across individuals into single data frames
all_bootstrap_data <- do.call(rbind, all_bootstrap_data)
all_consensus_data <- do.call(rbind, all_consensus_data)

# Print summaries of EffectivePopulationSize for verification
print(summary(all_bootstrap_data$EffectivePopulationSize))
print(summary(all_consensus_data$EffectivePopulationSize))

# Define axis limits for the plot
x_limits <- c(1e4, 1e7)  # X-axis (Years)
y_limits <- NULL          # Let ggplot automatically adjust Y-axis based on data

# Create the plot
ggplot() +
  # Plot bootstrap replicates as thin, semi-transparent lines
  geom_step(data = all_bootstrap_data, aes(x = Years, y = EffectivePopulationSize, group = BootstrapID, color = Individual), 
            alpha = 0.2, size = 0.5) +
  # Plot consensus as a thick, opaque line
  geom_step(data = all_consensus_data, aes(x = Years, y = EffectivePopulationSize, color = Individual), size = 1.5) +
  # Add labels and title
  labs(
    title = "",
    x = expression(paste("Years (log scale, g=10, ", mu, "=1.4", "\u00D7", "10"^-8, ")")),
    y = expression(paste("Effective Population Size (", 10^4, ")")),
    color = "Individual"
  ) +
  # Logarithmic scale for X-axis and dynamic Y-axis
  scale_x_log10(limits = x_limits, 
                labels = scales::trans_format("log10", 
                                              function(x) parse(text = paste0("10^", x)))) +
  scale_y_continuous(limits = y_limits) +
  scale_color_manual(values = c(
    "AnCo1" = "#3584e4",
    "AnCo2" = "#26a269",
    "CaCo1" = "#e01b24",
    "CaCo2" = "#ff7800"
  )) +
  theme_minimal() +
  theme(
    text = element_text(size = 20),
    legend.position = ""
  )
