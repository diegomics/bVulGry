# Load required libraries
library(tidyverse)
library(ggplot2)

# Control flags for ontology categories
show_BP <- TRUE   # biological_process
show_CC <- FALSE  # cellular_component
show_MF <- FALSE   # molecular_function

# Function to read and process enrichment files
process_enrichment_file <- function(filepath, category, species) {
  read.table(filepath, sep="\t", quote="", stringsAsFactors=FALSE,
             col.names=c("GO_ID", "gene_count", "description", "ontology", "pvalue")) %>%
    mutate(
      category = category,
      species = species,
      pvalue = as.numeric(pvalue),
      # Cap the -log10(p-value) at 200 for visualization purposes
      neg_log_pvalue = pmin(-log10(pvalue), 200)
    )
}

setwd("/Users/rwaterho/PEHUEL/orthovenn/expansions_contraction_analysis")
# Read the data for both species
speciesA_expanded <- process_enrichment_file("/Users/rwaterho/PEHUEL/orthovenn/expansions_contraction_analysis/4E/enrichment.txt", "Expanded", "AnCo")
speciesA_contracted <- process_enrichment_file("/Users/rwaterho/PEHUEL/orthovenn/expansions_contraction_analysis/4C/enrichment.txt", "Contracted", "AnCo")
speciesB_expanded <- process_enrichment_file("/Users/rwaterho/PEHUEL/orthovenn/expansions_contraction_analysis/5E/enrichment.txt", "Expanded", "CaCo")
speciesB_contracted <- process_enrichment_file("/Users/rwaterho/PEHUEL/orthovenn/expansions_contraction_analysis/5C/enrichment.txt", "Contracted", "CaCo")
# Combine all data
all_data <- bind_rows(
  speciesA_expanded,
  speciesA_contracted,
  speciesB_expanded,
  speciesB_contracted
)

# Filter ontologies based on control flags
selected_ontologies <- c(
  if(show_BP) "biological_process",
  if(show_CC) "cellular_component",
  if(show_MF) "molecular_function"
)

# Filter data by selected ontologies
all_data_filtered <- all_data %>%
  filter(ontology %in% selected_ontologies)

# Select top terms for each ontology based on p-value
top_terms <- all_data_filtered %>%
  group_by(ontology) %>%
  arrange(pvalue) %>%
  slice_head(n = 200) %>%
  pull(GO_ID) %>%
  unique()

# Filter data to include only top terms
plot_data <- all_data_filtered %>%
  filter(GO_ID %in% top_terms)

# Adjust plot height based on number of shown categories
n_categories <- length(selected_ontologies)
plot_height <- 5 * n_categories

# Create the visualization with adjusted color scale
ggplot(plot_data, aes(x = category, y = description)) +
  geom_point(aes(size = gene_count, color = neg_log_pvalue)) +
  scale_color_gradient2(
    low = "lightblue",
    mid = "blue",
    high = "darkblue",
    midpoint = 100,
    limits = c(0, 200),
    oob = scales::squish  # Handle out-of-bounds values
  ) +
  scale_size_continuous(range = c(2, 10)) +
  facet_grid(ontology ~ species, scales = "free_y") +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 8),
    axis.title.y = element_blank(),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 10, face = "bold")
  ) +
  labs(
    x = "",
    color = "-log10(p-value)",
    size = "Gene Count"
  )

# Save the plot with dynamic height
#ggsave("go_enrichment_plot.pdf", width = 15, height = plot_height)

ggsave("go_enrichment_plot.pdf", width = 15, height = plot_height)



library(readr)
write_delim(all_data, 
            file = "/Users/rwaterho/PEHUEL/orthovenn/all_species_data.txt", 
            delim = "\t")
