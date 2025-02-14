# Load required libraries
library(tidyverse)
library(ggplot2)

# Control flags for ontology categories
show_BP <- TRUE   # biological_process
show_CC <- FALSE  # cellular_component
show_MF <- TRUE   # molecular_function

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



##### Fig. 3

# Read the data for both species
speciesA_expanded <- process_enrichment_file("../exp_con/gry_E/enrichment.txt", "Expanded", "AnCo")
speciesA_contracted <- process_enrichment_file("../exp_con/gry_C/enrichment.txt", "Contracted", "AnCo")
speciesB_expanded <- process_enrichment_file("../exp_con/cal_E/enrichment.txt", "Expanded", "CaCo")
speciesB_contracted <- process_enrichment_file("../exp_con/cal_C/enrichment.txt", "Contracted", "CaCo")

# Combine all data
all_data <- bind_rows(
  speciesA_expanded,
  speciesA_contracted,
  speciesB_expanded,
  speciesB_contracted
)

dim(all_data)

# Filter ontologies based on control flags to get rid of CC at the moment
selected_ontologies <- c(
  if(show_BP) "biological_process",
  if(show_CC) "cellular_component",
  if(show_MF) "molecular_function"
)

# Filter data by selected ontologies
all_data_BPMF <- all_data %>%
  filter(ontology %in% selected_ontologies)

dim(all_data_BPMF)

# Interesting GOs for the plot
go_list <- read.table("../exp_con/GO_CONDORS", header=TRUE)

dim(go_list)

# Filter and prepare the data
filtered_data <- all_data_BPMF %>%
  filter(GO_ID %in% go_list$GO_ID) %>%
  mutate(
    # Create a factor with levels in reverse order of go_list
    GO_ID = factor(GO_ID, levels = rev(go_list$GO_ID))
  )

# Select top terms for each ontology based on p-value
top_terms <- filtered_data %>%
  #arrange(pvalue) %>%  # Remove group_by(ontology)
  #arrange(desc(gene_count)) %>%
  distinct(GO_ID) %>%  # Take unique GO_IDs
  slice_head(n = 100) %>%  # Take top
  pull(GO_ID) %>%
  unique()

data.frame(GO_ID = top_terms) %>% print()

# Filter data to include only top terms
plot_data <- filtered_data %>%
  filter(GO_ID %in% top_terms)


# Create the plot
ggplot(plot_data, aes(x = species, y = GO_ID)) +  # Using GO_ID instead of description
  geom_point(aes(size = gene_count, color = neg_log_pvalue)) +
  scale_color_gradient2(
    low = "#fdeaea",
    mid = "#fc4c4c",
    high = "#f40000",
#    low = "#e0e0fc",
#    mid = "#4c4cfc",
#    high = "#0000f4",
    midpoint = 100,
    limits = c(0, 200),
    oob = scales::squish
  ) +
  scale_size_continuous(range = c(2, 10), 
                        breaks = c(10, 25, 50, 100)) +
  facet_grid(. ~ category) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 9),
    axis.title.y = element_blank(),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 12, face = "bold"),
    plot.margin = margin(1, 1, 1, 1, "cm")
  ) +
  labs(
    x = "",
    color = "-log10(p-value)",
    size = "Gene Count"
  )



##### Fig. S4

speciesA_expanded <- process_enrichment_file("../exp_con/nw_E/enrichment.txt", "Expanded", "NW")
speciesA_contracted <- process_enrichment_file("../exp_con/nw_C/enrichment.txt", "Contracted", "NW")
speciesB_expanded <- process_enrichment_file("../exp_con/ow_E/enrichment.txt", "Expanded", "OW")
speciesB_contracted <- process_enrichment_file("../exp_con/ow_C/enrichment.txt", "Contracted", "OW")

# Combine all data
all_data <- bind_rows(
  speciesA_expanded,
  speciesA_contracted,
  speciesB_expanded,
  speciesB_contracted
)

dim(all_data)

# Filter ontologies based on control flags to get rid of CC at the moment
selected_ontologies <- c(
  if(show_BP) "biological_process",
  if(show_CC) "cellular_component",
  if(show_MF) "molecular_function"
)

# Filter data by selected ontologies
all_data_BPMF <- all_data %>%
  filter(ontology %in% selected_ontologies)

dim(all_data_BPMF)

# Select top terms for each ontology based on p-value
top_terms <- all_data_BPMF %>%
  group_by(ontology) %>%
  arrange(pvalue) %>%
  slice_head(n = 25) %>%
  pull(GO_ID) %>%
  unique()

# Filter data to include only top terms
plot_data <- all_data_BPMF %>%
  filter(GO_ID %in% top_terms)

# Adjust plot height based on number of shown categories
n_categories <- length(selected_ontologies)
plot_height <- 5 * n_categories


# Interesting GOs for the plot
go_list <- read.table("../exp_con/GO_NWOW", header=TRUE)

dim(go_list)

# Filter and prepare the data
filtered_data <- all_data_BPMF %>%
  filter(GO_ID %in% go_list$GO_ID) %>%
  mutate(
    # Create a factor with levels in reverse order of go_list
    GO_ID = factor(GO_ID, levels = rev(go_list$GO_ID))
  )

# Select top terms for each ontology based on p-value
top_terms <- filtered_data %>%
  #arrange(pvalue) %>%  # Remove group_by(ontology)
  #arrange(desc(gene_count)) %>%
  distinct(GO_ID) %>%  # Take unique GO_IDs
  slice_head(n = 100) %>%  # Take top 
  pull(GO_ID) %>%
  unique()

data.frame(GO_ID = top_terms) %>% print()

# Filter data to include only top terms
plot_data <- filtered_data %>%
  filter(GO_ID %in% top_terms)


# Create the plot
ggplot(plot_data, aes(x = species, y = GO_ID)) +  # Using GO_ID instead of description
  geom_point(aes(size = gene_count, color = neg_log_pvalue)) +
  scale_color_gradient2(
    low = "#fdeaea",
    mid = "#fc4c4c",
    high = "#f40000",
#        low = "#e0e0fc",
#        mid = "#4c4cfc",
#        high = "#0000f4",
    midpoint = 100,
    limits = c(0, 200),
    oob = scales::squish
  ) +
  scale_size_continuous(range = c(2, 10), 
                        breaks = c(10, 25, 50, 100)) +
  facet_grid(. ~ category) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 9),
    axis.title.y = element_blank(),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 12, face = "bold"),
    plot.margin = margin(1, 1, 1, 1, "cm")
  ) +
  labs(
    x = "",
    color = "-log10(p-value)",
    size = "Gene Count"
  )
