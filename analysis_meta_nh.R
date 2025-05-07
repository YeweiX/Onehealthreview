# --- Install and Load Required Packages ---
# Ensure you have these packages installed. Uncomment the line below and run if needed.
# install.packages(c("meta", "dplyr", "tidyr", "readxl", "openxlsx",
#                    "ggplot2", "viridisLite", "viridis", "reshape2",
#                    "scales", "officer", "rvg", "cowplot", "here", "metafor", "broom", "knitr", "ggsci", "ggthemes",
#                    "circlize", "RColorBrewer", "grid")) # Added packages for gene analysis/chord plot

# Load core data manipulation and meta-analysis packages
library(meta)     # For metaprop function (meta-analysis of proportions)
library(dplyr)    # For data manipulation (%>% pipe, filter, mutate, group_by, etc.)
library(tidyr)    # For data tidying (not explicitly used in provided code, but often useful)
library(readxl)   # For reading Excel files (.xlsx)
library(openxlsx) # For writing Excel files (.xlsx)
library(metafor)  # For meta-analysis and meta-regression (used in time-trend)
library(broom)    # For tidying model outputs (used with lm in time-trend)
library(knitr)    # For creating tables (used for trend_p table)

# Load plotting packages
library(ggplot2)  # For creating plots
library(viridisLite) # For viridis color palettes
library(viridis)  # For viridis color scales
library(reshape2) # For data reshaping (not explicitly used, but often associated with heatmaps)
library(scales)   # For formatting scales (e.g., percent_format)
library(ggsci)    # For scientific journal color palettes (used in time-trend plot)
library(ggthemes) # For additional plot themes (used in time-trend plot)
library(cowplot)  # For combining plots (used for country heatmaps)

# Load packages for saving plots to PowerPoint
library(officer)  # For creating and modifying PowerPoint documents
library(rvg)      # For exporting ggplot objects as vector graphics (dml)

# Load packages for gene analysis visualization (Chord Plot)
library(circlize)    # For creating chord diagrams
library(RColorBrewer) # For color palettes
library(grid)        # For grid graphics (used by circlize)


# Optional: Use 'here' package for reproducible file paths across different systems
# library(here)

# --- Set Data File Paths ---
# Define the path to your main Excel data file.
# Assuming 'Supplementary_file_analysis.xlsx' is in a 'data' subfolder in your project directory.
data_file_path <- "./data/Supplementary_file_analysis.xlsx"

data_file_path <- "C:/Users/e0977434/Desktop/OneDrive/PhD thesis/Prevalence review/Data analysis/Submission/CMI/Supplementary_file_yx_check.xlsx" # Update this path to your actual data file location
# Define the sheet names for different analyses from the main data file
overall_country_sheet <- "Table S1" # Sheet for overall and country meta-analysis
time_trend_sheet <- "Table S4"      # Sheet for time-trend analysis
gene_data_sheet <- "Table S7"       # Sheet for gene analysis


# --- Analysis Section 1: Time-Trend Analysis ---

cat("Performing time-trend analysis...\n")

# Load data specifically for time-trend analysis from the specified file/sheet
# Data columns expected: E, N, Class, Ctime, Sector ("Ctime" = "Midpoint of the year", rawp"= "Proportion", "E"=	"Resistant isolates (n)" and "N"= "Isolates (N)")
data_timetrend <- read_excel(data_file_path, sheet = time_trend_sheet)

# Filter data for specific drug classes for the trend analysis
filtered_data_timetrend <- data_timetrend %>%
  filter(Class %in% c("3GC", "Carbapenems", "Fluoroquinolones"))

# Reorder the 'Sector' factor for consistent plotting
filtered_data_timetrend$Sector <- factor(filtered_data_timetrend$Sector, levels = c("Human", "Animal", "Environment", "Food"))

# Define a function to calculate trend p-values using weighted linear regression (lm)
calculate_trend_p <- function(data) {
  # Perform linear regression of raw proportion (rawp) on Year (Ctime), weighted by sample size (N)
  # Ensure 'rawp' and 'Ctime' columns exist in the input 'data' dataframe

  # Assuming Table S4 contains pre-calculated 'rawp', 'Ctime', and 'N' for this analysis.
  model <- lm(rawp ~ Ctime, data = data, weights = N)
  summary_model <- summary(model)
  
  # Return key statistics from the model summary as a tibble
  return(tibble(
    Intercept = summary_model$coefficients[1, 1],
    Slope = summary_model$coefficients[2, 1],
    p_value = summary_model$coefficients[2, 4], # p-value for the trend (Ctime coefficient)
    R_squared = summary_model$r.squared,
    Adj_R_squared = summary_model$adj.r.squared,
    Residual_SE = summary_model$sigma
  ))
}

# Group the filtered data by Drug Class and Sector, then apply the trend calculation function to each group
# Note: This grouping does NOT include Country, so the trend analysis is aggregated across countries within each Class and Sector.
trend_p_results <- filtered_data_timetrend %>%
  group_by(Class, Sector) %>%
  group_modify(~ calculate_trend_p(.x))

# Display the trend p-values in a kable table (optional, good for R Markdown)
library(knitr) # Already loaded at the top
kable(trend_p_results)

# --- Figures - Time-Trend Plots ---

cat("Generating time-trend plots...\n")

# Create the time-trend plot using ggplot2
# This plot facets by Class and Sector, NOT by Country.
timetrend_plot <- ggplot(filtered_data_timetrend, aes(x = Ctime, y = rawp)) +
  # Add LOESS smooth curve (local regression)
  geom_smooth(aes(weight = N), method = "loess", se = FALSE,
              linewidth = 0.5, span = 1, color = "darkblue") +
  # Add linear model (LM) smooth line (weighted by N)
  geom_smooth(aes(weight = N), method = "lm", se = FALSE, linewidth = 0.8,
              formula = y ~ x, color = "#E31A1C", linetype = "dashed") + # Using a distinct color
  # Add data points, mapping sample size to alpha transparency
  geom_point(aes(alpha = N), color = "grey") +
  # Add text labels for trend p-values in the top-right corner of each facet
  geom_text(data = trend_p_results, # Use the trend_p_results data frame
            aes(x = Inf, y = Inf, # Position text at infinity (top-right corner)
                label = ifelse(p_value < 0.001, "p<0.001", sprintf("p=%.3f", p_value))), # Format p-value
            hjust = 1.1, vjust = 1.5, size = 3.5, color = "#4E79A7", # Adjust position, size, color
            check_overlap = TRUE) + # Prevent text overlap
  # Set alpha scale range for points
  scale_alpha_continuous(range = c(0.2, 1)) +
  # Set x-axis breaks
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
  # Set y-axis limits and format as percentage
  scale_y_continuous(limits = c(0, 1), labels = scales::percent,
                     expand = expansion(mult = c(0, 0.1))) + # Expand y-axis slightly
  # Set plot labels and title
  labs(
    x = "Year",
    y = "Prevalence (%)",
    alpha = "Sample size", # Label for the alpha legend
    title = "Temporal trend of prevalence by drug class and sector" # More descriptive title
  ) +
  # Apply a minimal theme with Times New Roman base font (ensure font is available)
  theme_minimal(base_family = "Times") +
  # Customize theme elements
  theme(
    plot.background = element_rect(fill = "white", color = NA), # Plot background
    panel.background = element_rect(fill = "white", color = "black"), # Panel background and border
    panel.grid.major = element_blank(), # Remove major grid lines
    panel.grid.minor = element_blank(), # Remove minor grid lines
    axis.text = element_text(size = 9, color = "black"), # Axis text style
    axis.title = element_text(size = 9, color = "black"), # Axis title style
    plot.title = element_text(size = 10, hjust = 0.5, face = "bold"), # Plot title style (centered)
    legend.position = "none", # Remove legend (sample size alpha legend)
    panel.spacing = unit(1, "lines"), # Spacing between facets
    strip.background = element_rect(fill = "grey90", color = "black"), # Background for facet labels
    strip.text = element_text(size = 9, face = "bold", color = "black") # Style facet labels
  ) +
  # Create facets for Class and Sector
  facet_grid(Class ~ Sector, scales = "free_y", switch = "y") # Facet by Class (rows) and Sector (columns)

# Print the time-trend plot
print(timetrend_plot)

# Save the time-trend plot to a file (e.g., PDF or PNG)
cat("Saving time-trend plot...\n")
ggsave("./results/temporal_trend_plot.png", timetrend_plot, width = 8, height = 6, dpi = 300)
# Or save as PDF: ggsave("./results/temporal_trend_plot.pdf", timetrend_plot, width = 8, height = 6)
cat("Time-trend plot saved to ./results/temporal_trend_plot.png\n")


# --- Analysis Section 2: Overall Pooled Analysis by Bacteria, Sector, Antibiotic ---
## for all following analysis, "rawp"= "Proportion", "E"=	"Resistant isolates (n)" and "N"= "Isolates (N)"

cat("Performing overall pooled meta-analysis...\n")

# Load data from the specified sheet (Table S1)
df_overall <- read_excel(data_file_path, sheet = overall_country_sheet)

# No filtering based on WHO classification as requested
df_filtered_overall <- df_overall # Use the entire dataframe

# Split data by Bacteria, Sector, and Antibiotic for meta-analysis
subsets_overall <- split(df_filtered_overall, list(df_filtered_overall$Bacteria, df_filtered_overall$Sector, df_filtered_overall$Antibiotic))

# Initialize an empty data frame to store results
results_df_overall <- data.frame()

# Loop through each subset and perform meta-analysis if enough studies exist
for (subset_data in subsets_overall) {
  # Check if there are enough studies (more than 1) for meta-analysis
  if (nrow(subset_data) > 1) {
    # Perform meta-analysis using Freeman-Tukey transformation (sm="PAS") and Inverse Variance method
    meta_result <- metaprop(event = E, n = N, data = subset_data, sm = "PAS", method = "Inverse")
    
    # Get the summary of meta-analysis results
    summary_meta <- summary(meta_result)
    
    # Extract pooled estimate (on transformed scale) and its confidence interval
    pooled_proportion_trans <- summary_meta$TE.random
    lower_ci_trans <- summary_meta$lower.random
    upper_ci_trans <- summary_meta$upper.random
    
    # Back-transform pooled estimate and CIs to the original proportion scale (using sin^2)
    pooled_proportion <- sin(pooled_proportion_trans)^2
    lower_ci <- sin(lower_ci_trans)^2
    upper_ci <- sin(upper_ci_trans)^2
    
    # Calculate summary statistics for the group (total events, total N, number of studies, I2 heterogeneity)
    event_n <- sum(subset_data$E, na.rm = TRUE)
    total_n <- sum(subset_data$N, na.rm = TRUE)
    study_n <- summary_meta$k
    I2 <- summary_meta$I2
    
    # Format the pooled result and its confidence interval as a string
    formatted_result <- paste0(
      sprintf("%.1f%%", pooled_proportion * 100), " (",
      sprintf("%.1f", lower_ci * 100), "-",
      sprintf("%.1f", upper_ci * 100), "%)"
    )
    
    # Store the results for this subset in the results data frame
    # Removed WHO column
    results_df_overall <- rbind(results_df_overall, data.frame(
      BacteriaType = subset_data$Bacteria[1], # Take the first value as representative for the group
      Sector = subset_data$Sector[1],
      # Removed WHO column
      DrugClass = subset_data$Antibiotic[1],
      FormattedResult = formatted_result,
      raw_proportion = pooled_proportion, # Store the raw pooled proportion for plotting
      event_n = event_n,
      total_n = total_n,
      study_n = study_n,
      I2 = I2
    ))
  }
}

# Save the overall meta-analysis results to an Excel file in the 'results' folder
write.xlsx(results_df_overall, "./results/meta_analysis_overall.xlsx")

cat("Overall pooled meta-analysis complete. Results saved to ./results/meta_analysis_overall.xlsx\n")

# --- Figures - Heatmap (Overall) ---

cat("Generating overall heatmap figure...\n")

# Adjust the order of 'Sector' for consistent plotting
results_df_overall$Sector <- factor(
  results_df_overall$Sector,
  levels = c("Human", "Animal", "Environment", "Food") # Define desired order
)

# Reorder the DrugClass factor for plotting (ordering based on unique values)
results_df_overall <- results_df_overall %>%
  mutate(
    DrugClass = factor(DrugClass, levels = unique(DrugClass)) # Reorder DrugClass based on unique values
  )

# Define labels for bacteria types, ensuring scientific names are italicized in the plot
bacteria_labels <- c(
  "Escherichia coli" = expression(paste(italic("Escherichia coli"))),
  "Salmonella spp." = expression(paste(italic("Salmonella"), " spp.")),
  "Klebsiella pneumoniae" = expression(paste(italic("Klebsiella pneumoniae"))),
  "Klebsiella spp." = expression(paste(italic("Klebsiella"), " spp.")),
  "Enterobacter spp." = expression(paste(italic("Enterobacter"), " spp.")),
  "Nontyphoidal Salmonella" = expression(paste("Nontyphoidal ", italic("Salmonella"))) # Add space for clarity
)

# Create the overall heatmap plot using ggplot2
overall_heatmap_plot <- ggplot(results_df_overall, aes(x = DrugClass, y = BacteriaType)) +
  # Use geom_point to create the heatmap effect with size and color mapping
  geom_point(aes(size = total_n, color = raw_proportion), alpha = 0.8) +
  # Set color scale using a viridis palette (magma, reversed)
  scale_color_viridis_c(option = "magma", direction = -1, labels = percent_format(scale = 100)) +
  # Set size scale for points based on total sample size
  scale_size(
    range = c(1, 6), # Adjusted size range slightly for potentially better visual
    breaks = c(100, 1000, 5000, 10000, 20000), # Adjusted breaks, add more if needed based on data
    limits = c(1, NA) # Set lower limit for size, upper is determined by data
  ) +
  # Create separate panels (facets) for each Sector, with free y-axis scales
  facet_wrap(~ Sector, scales = "free_y", ncol = 1) +
  # Apply a minimal theme
  theme_minimal() +
  # Customize theme elements for appearance
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8), # Rotate x-axis labels
    axis.text.y = element_text(size = 8),
    strip.background = element_blank(), # Remove background from facet labels
    strip.text = element_text(size = 8, face = "bold"), # Style facet labels
    panel.background = element_rect(fill = "white", color = "grey", size = 0.1), # Panel background and border
    panel.grid.major = element_line(color = "grey90", size = 0.1), # Major grid lines
    panel.grid.minor = element_line(color = "grey90", size = 0.1), # Minor grid lines
    legend.position = "right", # Position the legend
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9),
    plot.title = element_text(face = "bold", size = 10), # Plot title style
    plot.subtitle = element_text(size = 8) # Plot subtitle style
    # Font families removed for portability
  ) +
  # Set plot labels and title
  labs(
    size = "Sample size",
    color = "Pooled Proportion", # Label for the color legend
    x = "Antibiotics",
    y = "Bacteria",
    title = "Overall Antibiotic Resistance Heatmap"
  ) +
  # Customize legend appearance (ensure full alpha for size legend points)
  guides(size = guide_legend(override.aes = list(alpha = 1))) +
  # Apply the defined bacteria labels to the y-axis
  scale_y_discrete(labels = bacteria_labels)

# Print the plot (useful when running interactively)
print(overall_heatmap_plot)

# Save the overall heatmap plot to a PowerPoint file
cat("Saving overall heatmap figure to PowerPoint...\n")
ppt_overall <- read_pptx()
# Add a slide, assuming a standard layout with a title and content placeholder
slide_overall <- add_slide(ppt_overall, layout = "Title and Content", master = "Office Theme")

# Convert the ggplot object to a drawable object (dml) for PowerPoint
# Adjust width and height as needed to fit the slide
plot_svg_overall <- dml(ggobj = overall_heatmap_plot, width = 15, height = 7)
slide_overall <- ph_with(slide_overall, plot_svg_overall, location = ph_location(left = 1, top = 1.5, width = 15, height = 7))

# Print (save) the PowerPoint object
print(ppt_overall, target = "./results/antibiotic_resistance_heatmap_overall.pptx")

cat("Overall heatmap figure saved to ./results/antibiotic_resistance_heatmap_overall.pptx\n")

# Clean up objects from this section if memory is a concern (optional)
# rm(ppt_overall, slide_overall, plot_svg_overall)


# --- Analysis Section 3: Subgroup Analysis by Country ---

cat("Performing subgroup analysis by country...\n")

# Load data for country subgroup analysis (assuming same data file, same sheet as overall analysis)
df_country <- read_excel(data_file_path, sheet = overall_country_sheet) # Or load from a different sheet/file if needed

# No filtering based on WHO classification as requested
df_filtered_country <- df_country # Use the entire dataframe

# Split data by Country, Bacteria, Sector, and Antibiotic for meta-analysis
subsets_country <- split(df_filtered_country, list(df_filtered_country$Country, df_filtered_country$Bacteria, df_filtered_country$Sector, df_filtered_country$Antibiotic))

# Initialize an empty data frame to store results
results_df_country <- data.frame()

# Loop through each subset and perform meta-analysis if enough studies exist
for (subset_data in subsets_country) {
  # Check if there are enough studies (more than 1) for meta-analysis
  if (nrow(subset_data) > 1) {
    # Perform meta-analysis (Freeman-Tukey, Inverse)
    meta_result <- metaprop(event = E, n = N, data = subset_data, sm = "PAS", method = "Inverse")
    
    # Get the summary
    summary_meta <- summary(meta_result)
    
    # Extract and back-transform proportions and CIs
    pooled_proportion_trans <- summary_meta$TE.random
    lower_ci_trans <- summary_meta$lower.random
    upper_ci_trans <- summary_meta$upper.random
    
    pooled_proportion <- sin(pooled_proportion_trans)^2
    lower_ci <- sin(lower_ci_trans)^2
    upper_ci <- sin(upper_ci_trans)^2
    
    # Calculate summary stats
    event_n <- sum(subset_data$E, na.rm = TRUE)
    total_n <- sum(subset_data$N, na.rm = TRUE)
    study_n <- summary_meta$k
    I2 <- summary_meta$I2
    
    # Format result string
    formatted_result <- paste0(
      sprintf("%.1f%%", pooled_proportion * 100), " (",
      sprintf("%.1f", lower_ci * 100), "-",
      sprintf("%.1f", upper_ci * 100), "%)"
    )
    
    # Store results
    results_df_country <- rbind(results_df_country, data.frame(
      Country = subset_data$Country[1],
      BacteriaType = subset_data$Bacteria[1],
      Sector = subset_data$Sector[1],
      DrugClass = subset_data$Antibiotic[1],
      FormattedResult = formatted_result,
      raw_proportion = pooled_proportion, # Store raw proportion for plotting
      event_n = event_n,
      total_n = total_n,
      study_n = study_n,
      I2 = I2
    ))
  }
}

# Save results to Excel
write.xlsx(results_df_country, "./results/meta_analysis_country.xlsx")

cat("Subgroup analysis by country complete. Results saved to ./results/meta_analysis_country.xlsx\n")

# --- Analysis Section 4: Gene Analysis ---

cat("Performing pooled meta-analysis for resistance genes...\n")

# Load data specifically for gene analysis from the specified file and sheet
# Data columns expected: E (Resistant isolates), N (Isolates), gc (Gene Class/Type), Setting (Study Setting)
dfg <- read_excel(data_file_path, sheet = gene_data_sheet)

# Split data by Gene Class (gc) and Setting for meta-analysis
# This creates subsets for each unique combination of gc and Setting
subsets_gene <- split(dfg, list(dfg$gc, dfg$Setting))

# Initialize an empty data frame to store meta-analysis results for genes
results_dfg <- data.frame()

# Loop through each subset (combination of Gene Class and Setting)
for (subset_data in subsets_gene) {
  # Check if there are enough studies (more than 1 row) in the subset for meta-analysis
  if (nrow(subset_data) > 1) {
    # Perform meta-analysis using Freeman-Tukey transformation (sm="PAS") and Inverse Variance method
    # event = E (Number of resistant isolates)
    # n = N (Total number of isolates)
    meta_result <- metaprop(event = E, n = N, data = subset_data, sm = "PAS", method = "Inverse")
    
    # Get the summary of meta-analysis results
    summary_meta <- summary(meta_result)
    
    # Extract the pooled estimate (on transformed scale) and its confidence interval
    pooled_proportion_trans <- summary_meta$TE.random
    lower_ci_trans <- summary_meta$lower.random
    upper_ci_trans <- summary_meta$upper.random
    
    # Back-transform the pooled estimate and CIs to the original proportion scale (using sin^2)
    pooled_proportion <- sin(pooled_proportion_trans)^2
    lower_ci <- sin(lower_ci_trans)^2
    upper_ci <- sin(upper_ci_trans)^2
    
    # Calculate summary statistics for the group (total events, total N, number of studies, I2 heterogeneity)
    # Note: Summing E and N across the subset_data (which is for a specific gc and Setting combination)
    event_n <- sum(subset_data$E, na.rm = TRUE)
    total_n <- sum(subset_data$N, na.rm = TRUE)
    study_n <- summary_meta$k # Number of studies in this subset
    I2 <- summary_meta$I2 # Heterogeneity I2 statistic
    
    # Format the pooled result and its confidence interval as a string
    formatted_result <- sprintf("%.1f%% (%.1f-%.1f%%)", pooled_proportion * 100, lower_ci * 100, upper_ci * 100)
    
    # Store the results for this subset in the results data frame
    results_dfg <- rbind(results_dfg, data.frame(
      Gene = subset_data$gc[1], # Take the first value as representative for the group
      Setting = subset_data$Setting[1], # Take the first value as representative for the group
      FormattedResult = formatted_result,
      raw_proportion = pooled_proportion, # Store the raw pooled proportion for potential further use/plotting
      event_n = event_n,
      total_n = total_n,
      study_n = study_n,
      I2 = I2
    ))
  }
}

# Save the gene meta-analysis results to an Excel file in the 'results' folder
write.xlsx(results_dfg, "./results/meta_analysis_gene.xlsx") # Changed filename slightly for consistency

cat("Pooled meta-analysis for resistance genes complete. Results saved to ./results/meta_analysis_gene.xlsx\n")


# --- Figures - Chord Plot Visualization of Gene-Setting Links ---

cat("Generating chord plot for gene-setting links...\n")

# Load the data again for the chord plot (assuming the same data source and sheet)
# Data columns expected: gc (Gene Class/Type), Setting (Study Setting), E (Number of resistant isolates)
dfg_chord <- read_excel(data_file_path, sheet = gene_data_sheet)

# Create a data frame for the chord diagram
# Using E (Resistant isolates) to represent the link strength
chord_data <- dfg_chord %>%
  select(gc, Setting, E) %>%
  # Filter out rows where E is missing or zero, as they won't contribute to links
  filter(!is.na(E) & E > 0)

# Aggregate data if multiple entries exist for the same gene-setting pair
# Summing E for each unique combination of gc and Setting
chord_data_aggregated <- chord_data %>%
  group_by(gc, Setting) %>%
  summarise(Total_E = sum(E, na.rm = TRUE)) %>%
  ungroup() # Ungroup after summarizing

# Sort the data by the aggregated total number of resistant isolates (Total_E) in descending order
chord_data_aggregated <- chord_data_aggregated[order(-chord_data_aggregated$Total_E), ]

# Define the sectors (all unique gene classes and settings)
sectors <- unique(c(chord_data_aggregated$gc, chord_data_aggregated$Setting))

# Define colors for the sectors
num_colors <- length(sectors)
# Use a color palette, ensure enough colors for all sectors
mycolor <- brewer.pal(min(num_colors, 11), "Spectral") # Use up to 11 distinct colors from Spectral palette

# If there are more sectors than colors in the chosen palette, repeat the colors
if (num_colors > length(mycolor)) {
  warning("Not enough distinct colors in the palette for all sectors. Colors will be repeated.")
  mycolor <- rep(mycolor, length.out = num_colors)
}

# Assign colors to sectors (genes and settings)
# Create a named vector where names are sector names and values are colors
sector_colors <- setNames(mycolor, sectors)

# Initialize the chord diagram layout
circos.clear() # Clear any previous circos plots
par(family = "Times") # Set font family for plot labels

# Create the chord diagram
# x = data frame with links (from, to, value)
# grid.col = colors for each sector
# transparency = transparency level for links
# annotationTrack = add a track for sector names/labels
# preAllocateTracks = list of track heights (adjust to control ring width)
# directional = 1 means links are directional (from gc to Setting)
# direction.type = type of directional indicator (arrows and height difference)
# link.arr.type = style of the arrow
# link.sort = sort links within sectors
# link.largest.ontop = draw larger links on top
# diffHeight = difference in height between 'from' and 'to' sides of links
# link.lwd = line width of links, mapped to Total_E (number of resistant isolates)
chordDiagram(
  x = chord_data_aggregated,
  grid.col = sector_colors,
  transparency = 0.5,
  annotationTrack = "grid",
  preAllocateTracks = list(track.height = 0.05), # Adjust track height to control ring width
  directional = 1, # Links go from gene (gc) to setting
  direction.type = c("arrows", "diffHeight"),
  link.arr.type = "big.arrow",
  link.sort = TRUE,
  link.largest.ontop = TRUE,
  diffHeight = 0.04,
  link.lwd = chord_data_aggregated$Total_E / max(chord_data_aggregated$Total_E, na.rm = TRUE) * 5 # Scale line width based on Total_E (adjust multiplier 5 as needed)
)

# Add labels to the sectors (gene classes and settings)
circos.trackPlotRegion(
  track.index = 1, # Target the first track (where sectors are)
  bg.border = NA, # No border for the background
  # Define a function to draw labels within each sector's panel
  panel.fun = function(x, y) {
    sector.name = get.cell.meta.data("sector.index") # Get the name of the current sector
    circos.text(
      x = mean(get.cell.meta.data("xlim")), # X position at the center of the sector
      y = get.cell.meta.data("ylim")[1] + 0.1, # Y position slightly outside the inner edge
      labels = sector.name, # The label text
      facing = "clockwise", # Text orientation
      niceFacing = TRUE, # Adjust orientation for readability
      adj = c(0, 0.5), # Alignment of the text relative to the position
      cex = 0.6 # Font size (adjust as needed)
    )
  },
  track.height = 0.05 # Match the track height used in chordDiagram
)

# Add a color legend for the settings
# Extract colors and names specifically for settings from the sector_colors vector
setting_colors_legend <- sector_colors[unique(chord_data_aggregated$Setting)]

# Add the legend to the plot area
# Position: "topright"
# legend: names of the settings
# fill: corresponding colors for settings
# title: Legend title
# cex: Font size for legend text
# bty: "n" means no box around the legend
# text.font: 1 corresponds to plain text, often used with par(family)
legend(
  "topright",
  legend = names(setting_colors_legend),
  fill = setting_colors_legend,
  title = "Settings",
  cex = 0.6, # Font size for legend text (adjusted slightly)
  text.font = 1, # Use the font specified in par(family)
  bty = "n" # No box around the legend
)

# Note: Saving chord plots directly to file can be tricky with base R graphics and circos.
# You might need to use a graphics device like png(), pdf(), or svg() before calling circos.clear()
# and chordDiagram(), and then dev.off() after the plot is complete.

# Example of saving to PNG:
# png("./results/gene_setting_chord_plot.png", width = 800, height = 800, res = 150)
# circos.clear()
# par(family = "Times")
# # ... (your chordDiagram and circos.trackPlotRegion code here) ...
# legend(...)
# dev.off()


cat("Chord plot generation complete. Consider saving using a graphics device like png() or pdf().\n")


# --- End of Script ---

