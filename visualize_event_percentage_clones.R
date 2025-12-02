


library(dplyr)
library(purrr)

# Define sample names
sample_names <- c("C9", "B4", "D5", "D7", "D11", "G3")

# Load RDS files dynamically into a named list
df_list <- setNames(
  lapply(sample_names, function(sample) {
    readRDS(paste0("../../Volumes/Fl/scUnique_results/SLX-25394_", sample, "-15-qc_100/SLX-25394_", sample, "-15-qc_100.df_pass_post.RDS"))
  }),
  sample_names
)

# Load the parent dataset
df_pass_post <- readRDS("../../Volumes/Fl/scUnique_results/SLX-25393-15-qc_100/SLX-25393-15-qc_100.df_pass_post.RDS")

# Perform joins using loops
parent_joins <- map(df_list, ~ inner_join(df_pass_post, .x, by = c("start", "end", "chroms_l")))
exact_parent_joins <- map(df_list, ~ inner_join(df_pass_post, .x, by = c("start", "end", "chroms_l", "h0_l", "h1_l")))

# Pairwise joins between all samples
sample_combinations <- combn(sample_names, 2, simplify = FALSE)

exact_pairwise_joins <- map(sample_combinations, function(pair) {
  inner_join(df_list[[pair[1]]], df_list[[pair[2]]], by = c("start", "end", "chroms_l", "h0_l", "h1_l"))
})

pairwise_joins <- map(sample_combinations, function(pair) {
  inner_join(df_list[[pair[1]]], df_list[[pair[2]]], by = c("start", "end", "chroms_l"))
})

# Convert lists to named elements
names(parent_joins) <- paste0("parent_", sample_names)
names(exact_parent_joins) <- paste0("exact_parent_", sample_names)
names(exact_pairwise_joins) <- paste0("exact_", sapply(sample_combinations, paste, collapse = "_"))
names(pairwise_joins) <- paste0(sapply(sample_combinations, paste, collapse = "_"))

# Access results like:
# parent_joins$parent_C9
# exact_pairwise_joins$exact_C9_D5

# Summarize group counts for each join using pairwise_joins
summary_table <- map_dfr(parent_joins, ~ 
                           .x %>%
                           group_by(start, end, chroms_l) %>%
                           summarise(count = n(), .groups = "drop"), 
                         .id = "Sample"
)

# Print or save the summarized table
print(summary_table)





# Count the number of unique groups for each pair in pairwise_joins
parent_group_counts <- map_df(parent_joins, ~ 
                         tibble(num_groups = n_distinct(.x$start, .x$end, .x$chroms_l)), 
                       .id = "Pair"
)

num_groups_df_pass_post <- df_pass_post %>%
  distinct(start, end, chroms_l) %>%
  nrow()
parent_group_counts <- parent_group_counts %>%
  mutate(percentage = (num_groups / num_groups_df_pass_post) * 100)
# Print the result
print(parent_group_counts)




# Count the number of unique groups for each pair in pairwise_joins
pairwise_group_counts <- map_df(pairwise_joins, ~ 
                                tibble(num_groups = n_distinct(.x$start, .x$end, .x$chroms_l)), 
                              .id = "Pair"
)

# List of the datasets (C9, D7, ...)
datasets <- names(df_list)

# Function to calculate the number of unique groups for each dataset
get_num_groups <- function(dataset) {
  df_list[[dataset]] %>%
    distinct(start, end, chroms_l) %>%
    nrow()
}

# Calculate the number of unique groups for each dataset
num_groups_list <- sapply(datasets, get_num_groups, simplify = FALSE)

# Iterate over the datasets and add the percentage columns to pairwise_group_counts
for (dataset in datasets) {
  pairwise_group_counts <- pairwise_group_counts %>%
    mutate(!!paste0("percentage_", dataset) := (num_groups / num_groups_list[[dataset]]) * 100)
}

# Print the updated pairwise_group_counts to check the result
print(pairwise_group_counts)





library(ggplot2)
library(reshape2)

# Example: Similarity values between samples (use your actual similarity values)

#  p.   C9.     B4.      D5.    D7.    D11.  G3
similarity_matrix <- matrix(c( # NA,1.06,    0.579,  1.16,  3.47,  3.67,  3.76,
                               1.06, NA,   1.2,    2.46,  4.52 , 0.730 , 1.39,
                               0.579, 0.855, NA,  1.64,  3.23,  0.547 , 0.695,
                                1.16, 2.56,   2.41, NA,   3.23 ,  1.09 , 1.53,
                                3.47, 5.98,   6.02, 4.10,  NA,    2.19 ,2.64,
                                3.67, 3.42,   3.61,4.92,  7.74,   NA   ,5.70,
                                3.76, 8.55,   6.02, 9.02 , 12.3 , 7.48 ,NA), 
                            nrow = 6, byrow = TRUE)

rownames(similarity_matrix) <- c( "C9", "B4", "D5", "D7", "D11","G3")
colnames(similarity_matrix) <- c("Parent", "C9", "B4", "D5", "D7", "D11","G3")

# Convert the matrix to a data frame for ggplot
similarity_df <- melt(similarity_matrix)

# Plot heatmap
ggplot(similarity_df, aes(Var2, Var1, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 3) +
  theme_minimal() +
  labs(title = "Event Similarity Between Samples", x = "Sample", y = "Sample") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

