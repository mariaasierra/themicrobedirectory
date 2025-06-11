#Revisions for DATABASE 
#May 6th 2025

library(dplyr)
library(picante)
library(ape)
library(phytools)
library(stringr)
#TaxID from TMD
all_tmd_6.1=read.csv("TMD_V09.1.csv", stringsAsFactors = F)
taxid_tmd=select(all_tmd_6.1,scientific_name,taxonomy_id )
#writeLines(taxid_tmd$scientific_name, "tmd_taxa.txt")

#import tree
tree = read.newick("tree_tmd_taxa.txt")

# Create community matrix: all species present in one sample
comm_matrix <- matrix(1, nrow = 1, ncol = length(tree$tip.label))
colnames(comm_matrix) <- tree$tip.label
rownames(comm_matrix) <- "tmd"

# Compute Faith’s PD
pd_result <- pd(comm_matrix, tree, include.root = F)
print(pd_result)
# Add a second row to your community matrix
pd_result <- rbind(comm_matrix, na = rep(0, ncol(comm_matrix)))


ses.pd(
  samp = pd_result,
  tree = tree,
  null.model = "richness",
  runs = 999,
  include.root = F
)

run_pd_ses_by_sample_size <- function(tree, sizes = c(100, 500, 1000), runs = 999, null.model = "richness", seed = 42) {
  set.seed(seed)
  results <- list()
  
  for (size in sizes) {
    cat("Running for sample size:", size, "\n")
    
    # Randomly select species
    species_subset <- sample(tree$tip.label, size)
    
    # Build 1-row community matrix
    comm_matrix <- matrix(0, nrow = 1, ncol = length(tree$tip.label))
    colnames(comm_matrix) <- tree$tip.label
    rownames(comm_matrix) <- paste0("sample_", size)
    comm_matrix[1, species_subset] <- 1
    
    # Add dummy row to avoid ses.pd() bug
    dummy_matrix <- rbind(comm_matrix, dummy = rep(0, ncol(comm_matrix)))
    
    # Run ses.pd
    pd_res <- ses.pd(
      samp = dummy_matrix,
      tree = tree,
      null.model = null.model,
      runs = runs,
      include.root = FALSE
    )
    
    # Store just the real sample result
    results[[as.character(size)]] <- pd_res[1, ]
  }
  
  # Combine results into one data frame
  pd_summary <- do.call(rbind, results)
  pd_summary$sample_size <- as.integer(rownames(pd_summary))
  rownames(pd_summary) <- NULL
  return(pd_summary)
}

# Run for sizes: 100, 500, 1000, 5000
pd_summary <- run_pd_ses_by_sample_size(tree, sizes = c(100, 500, 1000, 5000))

# View results
print(pd_summary)


plot(pd_summary$sample_size, pd_summary$pd.obs.z,
     type = "b", pch = 16, col = "blue",
     xlab = "Sample Size", ylab = "PD Z-score",
     main = "Phylogenetic Diversity Z-score vs. Sample Size")
abline(h = c(-1.96, 1.96), lty = 2, col = "red")  # significance thresholds
#The consistent Z-scores suggest that this overdispersion is not an artifact of sample size it's a real pattern in our data

