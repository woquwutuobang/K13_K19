library(bio3d)
library(ggplot2)
library(cowplot)

# ========= 1. Compute Distance Matrix from PDB =========

#' Compute Distance Matrix from PDB Structure
#'
#' This function reads a PDB file and computes a pairwise distance matrix
#' between specified atoms (typically C-alpha atoms) for a given chain.
#' The distance matrix is used for spatial clustering analysis of protein residues.
#'
#' @param pdb_file Character string. Path to the PDB structure file.
#' @param chain Character string. Chain identifier to analyze (default: "A").
#' @param atom_type Character string. Type of atoms to use for distance calculation 
#'   (default: "CA" for C-alpha atoms).
#'
#' @return A symmetric distance matrix with residue numbers as row and column names.
#'   Distances are in Angstroms (Å).
#'
#' @examples
#' # Compute distance matrix for chain A using C-alpha atoms
#' dist_matrix <- compute_dist_matrix("protein.pdb", chain = "A", atom_type = "CA")
#' 
#' # Compute distance matrix for chain B using C-beta atoms
#' dist_matrix_b <- compute_dist_matrix("protein.pdb", chain = "B", atom_type = "CB")
#' 
#' # View matrix dimensions and first few residues
#' dim(dist_matrix)
#' head(rownames(dist_matrix))

compute_dist_matrix <- function(pdb_file, chain = "A", atom_type = "CA") {
  pdb <- read.pdb(pdb_file)
  inds <- atom.select(pdb, chain = chain, elety = atom_type)
  if (length(inds$xyz) == 0) {
    stop(paste("No", atom_type, "atoms found in chain", chain))
  }
  xyz_mat <- matrix(pdb$xyz[inds$xyz], ncol = 3, byrow = TRUE)
  dist_matrix <- as.matrix(dist(xyz_mat))
  residues <- pdb$atom[inds$atom, "resno"]
  if (length(residues) != nrow(dist_matrix)) {
    warning("Residue number mismatch, adjusting automatically.")
    residues <- residues[1:nrow(dist_matrix)]
  }
  rownames(dist_matrix) <- colnames(dist_matrix) <- as.character(residues)
  return(dist_matrix)
}

# ========= 2. Calculate Four Spatial Clustering Metrics =========

#' Calculate Spatial Clustering Metrics for Hotspot Residues
#'
#' This function computes four spatial clustering metrics by comparing actual
#' hotspot residue distances with random distributions. It calculates pairwise
#' distances, median distances, per-site minimum distances, and median per-site
#' minimum distances for hotspot residues versus random residue sets.
#'
#' @param hotspots Numeric vector. Residue numbers of hotspot residues to analyze.
#' @param dist_matrix Matrix. Distance matrix computed from PDB structure with
#'   residue numbers as row/column names.
#' @param n_sim Integer. Number of random simulations to perform for comparison
#'   (default: 1000).
#'
#' @return A list containing:
#'   \item{actual_distances}{Numeric vector of pairwise distances between hotspot residues}
#'   \item{actual_median}{Median of pairwise distances between hotspot residues}
#'   \item{per_site_min}{Numeric vector of minimum distances for each hotspot residue}
#'   \item{median_min}{Median of per-site minimum distances}
#'   \item{sim_dists}{Numeric vector of all pairwise distances from random simulations}
#'   \item{sim_medians}{Numeric vector of median distances from each random simulation}
#'   \item{sim_permin}{Numeric vector of all per-site minimum distances from random simulations}
#'   \item{sim_median_mins}{Numeric vector of median per-site minimum distances from each random simulation}
#'
#' @examples
#' # Load distance matrix
#' dist_matrix <- compute_dist_matrix("protein.pdb")
#' 
#' # Define hotspot residues
#' hotspots <- c(15, 16, 145, 10, 53, 77, 89, 151)
#' 
#' # Calculate clustering metrics
#' clustering_data <- compute_distributions(hotspots, dist_matrix, n_sim = 1000)
#' 
#' # Access results
#' actual_median <- clustering_data$actual_median
#' random_medians <- clustering_data$sim_medians
#' 
#' # Compare actual vs random
#' p_value <- sum(random_medians <= actual_median) / length(random_medians)

compute_distributions <- function(hotspots, dist_matrix, n_sim = 1000) {
  hotspot_ids <- as.character(intersect(hotspots, as.numeric(rownames(dist_matrix))))
  if (length(hotspot_ids) < 2) {
    stop(paste("Insufficient hotspot residues:", length(hotspot_ids)))
  }
  sub_dist <- dist_matrix[hotspot_ids, hotspot_ids]
  sub_dist[sub_dist == 0] <- NA
  actual_distances <- na.omit(as.vector(sub_dist))
  actual_median <- median(actual_distances)
  per_site_min <- apply(sub_dist, 1, function(x) min(x, na.rm = TRUE))
  median_min <- median(per_site_min)
  
  sim_medians <- numeric(n_sim)
  sim_median_mins <- numeric(n_sim)
  sim_dists <- list()
  sim_permin <- list()
  all_residues <- rownames(dist_matrix)
  
  for (i in 1:n_sim) {
    rand_ids <- sample(all_residues, length(hotspot_ids))
    rand_sub <- dist_matrix[rand_ids, rand_ids]
    rand_sub[rand_sub == 0] <- NA
    rand_dists <- na.omit(as.vector(rand_sub))
    sim_dists[[i]] <- rand_dists
    sim_medians[i] <- median(rand_dists)
    rand_permin <- apply(rand_sub, 1, function(x) min(x, na.rm = TRUE))
    sim_permin[[i]] <- rand_permin
    sim_median_mins[i] <- median(rand_permin)
  }
  
  return(list(
    actual_distances = actual_distances,
    actual_median = actual_median,
    per_site_min = per_site_min,
    median_min = median_min,
    sim_dists = unlist(sim_dists),
    sim_medians = sim_medians,
    sim_permin = unlist(sim_permin),
    sim_median_mins = sim_median_mins
  ))
}


# ========= 3. Plot Spatial Clustering Analysis Results =========


#' Plot Spatial Clustering Analysis Results
#'
#' This function creates a comprehensive visualization of spatial clustering analysis
#' results, generating four plots: pairwise distance distributions, median pairwise
#' distances, per-site minimum distances, and median per-site minimum distances.
#' Each plot compares actual hotspot distributions with random distributions.
#'
#' @param assay Character string. Name of the binding partner/assay being analyzed.
#' @param dist_data List. Output from \code{compute_distributions()} containing
#'   actual and simulated distance data.
#' @param outdir Character string. Output directory for saving plots (default:
#'   "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/figure_s6/20251011/allosteric_clustering_results").
#'
#' @return Invisibly returns the combined plot object. Saves a PDF file with
#'   four subplots showing spatial clustering analysis results.
#'
#' @examples
#' # Calculate clustering data first
#' dist_matrix <- compute_dist_matrix("protein.pdb")
#' hotspots <- c(15, 16, 145, 10, 53, 77, 89, 151)
#' clustering_data <- compute_distributions(hotspots, dist_matrix)
#' 
#' # Create plots
#' plot_connectivity("K13", clustering_data)
#' 
#' # Custom output directory
#' plot_connectivity("K19", clustering_data, 
#'                  outdir = "C:/path/to/custom/output")
#' 
#' # The function saves a PDF file named "spatial_clustering_[assay].pdf"

plot_spatial_clustering_analysis_results <- function(assay, dist_data, outdir = "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/figure_s6/20251022/allosteric_clustering_results") {
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  
  df_f <- data.frame(
    distance = c(dist_data$sim_dists, dist_data$actual_distances),
    group = rep(c("Random", "Hotspot"), 
                c(length(dist_data$sim_dists), length(dist_data$actual_distances)))
  )
  
  df_h <- data.frame(
    distance = c(dist_data$sim_permin, dist_data$per_site_min),
    group = rep(c("Random", "Hotspot"), 
                c(length(dist_data$sim_permin), length(dist_data$per_site_min)))
  )
  
  base_theme <- theme_classic(base_size = 8) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
      axis.text.y = element_text(size = 8),
      axis.title = element_text(size = 8),
      plot.title = element_text(size = 8),
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 8),
      legend.position = c(0.8, 0.8)
    )
  
  # a. Pairwise distances
  p_f <- ggplot(df_f, aes(x = distance, color = group)) +
    geom_density() +
    scale_color_manual(name = "Distribution", values = c("Random" = "black", "Hotspot" = "#FF6A56")) +
    labs(title = "a. Pairwise distances", x = "Distance (Å)", y = "Density") +
    base_theme
  
  # b. Median pairwise distance
  p_g <- ggplot() +
    geom_density(aes(x = dist_data$sim_medians, color = "Random")) +
    geom_vline(aes(xintercept = dist_data$actual_median, color = "Hotspot"),
               size = 0.8, linetype = "solid") +
    scale_color_manual(name = "Distribution",
                       values = c("Random" = "black", "Hotspot" = "#FF6A56"),
                       labels = c("Random", "Hotspot median")) +
    labs(title = "b. Median pairwise distance", x = "Median distance (Å)", y = "Density") +
    base_theme
  
  # c. Per-site minimum distance
  p_h <- ggplot(df_h, aes(x = distance, color = group)) +
    geom_density() +
    scale_color_manual(name = "Distribution", values = c("Random" = "black", "Hotspot" = "#FF6A56")) +
    labs(title = "c. Minimum distance per site", x = "Minimum distance (Å)", y = "Density") +
    base_theme
  
  # d. Median per-site minimum distance
  p_i <- ggplot() +
    geom_density(aes(x = dist_data$sim_median_mins, color = "Random")) +
    geom_vline(aes(xintercept = dist_data$median_min, color = "Hotspot"),
               size = 0.8, linetype = "solid") +
    scale_color_manual(name = "Distribution",
                       values = c("Random" = "black", "Hotspot" = "#FF6A56"),
                       labels = c("Random", "Hotspot median")) +
    labs(title = "d. Median minimum distance per site", x = "Median minimum distance (Å)", y = "Density") +
    base_theme
  
  p_all <- plot_grid(p_f, p_g, p_h, p_i, nrow = 1)
  outfile <- file.path(outdir, paste0("spatial_clustering_", assay, ".pdf"))
  ggsave(outfile, p_all, width = 16, height = 4)
  message("Saved:", outfile)
}

# ========= 4. Main Analysis Pipeline =========

#' Run Complete Allosteric Clustering Analysis
#'
#' This is the main function that orchestrates the complete allosteric clustering
#' analysis pipeline. It computes distance matrices from PDB structures and
#' analyzes spatial clustering patterns for multiple binding partners/assays.
#' For each assay, it calculates clustering metrics and generates visualization plots.
#'
#' @param hotspot_list Named list. Each element contains a numeric vector of
#'   hotspot residue numbers for a specific binding partner. Names should be
#'   assay/binding partner identifiers.
#' @param pdb_file Character string. Path to the PDB structure file.
#' @param chain Character string. Chain identifier to analyze (default: "A").
#' @param atom_type Character string. Type of atoms to use for distance calculation
#'   (default: "CA" for C-alpha atoms).
#' @param n_sim Integer. Number of random simulations to perform for comparison
#'   (default: 1000).
#'
#' @return Invisibly returns NULL. The function processes each assay in the
#'   hotspot_list, computes clustering metrics, and saves visualization plots
#'   to the output directory.
#'
#' @examples
#' # Define hotspot residues for different binding partners
#' hotspot_list <- list(
#'   K13 = c(15, 16, 145, 10, 53, 77, 89, 151),
#'   K19 = c(15, 16, 145, 10, 77, 89, 134, 151),
#'   RAF1 = c(15, 16, 17, 28, 32, 35, 57, 60, 145, 146, 10, 54, 58, 59, 144, 163)
#' )
#' 
#' # Run analysis
#' run_allosteric_clustering_analysis(
#'   hotspot_list = hotspot_list,
#'   pdb_file = "protein.pdb",
#'   chain = "A",
#'   atom_type = "CA",
#'   n_sim = 1000
#' )
#' 
#' # Analyze specific chain with different atom types
#' run_allosteric_clustering_analysis(
#'   hotspot_list = hotspot_list,
#'   pdb_file = "protein.pdb",
#'   chain = "B",
#'   atom_type = "CB",
#'   n_sim = 500
#' )

run_allosteric_clustering_analysis <- function(hotspot_list, pdb_file, chain = "A", atom_type = "CA", n_sim = 1000) {
  dist_matrix <- compute_dist_matrix(pdb_file, chain = chain, atom_type = atom_type)
  cat("Distance matrix dimensions:", dim(dist_matrix), "\n")
  cat("First 10 residues:", head(rownames(dist_matrix), 10), "\n")
  
  for (assay in names(hotspot_list)) {
    cat("Processing assay:", assay, "\n")
    hotspots <- hotspot_list[[assay]]
    cat("Hotspot residues:", hotspots, "\n")
    
    tryCatch({
      dist_data <- compute_distributions(hotspots, dist_matrix, n_sim = n_sim)
      plot_spatial_clustering_analysis_results(assay, dist_data)
    }, error = function(e) {
      message("ERROR:Error processing ", assay, ": ", e$message)
    })
  }
}



# ========= 5. Usage Example =========
# Define allosteric hotspot residues for different binding partners
hotspot_list <- list(
  K13    = c(15, 16, 145, 10, 53, 77, 89, 151),
  K19    = c(15, 16, 145, 10, 77, 89, 134, 151),
  RAF1   = c(15, 16, 17, 28, 32, 35, 57, 60, 145, 146, 10, 54, 58, 59, 144, 163),
  K55    = c(15, 16, 17, 28, 32, 35, 57, 60, 117, 146, 58, 59, 63, 68, 69, 71, 72, 77, 78),
  K27    = c(15, 16, 17, 57, 60, 117, 119, 145, 10, 37, 63, 76, 77, 83, 85, 112, 156),
  RALGDS = c(15, 16, 17, 28, 32, 34, 35, 57, 60, 146, 10, 20, 58, 59, 85),
  PIK3CG = c(16, 17, 28, 32, 34, 35, 57, 60, 146, 20, 58, 59, 68, 77, 85),
  SOS1   = c(16, 17, 28, 32, 35, 57, 146, 40, 54, 55, 69, 85, 144, 148)
)

# PDB structure file
pdb_file <- "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/6vjj.pdb"

# Run analysis for all binding partners
run_allosteric_clustering_analysis(
  hotspot_list = hotspot_list, 
  pdb_file = pdb_file, 
  chain = "A", 
  atom_type = "CA", 
  n_sim = 1000
)




# Analysis Summary
cat("\nAllosteric Clustering Analysis Summary:\n")
cat("========================================\n")
cat("Analyzed", length(hotspot_list), "binding partners\n")
cat("Structure file:", pdb_file, "\n")
cat("Spatial metrics calculated:\n")
cat("  • Pairwise distance distributions\n")
cat("  • Median pairwise distances\n")
cat("  • Per-site minimum distances\n")
cat("  • Median per-site minimum distances\n")
cat("Random simulations:", 1000, "iterations\n")
cat("Output directory: path/to/allosteric_clustering_results\n")
