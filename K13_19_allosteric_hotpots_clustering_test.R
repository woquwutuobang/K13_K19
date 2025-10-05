library(bio3d)
library(ggplot2)
library(cowplot)

# ========= 1. Compute Distance Matrix from PDB =========
compute_dist_matrix <- function(pdb_file, chain = "A", atom_type = "CA") {
  # Read PDB file
  pdb <- read.pdb(pdb_file)
  
  # Select specified atoms
  inds <- atom.select(pdb, chain = chain, elety = atom_type)
  if (length(inds$xyz) == 0) {
    stop(paste("No", atom_type, "atoms found in chain", chain))
  }
  
  # Extract coordinate matrix
  xyz_mat <- matrix(pdb$xyz[inds$xyz], ncol = 3, byrow = TRUE)
  
  # Compute distance matrix
  dist_matrix <- as.matrix(dist(xyz_mat))
  
  # Get residue numbers
  residues <- pdb$atom[inds$atom, "resno"]
  
  # Ensure residue numbers match matrix dimensions
  if (length(residues) != nrow(dist_matrix)) {
    warning("Number of residues (", length(residues), 
            ") doesn't match distance matrix size (", nrow(dist_matrix), 
            "), automatically truncating/adjusting")
    residues <- residues[1:nrow(dist_matrix)]
  }
  
  # Set row and column names
  rownames(dist_matrix) <- colnames(dist_matrix) <- as.character(residues)
  return(dist_matrix)
}

# ========= 2. Calculate Four Spatial Clustering Metrics =========
compute_distributions <- function(hotspots, dist_matrix, n_sim = 1000) {
  # Ensure hotspot IDs are character type to match matrix row names
  hotspot_ids <- as.character(intersect(hotspots, as.numeric(rownames(dist_matrix))))
  
  if (length(hotspot_ids) < 2) {
    stop(paste("Insufficient hotspot residues:", length(hotspot_ids)))
  }
  
  # Extract sub-matrix for hotspots
  sub_dist <- dist_matrix[hotspot_ids, hotspot_ids]
  sub_dist[sub_dist == 0] <- NA  # Remove self-distances
  
  # f: All pairwise distances between hotspots
  actual_distances <- na.omit(as.vector(sub_dist))
  
  # g: Median pairwise distance
  actual_median <- median(actual_distances)
  
  # h: Per-site minimum distance to other hotspots
  per_site_min <- apply(sub_dist, 1, function(x) min(x, na.rm = TRUE))
  
  # i: Median per-site minimum distance
  median_min <- median(per_site_min)
  
  # ========== Random Sampling Simulation ==========
  sim_medians <- numeric(n_sim)
  sim_median_mins <- numeric(n_sim)
  sim_dists <- list()
  sim_permin <- list()
  
  all_residues <- rownames(dist_matrix)
  
  for (i in 1:n_sim) {
    # Randomly sample residues matching hotspot count
    rand_ids <- sample(all_residues, length(hotspot_ids))
    rand_sub <- dist_matrix[rand_ids, rand_ids]
    rand_sub[rand_sub == 0] <- NA
    
    # Calculate metrics for random sample
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

# ========= 3. Plot Connectivity Analysis Results =========
plot_connectivity <- function(assay, dist_data, outdir = "path/to/allosteric_clustering_results") {
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  
  # Prepare data frames for plotting
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
  
  # a: Pairwise distances distribution
  p_f <- ggplot(df_f, aes(x = distance, color = group)) +
    geom_density() +
    scale_color_manual(
      name = "Distribution", 
      values = c("Random" = "black", "Hotspot" = "red")
    ) +
    labs(
      title = "a. Pairwise distances between sites",
      x = "Distance (Å)", 
      y = "Density"
    ) +
    theme_classic() +
    theme(legend.position = c(0.8, 0.8))
  
  # b: Median pairwise distance comparison
  p_g <- ggplot() +
    geom_density(aes(x = dist_data$sim_medians, color = "Random")) +
    geom_vline(
      aes(xintercept = dist_data$actual_median, color = "Hotspot"), 
      size = 1, 
      linetype = "solid"
    ) +
    scale_color_manual(
      name = "Distribution", 
      values = c("Random" = "black", "Hotspot" = "red"),
      labels = c("Random", "Hotspot median")
    ) +
    labs(
      title = "b. Median pairwise distance",
      x = "Median distance (Å)", 
      y = "Density"
    ) +
    theme_classic() +
    theme(legend.position = c(0.8, 0.8))
  
  # c: Per-site minimum distance distribution
  p_h <- ggplot(df_h, aes(x = distance, color = group)) +
    geom_density() +
    scale_color_manual(
      name = "Distribution", 
      values = c("Random" = "black", "Hotspot" = "red")
    ) +
    labs(
      title = "c. Minimum distance per site",
      x = "Minimum distance (Å)", 
      y = "Density"
    ) +
    theme_classic() +
    theme(legend.position = c(0.8, 0.8))
  
  # d: Median per-site minimum distance comparison
  p_i <- ggplot() +
    geom_density(aes(x = dist_data$sim_median_mins, color = "Random")) +
    geom_vline(
      aes(xintercept = dist_data$median_min, color = "Hotspot"), 
      size = 1, 
      linetype = "solid"
    ) +
    scale_color_manual(
      name = "Distribution", 
      values = c("Random" = "black", "Hotspot" = "red"),
      labels = c("Random", "Hotspot median")
    ) +
    labs(
      title = "d. Median minimum distance per site",
      x = "Median minimum distance (Å)", 
      y = "Density"
    ) +
    theme_classic() +
    theme(legend.position = c(0.8, 0.8))
  
  # Combine all plots
  p_all <- plot_grid(p_f, p_g, p_h, p_i, nrow = 1)
  
  # Save combined plot
  outfile <- file.path(outdir, paste0("spatial_clustering_", assay, ".pdf"))
  ggsave(outfile, p_all, width = 16, height = 4)
  message("✅ Saved: ", outfile)
}

# ========= 4. Main Analysis Pipeline =========
run_allosteric_clustering_analysis <- function(hotspot_list, pdb_file, chain = "A", atom_type = "CA", n_sim = 1000) {
  # Compute distance matrix from PDB structure
  dist_matrix <- compute_dist_matrix(pdb_file, chain = chain, atom_type = atom_type)
  
  # Print matrix information
  cat("Distance matrix dimensions:", dim(dist_matrix), "\n")
  cat("First 10 residues:", head(rownames(dist_matrix), 10), "\n")
  
  # Analyze each binding partner
  for (assay in names(hotspot_list)) {
    cat("Processing assay:", assay, "\n")
    hotspots <- hotspot_list[[assay]]
    cat("Hotspot residues:", hotspots, "\n")
    
    tryCatch({
      # Compute spatial distributions
      dist_data <- compute_distributions(hotspots, dist_matrix, n_sim = n_sim)
      
      # Generate and save plots
      plot_connectivity(assay, dist_data)
      
    }, error = function(e) {
      message("❌ Error processing ", assay, ": ", e$message)
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
pdb_file <- "path/to/kras_wt_structure.pdb"

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