#!/usr/bin/env Rscript
# =============================================================================
# Module 6 — Clustering comparison: function-based vs identity-based
# Produces functional clusters (GFC) and identity clusters, then quantifies
# divergence between the two topologies.
# =============================================================================
suppressPackageStartupMessages({
  library(argparse)
  library(cluster)
  library(stats)
})

parser <- ArgumentParser(description = "Compare functional vs identity clustering")
parser$add_argument("--trait-matrix", required = TRUE)
parser$add_argument("--orthogroups", required = TRUE)
parser$add_argument("--method", default = "apcluster")
parser$add_argument("--mantel-perms", type = "integer", default = 999)
parser$add_argument("--out-functional", required = TRUE)
parser$add_argument("--out-identity", required = TRUE)
parser$add_argument("--out-report", required = TRUE)
parser$add_argument("--out-metrics", required = TRUE)
args <- parser$parse_args()

# ---- Load trait matrix (genome x trait) ----
trait_mat <- read.delim(args$trait_matrix, row.names = 1, check.names = FALSE)
cat(sprintf("Trait matrix: %d genomes x %d traits\n", nrow(trait_mat), ncol(trait_mat)))

# ---- Load orthogroups presence/absence ----
og_mat <- read.delim(args$orthogroups, row.names = 1, check.names = FALSE)
cat(sprintf("Orthogroups: %d genomes x %d OGs\n", nrow(og_mat), ncol(og_mat)))

# ---- Distance matrices ----
# Use correlation-based distance for functional (GFC-style)
trait_dist <- as.dist(1 - cor(trait_mat, use = "pairwise.complete.obs"))
# Jaccard/binary distance for identity (OG presence/absence)
og_dist <- dist(og_mat, method = "binary")

# Replace NA/NaN/Inf
trait_dist[is.na(trait_dist) | is.infinite(trait_dist)] <- 1
og_dist[is.na(og_dist) | is.infinite(og_dist)] <- 1

# ---- Functional clustering (Affinity Propagation on correlation) ----
# Use negative squared Euclidean distance as similarity for apcluster
library(apcluster)

sim_trait <- negDistMat(as.matrix(trait_dist), r = 2)
ap_func <- apcluster(sim_trait, q = 0.5)
func_clusters <- data.frame(
  genome = rownames(trait_mat),
  functional_cluster = sapply(1:nrow(trait_mat), function(i) {
    cl <- which(sapply(ap_func@clusters, function(c) rownames(trait_mat)[i] %in% c))
    if (length(cl) > 0) cl else NA
  })
)
write.table(func_clusters, args$out_functional, sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("Functional clusters: %d clusters\n", length(ap_func@clusters)))

# ---- Identity clustering (hierarchical on OG distances) ----
hc_identity <- hclust(og_dist, method = "ward.D2")
# Cut tree at height that gives comparable number of clusters
n_cl <- max(2, length(ap_func@clusters))
id_clusters <- cutree(hc_identity, k = n_cl)
identity_df <- data.frame(
  genome = names(id_clusters),
  identity_cluster = as.integer(id_clusters)
)
write.table(identity_df, args$out_identity, sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("Identity clusters: %d clusters\n", n_cl))

# ---- Divergence metrics ----
# Adjusted Rand Index
if (requireNamespace("mclust", quietly = TRUE)) {
  ari <- mclust::adjustedRandIndex(func_clusters$functional_cluster,
                                    identity_df$identity_cluster)
} else {
  # Manual ARI computation
  contingency <- table(func_clusters$functional_cluster, identity_df$identity_cluster)
  ari <- {
    n <- sum(contingency)
    sum_comb <- sum(choose(contingency, 2))
    row_combs <- sum(choose(rowSums(contingency), 2))
    col_combs <- sum(choose(colSums(contingency), 2))
    expected <- row_combs * col_combs / choose(n, 2)
    max_possible <- (row_combs + col_combs) / 2
    if (max_possible == expected) 0 else (sum_comb - expected) / (max_possible - expected)
  }
}

# Mantel test between distance matrices
mantel_result <- if (requireNamespace("vegan", quietly = TRUE)) {
  vegan::mantel(trait_dist, og_dist, permutations = args$mantel_perms)
} else {
  # Simple correlation-based mantel proxy
  m1 <- as.vector(as.matrix(trait_dist))
  m2 <- as.vector(as.matrix(og_dist))
  list(statistic = cor(m1, m2, use = "complete.obs"), signif = NA)
}

# Robinson-Foulds-like comparison (cluster overlap)
# Identify discordant genome pairs
discordant <- data.frame(genome = character(), func_cl = integer(),
                          id_cl = integer(), discordant = logical())

merged <- merge(func_clusters, identity_df, by = "genome")
for (i in 1:(nrow(merged)-1)) {
  for (j in (i+1):nrow(merged)) {
    same_func <- merged$functional_cluster[i] == merged$functional_cluster[j]
    same_id <- merged$identity_cluster[i] == merged$identity_cluster[j]
    if (same_func != same_id) {
      discordant <- rbind(discordant, data.frame(
        genome = paste(merged$genome[i], merged$genome[j], sep = " vs "),
        func_cl = as.integer(same_func),
        id_cl = as.integer(same_id),
        discordant = TRUE
      ))
    }
  }
}

metrics <- list(
  adjusted_rand_index = round(ari, 4),
  mantel_correlation = round(as.numeric(mantel_result$statistic), 4),
  mantel_significance = ifelse(is.null(mantel_result$signif), NA, round(mantel_result$signif, 4)),
  n_functional_clusters = length(ap_func@clusters),
  n_identity_clusters = n_cl,
  n_discordant_pairs = nrow(discordant),
  interpretation = ifelse(ari > 0.7, "high concordance",
                   ifelse(ari > 0.3, "moderate divergence", "high divergence"))
)

jsonlite::write_json(metrics, args$out_metrics, pretty = TRUE, auto_unbox = TRUE)

# ---- Report ----
report <- sprintf(
"# Clustering Comparison Report

## Summary
- Functional clusters (GFC): %d
- Identity clusters (core genome): %d
- Adjusted Rand Index: %.4f (%s)
- Mantel correlation: %.4f (p = %s)
- Discordant genome pairs: %d

## Interpretation
%s

## Discordant Cases (functional clustering != identity clustering)
Genomes that cluster together functionally but not taxonomically (or vice versa)
are candidates for functional convergence events.

| Genome Pair | Same Functional | Same Identity |
|-------------|-----------------|---------------|
%s
",
  metrics$n_functional_clusters,
  metrics$n_identity_clusters,
  metrics$adjusted_rand_index,
  metrics$interpretation,
  metrics$mantel_correlation,
  ifelse(is.na(metrics$mantel_significance), "NA", sprintf("%.4f", metrics$mantel_significance)),
  metrics$n_discordant_pairs,
  ifelse(metrics$adjusted_rand_index > 0.7,
         "The functional and identity clusterings are highly concordant. Convergence events, if present, affect only a small fraction of the genome.",
         ifelse(metrics$adjusted_rand_index > 0.3,
                "Moderate divergence detected. Some organisms cluster by function rather than taxonomy, suggesting convergent trait acquisition.",
                "High divergence detected. Functional and taxonomic clustering disagree substantially, indicating widespread functional convergence or horizontal gene transfer.")),
  ifelse(nrow(discordant) > 0,
         paste(sprintf("| %s | %s | %s |", discordant$genome,
                        ifelse(discordant$func_cl, "YES", "NO"),
                        ifelse(discordant$id_cl, "YES", "NO")), collapse = "\n"),
         "| (none) |  |  |")
)

writeLines(report, args$out_report)
cat(sprintf("Clustering comparison complete. ARI = %.4f\n", metrics$adjusted_rand_index))
