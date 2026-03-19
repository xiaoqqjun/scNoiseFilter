#' Filter noise genes from a Seurat object or gene list
#'
#' Identifies and filters noise genes based on species-specific patterns.
#' Returns detailed information about filtered genes in each category.
#'
#' @param seurat_obj A Seurat object (preferred) - will auto-extract genes and UMI matrix
#' @param genes Character vector of gene names (optional if seurat_obj provided)
#' @param species Species name: "mouse", "rat", or "human" (required)
#' @param umi_matrix UMI count matrix (optional, for sampling scenarios)
#' @param sample_rate Sampling rate for high-umi calculation (default: 1.0 = 100%, 0.1 = 10%)
#' @param high_umi_threshold Median percent threshold for high-abundance genes (default: 1, meaning 1%)
#' @param filter_riken Filter Riken genes (default: TRUE)
#' @param filter_predicted Filter predicted genes Gm/RGD (default: TRUE)
#' @param filter_digits Filter 5-digit genes (default: TRUE)
#' @param filter_ribosomal Filter ribosomal genes (default: TRUE)
#' @param filter_mito Filter mitochondrial genes (default: TRUE)
#' @param filter_hemoglobin Filter hemoglobin genes (default: TRUE)
#' @param filter_high_umi Filter high-abundance genes (default: TRUE)
#' @param verbose Print progress messages (default: TRUE)
#'
#' @return A list containing:
#' \itemize{
#'   \item total_genes: Total number of input genes
#'   \item noise_genes: Named list of filtered genes by category
#'   \item filtered_genes: Gene vector after filtering
#'   \item removed_genes: All removed genes combined
#'   \item removed_count: Total number of removed genes
#'   \item retained_count: Number of retained genes
#'   \item species: Species used for filtering
#'   \item sampled_cells: Number of cells used (after sampling)
#' }
#'
#' @export
#'
#' @examples
#' # Simple usage with Seurat object
#' result <- filter_noise_genes(seurat_obj, species = "mouse")
#'
#' # With sampling for large datasets
#' result <- filter_noise_genes(seurat_obj, species = "mouse", sample_rate = 0.1)
#'
#' # With custom threshold
#' result <- filter_noise_genes(seurat_obj, species = "mouse", high_umi_threshold = 0.5)
#'
#' # Gene list only (no high-umi filtering)
#' result <- filter_noise_genes(genes = c("Gapdh", "Gm12345"), species = "mouse")
#'
filter_noise_genes <- function(seurat_obj = NULL,
                                genes = NULL,
                                species,
                                umi_matrix = NULL,
                                sample_rate = 1.0,
                                high_umi_threshold = 1,
                                filter_riken = TRUE,
                                filter_predicted = TRUE,
                                filter_digits = TRUE,
                                filter_ribosomal = TRUE,
                                filter_mito = TRUE,
                                filter_hemoglobin = TRUE,
                                filter_high_umi = TRUE,
                                verbose = TRUE) {

  # Validate species first
  if (missing(species)) {
    stop("'species' is required. Use 'mouse', 'rat', or 'human'")
  }
  species <- validate_species(species)

  # Validate sample_rate
  if (!is.numeric(sample_rate) || sample_rate <= 0 || sample_rate > 1) {
    stop("'sample_rate' must be a number between 0 and 1")
  }

  # Track sampled cells
  sampled_cells <- 0

  # Extract genes and UMI matrix from Seurat object if provided
  if (!is.null(seurat_obj)) {
    if (!requireNamespace("Seurat", quietly = TRUE)) {
      stop("Package 'Seurat' is required when using seurat_obj parameter")
    }

    # Auto-extract genes
    if (is.null(genes)) {
      genes <- rownames(seurat_obj)
    }

    # Auto-extract UMI matrix if not provided and high-umi filtering is enabled
    if (is.null(umi_matrix) && filter_high_umi) {
      umi_matrix <- get_umi_matrix(seurat_obj)
      if (verbose) {
        message("Auto-extracted UMI matrix from Seurat object")
      }
    }
  }

  # Validate genes
  if (is.null(genes) || !is.character(genes)) {
    stop("'genes' must be provided (either directly or via seurat_obj)")
  }

  # Apply sampling if sample_rate < 1 and umi_matrix is available
  if (!is.null(umi_matrix) && sample_rate < 1.0) {
    total_cells <- ncol(umi_matrix)
    n_sample <- max(100, round(total_cells * sample_rate))
    n_sample <- min(n_sample, total_cells)

    set.seed(42)  # Fixed seed for reproducibility
    sample_cols <- sample(colnames(umi_matrix), n_sample)
    umi_matrix <- umi_matrix[, sample_cols]
    sampled_cells <- n_sample

    if (verbose) {
      message(sprintf("Sampling: using %d of %d cells (%.0f%%)", n_sample, total_cells, sample_rate * 100))
    }
  } else if (!is.null(umi_matrix)) {
    sampled_cells <- ncol(umi_matrix)
  }

  # Get species-specific patterns
  patterns <- get_noise_patterns(species)

  if (verbose) {
    message(sprintf("Filtering noise genes for species: %s", species))
    message(sprintf("Total input genes: %d", length(genes)))
  }

  # Initialize results
  noise_genes <- list()
  removed_set <- character(0)

  # Helper function to process a filter category
  process_filter <- function(pattern_info, filter_enabled, filter_name) {
    if (!filter_enabled || is.null(pattern_info)) {
      return(NULL)
    }

    result <- filter_genes_by_pattern(
      genes = genes,
      pattern = pattern_info$pattern,
      name = pattern_info$name,
      use_perl = pattern_info$use_perl,
      ignore_case = pattern_info$ignore_case
    )

    if (result$count > 0) {
      if (verbose) {
        message(sprintf("  %s: %d genes", result$name, result$count))
      }
      result
    } else {
      NULL
    }
  }

  # Apply filters based on species and user settings

  # Riken (mouse and rat only)
  if (species %in% c("mouse", "rat") && filter_riken) {
    result <- process_filter(patterns$riken, TRUE, "riken")
    if (!is.null(result)) {
      noise_genes$riken <- result
      removed_set <- union(removed_set, result$genes)
    }
  }

  # Predicted genes (Gm for mouse, RGD for rat)
  if (species %in% c("mouse", "rat") && filter_predicted) {
    result <- process_filter(patterns$predicted, TRUE, "predicted")
    if (!is.null(result)) {
      noise_genes$predicted <- result
      removed_set <- union(removed_set, result$genes)
    }
  }

  # Digit genes (all species)
  if (filter_digits) {
    result <- process_filter(patterns$digits, TRUE, "digits")
    if (!is.null(result)) {
      noise_genes$digits <- result
      removed_set <- union(removed_set, result$genes)
    }
  }

  # Ribosomal genes (all species)
  if (filter_ribosomal) {
    result <- process_filter(patterns$ribosomal, TRUE, "ribosomal")
    if (!is.null(result)) {
      noise_genes$ribosomal <- result
      removed_set <- union(removed_set, result$genes)
    }
  }

  # Mitochondrial genes (all species)
  if (filter_mito) {
    result <- process_filter(patterns$mito, TRUE, "mito")
    if (!is.null(result)) {
      noise_genes$mito <- result
      removed_set <- union(removed_set, result$genes)
    }
  }

  # Hemoglobin genes (all species)
  if (filter_hemoglobin) {
    result <- process_filter(patterns$hemoglobin, TRUE, "hemoglobin")
    if (!is.null(result)) {
      noise_genes$hemoglobin <- result
      removed_set <- union(removed_set, result$genes)
    }
  }

  # High abundance genes using median percent method
  if (filter_high_umi && !is.null(umi_matrix)) {

    # Get common genes between input and matrix
    matrix_genes <- rownames(umi_matrix)
    common_genes <- intersect(matrix_genes, genes)

    if (length(common_genes) > 0) {
      # Subset matrix to common genes
      umi_subset <- umi_matrix[common_genes, , drop = FALSE]

      # Exclude genes already matched by other patterns
      genes_for_high_umi <- setdiff(common_genes, removed_set)

      if (length(genes_for_high_umi) > 0) {
        umi_for_high_umi <- umi_subset[genes_for_high_umi, , drop = FALSE]

        # Calculate percentage per cell: each gene / cell total * 100
        cell_totals <- Matrix::colSums(umi_for_high_umi)
        cell_totals[cell_totals == 0] <- 1

        gene_percent_per_cell <- Matrix::t(Matrix::t(umi_for_high_umi) / cell_totals) * 100

        # Calculate median percent across all cells for each gene
        gene_median_percent <- apply(gene_percent_per_cell, 1, median)

        # Find genes exceeding threshold
        high_umi_genes <- names(gene_median_percent[gene_median_percent >= high_umi_threshold])

        if (length(high_umi_genes) > 0) {
          high_umi_info <- data.frame(
            gene = high_umi_genes,
            median_percent = gene_median_percent[high_umi_genes],
            stringsAsFactors = FALSE
          )
          high_umi_info <- high_umi_info[order(-high_umi_info$median_percent), ]

          noise_genes$high_umi <- list(
            name = "高丰度基因",
            pattern = sprintf("Median percent >= %.2f%%", high_umi_threshold),
            genes = high_umi_info$gene,
            count = nrow(high_umi_info),
            details = high_umi_info
          )
          removed_set <- union(removed_set, high_umi_info$gene)

          if (verbose) {
            message(sprintf("  高丰度基因: %d genes (median percent >= %.2f%%)", nrow(high_umi_info), high_umi_threshold))
          }
        } else {
          if (verbose) {
            message(sprintf("  高丰度基因: 0 genes (no genes with median percent >= %.2f%%)", high_umi_threshold))
          }
        }
      }
    }
  } else if (filter_high_umi && is.null(umi_matrix)) {
    if (verbose) {
      message("  高丰度基因: skipped (no UMI matrix available)")
    }
  }

  # Calculate final results
  filtered_genes <- setdiff(genes, removed_set)

  if (verbose) {
    message(sprintf("\nSummary:"))
    message(sprintf("  Removed: %d genes", length(removed_set)))
    message(sprintf("  Retained: %d genes", length(filtered_genes)))
  }

  # Return comprehensive result
  list(
    total_genes = length(genes),
    noise_genes = noise_genes,
    filtered_genes = filtered_genes,
    removed_genes = removed_set,
    removed_count = length(removed_set),
    retained_count = length(filtered_genes),
    species = species,
    sampled_cells = sampled_cells,
    params = list(
      sample_rate = sample_rate,
      high_umi_threshold = high_umi_threshold,
      filter_riken = filter_riken,
      filter_predicted = filter_predicted,
      filter_digits = filter_digits,
      filter_ribosomal = filter_ribosomal,
      filter_mito = filter_mito,
      filter_hemoglobin = filter_hemoglobin,
      filter_high_umi = filter_high_umi
    )
  )
}
