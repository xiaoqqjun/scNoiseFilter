#' Generate a summary report of noise genes
#'
#' Creates a formatted report showing the breakdown of noise genes by category.
#' Useful for documentation and QC reporting.
#'
#' @param result Output from filter_noise_genes()
#' @param format Output format: "text", "data.frame", or "list"
#' @return Report in the specified format
#' @export
#' @examples
#' genes <- c("Gapdh", "Actb", "Gm12345", "Rpl13a", "mt-Co1")
#' result <- filter_noise_genes(genes, species = "mouse", verbose = FALSE)
#' report <- get_noise_report(result)
#' print(report)
get_noise_report <- function(result, format = c("text", "data.frame", "list")) {
  format <- match.arg(format)

  if (!is.list(result) || is.null(result$noise_genes)) {
    stop("Invalid input. Expected output from filter_noise_genes()")
  }

  # Build summary data
  categories <- names(result$noise_genes)

  summary_data <- lapply(categories, function(cat) {
    info <- result$noise_genes[[cat]]
    data.frame(
      category = cat,
      name = info$name,
      count = info$count,
      stringsAsFactors = FALSE
    )
  })

  summary_df <- do.call(rbind, summary_data)

  # Add totals row
  total_row <- data.frame(
    category = "TOTAL",
    name = "总计",
    count = result$removed_count,
    stringsAsFactors = FALSE
  )
  summary_df <- rbind(summary_df, total_row)

  if (format == "data.frame") {
    return(summary_df)
  }

  if (format == "list") {
    return(list(
      species = result$species,
      total_input = result$total_genes,
      total_removed = result$removed_count,
      total_retained = result$retained_count,
      retention_rate = round(result$retained_count / result$total_genes * 100, 2),
      breakdown = summary_df
    ))
  }

  # Text format
  lines <- c(
    "=" %>% paste(rep(., 50), collapse = ""),
    sprintf("Noise Gene Filtering Report - %s", toupper(result$species)),
    "=" %>% paste(rep(., 50), collapse = ""),
    "",
    sprintf("Input genes: %d", result$total_genes),
    sprintf("Removed genes: %d", result$removed_count),
    sprintf("Retained genes: %d", result$retained_count),
    sprintf("Retention rate: %.1f%%", result$retained_count / result$total_genes * 100),
    "",
    "Breakdown by category:",
    "-" %>% paste(rep(., 40), collapse = "")
  )

  for (cat in categories) {
    info <- result$noise_genes[[cat]]
    lines <- c(lines, sprintf("  %s: %d", info$name, info$count))
  }

  lines <- c(lines, "")

  paste(lines, collapse = "\n")
}

#' Get gene list for a specific noise category
#'
#' Extracts the gene list from a filter result for a specific category.
#'
#' @param result Output from filter_noise_genes()
#' @param category Category name: "riken", "predicted", "digits", "ribosomal", "mito", "hemoglobin", "high_umi"
#' @return Character vector of gene names, or NULL if category not found
#' @export
#' @examples
#' genes <- c("Gapdh", "Gm12345", "1700001C02Rik", "Rpl13a")
#' result <- filter_noise_genes(genes, species = "mouse", verbose = FALSE)
#' riken_genes <- get_category_genes(result, "riken")
get_category_genes <- function(result, category) {
  if (!is.list(result) || is.null(result$noise_genes)) {
    stop("Invalid input. Expected output from filter_noise_genes()")
  }

  if (!(category %in% names(result$noise_genes))) {
    warning(sprintf("Category '%s' not found in result", category))
    return(NULL)
  }

  result$noise_genes[[category]]$genes
}
