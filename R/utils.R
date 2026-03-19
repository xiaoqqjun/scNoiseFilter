#' Filter genes by pattern
#'
#' @param genes Character vector of gene names
#' @param pattern Regular expression pattern
#' @param name Category name for the filter
#' @param use_perl Use Perl-compatible regex
#' @param ignore_case Case sensitivity
#' @return List with name, pattern, genes, and count
#' @keywords internal
filter_genes_by_pattern <- function(genes, pattern, name, use_perl = FALSE, ignore_case = FALSE) {
  matched <- grep(pattern, genes, value = TRUE, ignore.case = ignore_case, perl = use_perl)
  list(
    name = name,
    pattern = pattern,
    genes = matched,
    count = length(matched)
  )
}

#' Validate species parameter
#'
#' @param species Species name
#' @return Normalized species name
#' @keywords internal
validate_species <- function(species) {
  species <- tolower(species)

  species_map <- c(
    "mouse" = "mouse",
    "mus musculus" = "mouse",
    "mmu" = "mouse",
    "rat" = "rat",
    "rattus norvegicus" = "rat",
    "rno" = "rat",
    "human" = "human",
    "homo sapiens" = "human",
    "hsa" = "human"
  )

  if (!(species %in% names(species_map))) {
    stop("Invalid species. Must be one of: mouse, rat, human")
  }

  unname(species_map[species])
}

#' Get noise patterns for a species
#'
#' @param species Species name
#' @return List of pattern info
#' @keywords internal
get_noise_patterns <- function(species) {
  species <- validate_species(species)

  patterns <- list()

  # Digit genes (5+ consecutive digits)
  if (species %in% c("mouse", "rat")) {
    patterns$digits <- list(
      pattern = "^(?=.*\\d{5,})(?!.*Rik[0-9]*$)(?!^Gm\\d+)(?!^RGD\\d+).+$",
      name = "数字基因",
      use_perl = TRUE,
      ignore_case = FALSE
    )
    patterns$riken <- list(
      pattern = "Rik[0-9]*$",
      name = "Riken基因",
      use_perl = FALSE,
      ignore_case = FALSE
    )
    if (species == "mouse") {
      patterns$predicted <- list(
        pattern = "^Gm\\d+$",
        name = "Gm基因",
        use_perl = FALSE,
        ignore_case = FALSE
      )
    } else {
      patterns$predicted <- list(
        pattern = "^RGD\\d+$",
        name = "RGD基因",
        use_perl = FALSE,
        ignore_case = FALSE
      )
    }
  } else {
    # Human
    patterns$digits <- list(
      pattern = "^(?=.*\\d{5,}).+$",
      name = "数字基因",
      use_perl = TRUE,
      ignore_case = FALSE
    )
  }

  # Ribosomal genes
  if (species == "human") {
    patterns$ribosomal <- list(
      pattern = "^RPL[0-9]|^RPS[0-9]",
      name = "核糖体基因",
      use_perl = FALSE,
      ignore_case = FALSE
    )
  } else {
    patterns$ribosomal <- list(
      pattern = "^Rpl[0-9]|^Rps[0-9]",
      name = "核糖体基因",
      use_perl = FALSE,
      ignore_case = FALSE
    )
  }

  # Mitochondrial genes
  if (species == "human") {
    patterns$mito <- list(
      pattern = "^MT-|^MTND|^MTCO|^MTATP|^MTCYB",
      name = "线粒体基因",
      use_perl = FALSE,
      ignore_case = FALSE
    )
  } else {
    patterns$mito <- list(
      pattern = "^mt-|^MT-",
      name = "线粒体基因",
      use_perl = FALSE,
      ignore_case = FALSE
    )
  }

  # Hemoglobin genes
  if (species == "human") {
    patterns$hemoglobin <- list(
      pattern = "^HBA|^HBB|^HBQ|^HBD|^HBE|^HBG|^HBM|^HBZ",
      name = "血红蛋白基因",
      use_perl = FALSE,
      ignore_case = FALSE
    )
  } else {
    patterns$hemoglobin <- list(
      pattern = "^Hba|^Hbb|^Hbq",
      name = "血红蛋白基因",
      use_perl = FALSE,
      ignore_case = FALSE
    )
  }

  patterns
}

#' Get UMI count matrix from Seurat object
#'
#' Extracts the UMI count matrix from a Seurat object, compatible with both v4 and v5.
#'
#' @param seurat_obj A Seurat object
#' @return A sparse matrix of UMI counts (genes x cells)
#' @export
#' @examples
#' \dontrun{
#' umi_matrix <- get_umi_matrix(seurat_obj)
#' }
get_umi_matrix <- function(seurat_obj) {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' is required for this function")
  }

  if (is.null(seurat_obj@assays$RNA)) {
    stop("Seurat object does not have 'RNA' assay")
  }

  # Seurat v5: Assay5 class
  if (inherits(seurat_obj@assays$RNA, "Assay5")) {
    Seurat::LayerData(seurat_obj, layer = "counts")
  } else {
    # Seurat v4: Assay class
    Seurat::GetAssayData(seurat_obj, slot = "counts")
  }
}
