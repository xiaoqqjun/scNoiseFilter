#' Get noise gene patterns for a specific species
#'
#' Returns species-specific patterns for identifying noise genes.
#' Mouse and rat include Riken and predicted genes (Gm/RGD).
#' Human starts from digit genes only.
#'
#' @param species Species name: "mouse", "rat", or "human"
#' @return Named list of noise gene patterns with metadata
#' @export
#' @examples
#' patterns <- get_noise_patterns("mouse")
#' patterns <- get_noise_patterns("human")
get_noise_patterns <- function(species) {
  species <- validate_species(species)

  # Base patterns common to all species
  patterns <- list()

  # Digit genes pattern (5+ consecutive digits, excluding Riken and predicted)
  if (species %in% c("mouse", "rat")) {
    # Exclude Riken and Gm/RGD from digit genes
    # Riken pattern: Rik followed by digits at end
    patterns$digits <- list(
      pattern = "^(?=.*\\d{5,})(?!.*Rik[0-9]*$)(?!^Gm\\d+)(?!^RGD\\d+).+$",
      name = "数字基因",
      use_perl = TRUE,
      ignore_case = FALSE
    )
  } else {
    # Human: exclude Rik and GM (though they shouldn't exist in human)
    patterns$digits <- list(
      pattern = "^(?=.*\\d{5,})(?!.*Rik[0-9]*$)(?!^GM\\d+).+$",
      name = "数字基因",
      use_perl = TRUE,
      ignore_case = FALSE
    )
  }

  # Ribosomal genes - must be followed by digit
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
    # Human: comprehensive MT patterns including MTND, MTCO, MTATP, MTCYB
    patterns$mito <- list(
      pattern = "^MT-|^MTND|^MTCO|^MTATP|^MTCYB",
      name = "线粒体基因",
      use_perl = FALSE,
      ignore_case = FALSE
    )
  } else {
    # Mouse/Rat: support both lowercase mt- and uppercase MT-
    patterns$mito <- list(
      pattern = "^mt-|^MT-",
      name = "线粒体基因",
      use_perl = FALSE,
      ignore_case = FALSE
    )
  }

  # Hemoglobin genes
  if (species == "human") {
    # Human: comprehensive HB patterns
    patterns$hemoglobin <- list(
      pattern = "^HBA|^HBB|^HBQ|^HBD|^HBE|^HBG|^HBM|^HBZ",
      name = "血红蛋白基因",
      use_perl = FALSE,
      ignore_case = FALSE
    )
  } else {
    # Mouse/Rat: Hba, Hbb, Hbq
    patterns$hemoglobin <- list(
      pattern = "^Hba|^Hbb|^Hbq",
      name = "血红蛋白基因",
      use_perl = FALSE,
      ignore_case = FALSE
    )
  }

  # Riken genes (mouse and rat only)
  # Pattern: ends with Rik followed by optional digits
  if (species %in% c("mouse", "rat")) {
    patterns$riken <- list(
      pattern = "Rik[0-9]*$",
      name = "Riken基因",
      use_perl = FALSE,
      ignore_case = FALSE
    )
  }

  # Predicted genes (mouse: Gm, rat: RGD)
  if (species == "mouse") {
    patterns$predicted <- list(
      pattern = "^Gm\\d+$",
      name = "预测基因",
      use_perl = FALSE,
      ignore_case = FALSE
    )
  } else if (species == "rat") {
    patterns$predicted <- list(
      pattern = "^RGD\\d+$",
      name = "预测基因",
      use_perl = FALSE,
      ignore_case = FALSE
    )
  }
  # Human has no predicted genes category

  attr(patterns, "species") <- species
  patterns
}
