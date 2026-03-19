# Tests for scNoiseFilter

test_that("validate_species works correctly", {
  expect_equal(scNoiseFilter:::validate_species("mouse"), "mouse")
  expect_equal(scNoiseFilter:::validate_species("MOUSE"), "mouse")
  expect_equal(scNoiseFilter:::validate_species("Mus musculus"), "mouse")
  expect_equal(scNoiseFilter:::validate_species("mmu"), "mouse")

  expect_equal(scNoiseFilter:::validate_species("rat"), "rat")
  expect_equal(scNoiseFilter:::validate_species("RAT"), "rat")
  expect_equal(scNoiseFilter:::validate_species("rno"), "rat")

  expect_equal(scNoiseFilter:::validate_species("human"), "human")
  expect_equal(scNoiseFilter:::validate_species("HUMAN"), "human")
  expect_equal(scNoiseFilter:::validate_species("hsa"), "human")

  expect_error(scNoiseFilter:::validate_species("invalid"))
})

test_that("get_noise_patterns returns correct patterns for mouse", {
  patterns <- get_noise_patterns("mouse")

  expect_true("riken" %in% names(patterns))
  expect_true("predicted" %in% names(patterns))
  expect_true("digits" %in% names(patterns))
  expect_true("ribosomal" %in% names(patterns))
  expect_true("mito" %in% names(patterns))
  expect_true("hemoglobin" %in% names(patterns))

  expect_equal(patterns$predicted$pattern, "^Gm\\d+$")
})

test_that("get_noise_patterns returns correct patterns for rat", {
  patterns <- get_noise_patterns("rat")

  expect_true("riken" %in% names(patterns))
  expect_true("predicted" %in% names(patterns))
  expect_equal(patterns$predicted$pattern, "^RGD\\d+$")
})

test_that("get_noise_patterns returns correct patterns for human", {
  patterns <- get_noise_patterns("human")

  expect_false("riken" %in% names(patterns))
  expect_false("predicted" %in% names(patterns))
  expect_true("digits" %in% names(patterns))
  expect_true("ribosomal" %in% names(patterns))
  expect_equal(patterns$ribosomal$pattern, "^RPL|^RPS")
  expect_equal(patterns$mito$pattern, "^MT-")
})

test_that("filter_noise_genes filters mouse genes correctly", {
  genes <- c(
    "Gapdh", "Actb",  # Normal genes
    "Gm12345", "Gm67890",  # Predicted genes
    "1700001C02Rik", "1200002C02Rik",  # Riken genes
    "AI597479", "AA986860",  # Digit genes
    "Rpl13a", "Rps18",  # Ribosomal
    "mt-Co1", "mt-Nd1",  # Mitochondrial
    "Hba-a1", "Hbb-bt"  # Hemoglobin
  )

  result <- filter_noise_genes(genes, species = "mouse", verbose = FALSE)

  expect_equal(result$total_genes, 14)
  expect_true(result$removed_count > 0)
  expect_true("Gapdh" %in% result$filtered_genes)
  expect_true("Actb" %in% result$filtered_genes)
  expect_false("Gm12345" %in% result$filtered_genes)
  expect_false("1700001C02Rik" %in% result$filtered_genes)
})

test_that("filter_noise_genes filters human genes correctly", {
  genes <- c(
    "GAPDH", "ACTB",  # Normal genes
    "AI597479",  # Digit gene
    "RPL13A", "RPS18",  # Ribosomal
    "MT-CO1", "MT-ND1",  # Mitochondrial
    "HBA1", "HBB"  # Hemoglobin
  )

  result <- filter_noise_genes(genes, species = "human", verbose = FALSE)

  expect_true("GAPDH" %in% result$filtered_genes)
  expect_true("ACTB" %in% result$filtered_genes)
  expect_false("RPL13A" %in% result$filtered_genes)
  expect_false("MT-CO1" %in% result$filtered_genes)
})

test_that("filter_noise_genes respects filter flags", {
  genes <- c("Gm12345", "Rpl13a", "mt-Co1", "Gapdh")

  # All filters on
  result1 <- filter_noise_genes(genes, species = "mouse",
                                 filter_predicted = TRUE,
                                 filter_ribosomal = TRUE,
                                 filter_mito = TRUE,
                                 verbose = FALSE)
  expect_equal(length(result1$filtered_genes), 1)  # Only Gapdh

  # Turn off predicted filter
  result2 <- filter_noise_genes(genes, species = "mouse",
                                 filter_predicted = FALSE,
                                 filter_ribosomal = TRUE,
                                 filter_mito = TRUE,
                                 verbose = FALSE)
  expect_true("Gm12345" %in% result2$filtered_genes)
})

test_that("filter_noise_genes handles empty input", {
  result <- filter_noise_genes(character(0), species = "mouse", verbose = FALSE)
  expect_equal(result$total_genes, 0)
  expect_equal(result$removed_count, 0)
})

test_that("filter_noise_genes requires species parameter", {
  genes <- c("Gapdh", "Actb")
  expect_error(filter_noise_genes(genes), "species")
})

test_that("digit genes pattern excludes Riken and Gm", {
  # Grik1 should NOT be filtered as digit gene (it's not 5 digits)
  genes <- c("Grik1", "Gapdh")
  result <- filter_noise_genes(genes, species = "mouse", verbose = FALSE)

  expect_true("Grik1" %in% result$filtered_genes)

  # AI597479 should be filtered (5 consecutive digits)
  genes2 <- c("AI597479", "Gapdh")
  result2 <- filter_noise_genes(genes2, species = "mouse", verbose = FALSE)

  expect_false("AI597479" %in% result2$filtered_genes)
})
