library(testthat)
library(scNoiseFilter)

cat('=== Testing scNoiseFilter ===\n\n')

# Test 1: validate_species
cat('1. Testing validate_species...\n')
expect_equal(scNoiseFilter:::validate_species('mouse'), 'mouse')
expect_equal(scNoiseFilter:::validate_species('rat'), 'rat')
expect_equal(scNoiseFilter:::validate_species('human'), 'human')
expect_equal(scNoiseFilter:::validate_species('mmu'), 'mouse')
expect_equal(scNoiseFilter:::validate_species('hsa'), 'human')
cat('   PASSED\n')

# Test 2: get_noise_patterns
cat('2. Testing get_noise_patterns...\n')
mouse_patterns <- get_noise_patterns('mouse')
expect_true('riken' %in% names(mouse_patterns))
expect_true('predicted' %in% names(mouse_patterns))

rat_patterns <- get_noise_patterns('rat')
expect_true('riken' %in% names(rat_patterns))
expect_true('predicted' %in% names(rat_patterns))

human_patterns <- get_noise_patterns('human')
expect_false('riken' %in% names(human_patterns))
expect_false('predicted' %in% names(human_patterns))
cat('   PASSED\n')

# Test 3: filter_noise_genes mouse
cat('3. Testing filter_noise_genes (mouse)...\n')
genes <- c('Gapdh', 'Actb', 'Gm12345', '1700001C02Rik', 'AI597479', 'Rpl13a', 'mt-Co1', 'Hba-a1')
result <- filter_noise_genes(genes, species = 'mouse', verbose = FALSE)
expect_true('Gapdh' %in% result$filtered_genes)
expect_true('Actb' %in% result$filtered_genes)
expect_false('Gm12345' %in% result$filtered_genes)
expect_false('1700001C02Rik' %in% result$filtered_genes)
expect_false('AI597479' %in% result$filtered_genes)
expect_false('Rpl13a' %in% result$filtered_genes)
expect_false('mt-Co1' %in% result$filtered_genes)
expect_false('Hba-a1' %in% result$filtered_genes)
cat('   PASSED\n')

# Test 4: filter_noise_genes human
cat('4. Testing filter_noise_genes (human)...\n')
genes <- c('GAPDH', 'ACTB', 'RPL13A', 'MT-CO1', 'HBA1')
result <- filter_noise_genes(genes, species = 'human', verbose = FALSE)
expect_true('GAPDH' %in% result$filtered_genes)
expect_true('ACTB' %in% result$filtered_genes)
expect_false('RPL13A' %in% result$filtered_genes)
expect_false('MT-CO1' %in% result$filtered_genes)
expect_false('HBA1' %in% result$filtered_genes)
cat('   PASSED\n')

# Test 5: Grik1 NOT filtered
cat('5. Testing Grik1 is NOT filtered...\n')
genes <- c('Grik1', 'Gapdh')
result <- filter_noise_genes(genes, species = 'mouse', verbose = FALSE)
expect_true('Grik1' %in% result$filtered_genes)
cat('   PASSED\n')

# Test 6: Verbose output
cat('6. Testing verbose output...\n')
genes <- c('Gapdh', 'Gm12345', 'Rpl13a')
result <- filter_noise_genes(genes, species = 'mouse', verbose = TRUE)
cat('   PASSED\n')

# Test 7: Rat patterns
cat('7. Testing rat patterns...\n')
genes <- c('Gapdh', 'RGD12345', '1200002C02Rik', 'Rpl13a')
result <- filter_noise_genes(genes, species = 'rat', verbose = FALSE)
expect_true('Gapdh' %in% result$filtered_genes)
expect_false('RGD12345' %in% result$filtered_genes)
expect_false('1200002C02Rik' %in% result$filtered_genes)
cat('   PASSED\n')

# Test 8: Report generation
cat('8. Testing get_noise_report...\n')
genes <- c('Gapdh', 'Gm12345', 'Rpl13a', 'mt-Co1')
result <- filter_noise_genes(genes, species = 'mouse', verbose = FALSE)
report <- get_noise_report(result, format = 'list')
expect_true(!is.null(report$breakdown))
cat('   PASSED\n')

cat('\n=== All tests passed! ===\n')
