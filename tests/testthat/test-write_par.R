# Tests for write_par(): assembling RENUMF90 parameter (.par) files.
# write_par() returns the written lines invisibly, so we test the return value
# (and, in one case, the round-tripped file) instead of eyeballing text.

test_that("single-trait genomic model produces the expected .par lines", {
  tmp <- tempfile(fileext = ".par")
  res <- write_par(
    file     = tmp,
    datafile = "sealice_match.txt",
    traits   = 4,
    residual_variance = 0.157,
    effects  = list(mean   = list(col = 2, type = "cross", class = "numer"),
                    tank   = list(col = 3, type = "cross", class = "numer"),
                    animal = list(col = 1, type = "cross", class = "alpha")),
    random        = 1,
    pedigree_file = "sealice.ped",
    snp_file      = "sealice_match.geno",
    missing_value = -999,
    options = c("callrate 0.9", "saveGInverse")
  )

  expected <- c(
    "DATAFILE", "sealice_match.txt",
    "TRAITS", "4",
    "FIELDS_PASSED TO OUTPUT", "",
    "WEIGHT(S)", "",
    "RESIDUAL_VARIANCE", "0.157",
    "EFFECT #mean", "2 cross numer",
    "EFFECT #tank", "3 cross numer",
    "EFFECT #animal", "1 cross alpha",
    "RANDOM", "animal",
    "FILE", "sealice.ped",
    "FILE_POS", "1 2 3 0 0",
    "SNP_FILE", "sealice_match.geno",
    "OPTION missing -999",
    "OPTION callrate 0.9",
    "OPTION saveGInverse"
  )
  expect_identical(res, expected)             # returned lines
  expect_identical(readLines(tmp), expected)  # round-trips to disk
})

test_that("random matches by column, independent of the effect label", {
  res <- write_par(
    file = tempfile(fileext = ".par"), datafile = "d.txt", traits = 1,
    residual_variance = 1,
    effects = list(mu = list(col = 2, class = "numer"),
                   id = list(col = 1, class = "alpha")),   # label 'id', not 'animal'
    random = 1, pedigree_file = "ped.txt"
  )
  i <- match("EFFECT #id", res)               # the col-1 effect, whatever its label
  expect_identical(res[i + 1], "1 cross alpha")
  expect_identical(res[i + 2], "RANDOM")
  expect_identical(res[i + 3], "animal")
  expect_identical(res[match("FILE", res) + 1], "ped.txt")
})

test_that("a random effect not listed last is moved to the end (with a message)", {
  expect_message(
    res <- write_par(
      file = tempfile(fileext = ".par"), datafile = "d.txt", traits = 1,
      residual_variance = 1,
      effects = list(plant = list(col = 2, class = "numer"),   # random, listed FIRST
                     sky   = list(col = 3, class = "numer"),
                     fish  = list(col = 1, class = "alpha")),
      random = 2, pedigree_file = "ped.txt"
    ),
    "moving random effect"
  )
  expect_identical(res[grepl("^EFFECT", res)],
                   c("EFFECT #sky", "EFFECT #fish", "EFFECT #plant"))  # plant now last
  i <- match("EFFECT #plant", res)
  expect_identical(res[i + 1], "2 cross numer")
  expect_identical(res[i + 2], "RANDOM")      # RANDOM directly follows its effect
})

test_that("'pos' is an alias for 'col'", {
  a <- write_par(file = tempfile(fileext = ".par"), datafile = "d.txt", traits = 1,
                 residual_variance = 1,
                 effects = list(mu = list(pos = 2, class = "numer"),
                                id = list(pos = 1, class = "alpha")),
                 random = 1, pedigree_file = "ped.txt")
  b <- write_par(file = tempfile(fileext = ".par"), datafile = "d.txt", traits = 1,
                 residual_variance = 1,
                 effects = list(mu = list(col = 2, class = "numer"),
                                id = list(col = 1, class = "alpha")),
                 random = 1, pedigree_file = "ped.txt")
  expect_identical(a, b)
})

test_that("'comment' overrides the label on the EFFECT line", {
  res <- write_par(file = tempfile(fileext = ".par"), datafile = "d.txt", traits = 1,
                   residual_variance = 1,
                   effects = list(a = list(col = 1, class = "numer", comment = "my label")))
  expect_true("EFFECT #my label" %in% res)
})

test_that("multi-trait: matrix variances, per-trait columns, OPTIONAL and (CO)VARIANCES", {
  res <- write_par(
    file = tempfile(fileext = ".par"), datafile = "act.csv", traits = c(6, 7),
    residual_variance = matrix(c(1, 0, 0, 1), 2, 2),
    effects = list(cg     = list(col = 5, class = "alpha"),
                   animal = list(col = 1, class = "alpha")),
    random = 1, optional = "pe", pedigree_file = "p.ped",
    snp_file = "g.geno", covariances = matrix(c(1, 0.1, 0.1, 1), 2, 2),
    options = "cat 0 0"
  )
  expect_identical(res[match("TRAITS", res) + 1], "6 7")
  rv <- match("RESIDUAL_VARIANCE", res)
  expect_identical(res[c(rv + 1, rv + 2)], c("1 0", "0 1"))   # one matrix row per line
  expect_true("5 5 cross alpha" %in% res)                     # column repeated per trait
  expect_true("1 1 cross alpha" %in% res)
  expect_identical(res[match("RANDOM", res) + 1], "animal")
  expect_identical(res[match("OPTIONAL", res) + 1], "pe")
  cv <- match("(CO)VARIANCES", res)
  expect_identical(res[c(cv + 1, cv + 2)], c("1 0.1", "0.1 1"))
  expect_true("OPTION cat 0 0" %in% res)
})

test_that("fixed-effects-only model writes no RANDOM/FILE block", {
  res <- write_par(file = tempfile(fileext = ".par"), datafile = "d.txt", traits = 1,
                   residual_variance = 1,
                   effects = list(mu = list(col = 2, class = "numer")))
  expect_false(any(c("RANDOM", "FILE", "FILE_POS") %in% res))
})

test_that("informative errors on bad input", {
  expect_error(                                              # missing column
    write_par(file = tempfile(), datafile = "d", traits = 1, residual_variance = 1,
              effects = list(a = list(class = "numer"))),
    "data-file column"
  )
  expect_error(                                              # random matches nothing
    write_par(file = tempfile(), datafile = "d", traits = 1, residual_variance = 1,
              effects = list(a = list(col = 1, class = "numer")),
              random = 9, pedigree_file = "p"),
    "No effect matches"
  )
  expect_error(                                              # pedigree required with random
    write_par(file = tempfile(), datafile = "d", traits = 1, residual_variance = 1,
              effects = list(a = list(col = 1, class = "alpha")),
              random = 1),
    "pedigree_file"
  )
  expect_error(                                              # bad class
    write_par(file = tempfile(), datafile = "d", traits = 1, residual_variance = 1,
              effects = list(a = list(col = 1, class = "bogus"))),
    "'class'"
  )
})

test_that("refuses to overwrite an existing file unless overwrite = TRUE", {
  tmp <- tempfile(fileext = ".par")
  writeLines("x", tmp)                                       # file now exists
  expect_error(
    write_par(file = tmp, datafile = "d", traits = 1, residual_variance = 1,
              effects = list(a = list(col = 1, class = "numer"))),
    "already exists"
  )
  res <- write_par(file = tmp, datafile = "d", traits = 1, residual_variance = 1,
                   effects = list(a = list(col = 1, class = "numer")), overwrite = TRUE)
  expect_true("DATAFILE" %in% res)
})
