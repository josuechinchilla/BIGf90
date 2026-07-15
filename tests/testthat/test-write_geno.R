# Tests for write_geno(): PLINK .ped -> BLUPf90 .geno.
# write_geno() returns the written IDs invisibly and writes the .geno to disk,
# so we test the return value and the round-tripped file. Each test builds a
# tiny ped with writeLines instead of shipping fixtures.

# helper: write toy lines to a temp file and return its path
tmp_lines <- function(x, ext = ".ped") {
  p <- tempfile(fileext = ext)
  writeLines(x, p)
  p
}

test_that("--recode12 ped (two alleles/locus) collapses to a dosage .geno", {
  # loci per row:  (1 1)->0  (1 2)->1  (2 2)->2   and   (2 2)->2 (1 2)->1 (1 1)->0
  ped <- tmp_lines(c("fam1 A1 0 0 1 -9 1 1 1 2 2 2",
                     "fam2 A2 0 0 2 -9 2 2 1 2 1 1"))
  out <- tempfile(fileext = ".geno")
  ids <- suppressMessages(write_geno(ped = ped, file = out))

  expect_identical(ids, c("A1", "A2"))              # returned IDs, in file order
  expect_identical(readLines(out), c("A1 012", "A2 210"))
})

test_that("IDs are padded so genotypes start at the same column (fixed format)", {
  ped <- tmp_lines(c("f short 0 0 1 -9 1 1 1 2 2 2",
                     "f longername 0 0 1 -9 2 2 1 2 1 1"))
  out <- tempfile(fileext = ".geno")
  suppressMessages(write_geno(ped = ped, file = out))
  lines  <- readLines(out)
  starts <- regexpr("[012]{3}$", lines)               # where the dosages begin
  expect_identical(starts[1], starts[2])              # same column on every row
  expect_identical(sub("^.* ", "", lines), c("012", "210"))  # correct dosages
})

test_that("missing alleles (0 0) and a custom missing code become 5", {
  ped <- tmp_lines("fam1 A1 0 0 1 -9 1 1 0 0 2 2")  # middle locus is 0 0 -> missing
  out <- tempfile(fileext = ".geno")
  suppressMessages(write_geno(ped = ped, file = out))
  expect_identical(readLines(out), "A1 052")

  # a stray custom missing allele code (here "9") is also honoured
  ped9 <- tmp_lines("fam1 A1 0 0 1 -9 1 1 9 9 2 2")
  out9 <- tempfile(fileext = ".geno")
  suppressMessages(write_geno(ped = ped9, file = out9, missing = 9))
  expect_identical(readLines(out9), "A1 052")
})

test_that("count_allele switches which homozygote is 0 vs 2", {
  ped <- tmp_lines("fam1 A1 0 0 1 -9 1 1 1 2 2 2")
  out <- tempfile(fileext = ".geno")
  suppressMessages(write_geno(ped = ped, file = out, count_allele = 1))
  expect_identical(readLines(out), "A1 210")        # counts allele 1 instead of 2
})

test_that("alleles_per_locus = 1 accepts an already-dosage matrix", {
  mat <- tmp_lines(c("IND1 0 1 2", "IND2 2 1 0"))
  out <- tempfile(fileext = ".geno")
  suppressMessages(write_geno(ped = mat, file = out, id_col = 1, n_lead_cols = 1,
                             alleles_per_locus = 1))
  expect_identical(readLines(out), c("IND1 012", "IND2 210"))
})

test_that("map row count is validated and map_out is written", {
  ped  <- tmp_lines(c("fam1 A1 0 0 1 -9 1 1 1 2 2 2",
                      "fam2 A2 0 0 2 -9 2 2 1 2 1 1"))   # 3 loci
  good <- tmp_lines(c("1 rs1 0 100", "1 rs2 0 200", "1 rs3 0 300"), ext = ".map")
  bad  <- tmp_lines(c("1 rs1 0 100", "1 rs2 0 200"), ext = ".map")

  out  <- tempfile(fileext = ".geno")
  mout <- tempfile(fileext = ".map")
  suppressMessages(write_geno(ped = ped, file = out, map = good, map_out = mout))
  expect_identical(readLines(mout), c("rs1 1 100", "rs2 1 200", "rs3 1 300"))

  expect_error(
    suppressMessages(write_geno(ped = ped, file = tempfile(fileext = ".geno"), map = bad)),
    "Locus count mismatch"
  )
})

test_that("non-recode12 allele codes and odd column counts are rejected", {
  letters_ped <- tmp_lines("fam1 A1 0 0 1 -9 A A C C")     # letters, not 1/2
  expect_error(
    suppressMessages(write_geno(ped = letters_ped, file = tempfile(fileext = ".geno"))),
    "recode12"
  )
  odd_ped <- tmp_lines("fam1 A1 0 0 1 -9 1 1 2")           # odd genotype-column count
  expect_error(
    suppressMessages(write_geno(ped = odd_ped, file = tempfile(fileext = ".geno"))),
    "two allele columns per locus"
  )
})

test_that("informative errors on missing/existing files", {
  expect_error(write_geno(file = tempfile()), "'ped'")            # no ped
  expect_error(write_geno(ped = tempfile()), "'file'")           # no file

  ped <- tmp_lines("fam1 A1 0 0 1 -9 1 1 1 2 2 2")
  out <- tmp_lines("already here", ext = ".geno")                # file exists
  expect_error(suppressMessages(write_geno(ped = ped, file = out)), "already exists")
  ids <- suppressMessages(write_geno(ped = ped, file = out, overwrite = TRUE))
  expect_identical(ids, "A1")
})
