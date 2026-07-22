# Input-validation tests for run_pregs() / run_postgs().
# These exercise the argument checks only (no BLUPf90 binaries are run), so they
# are safe on CI / CRAN.

test_that("run_pregs / run_postgs stop when the output directory does not exist", {
  nodir <- file.path(tempdir(), "does_not_exist_xyz123")
  expect_error(run_pregs(path_2_execs = tempdir(), output_files_dir = nodir),  "does not exist")
  expect_error(run_postgs(path_2_execs = tempdir(), output_files_dir = nodir), "does not exist")
})

test_that("run_pregs / run_postgs stop when the executable is missing", {
  d <- tempfile(); dir.create(d)
  expect_error(run_pregs(path_2_execs = d, input_files_dir = d, output_files_dir = d),  "Executable not found")
  expect_error(run_postgs(path_2_execs = d, input_files_dir = d, output_files_dir = d), "Executable not found")
})

test_that("run_pregs / run_postgs stop when the parameter file is missing", {
  execs <- tempfile(); dir.create(execs)
  exe <- if (.Platform$OS.type == "windows") c("preGSf90.exe", "postGSf90.exe") else c("preGSf90", "postGSf90")
  file.create(file.path(execs, exe))                 # fake executables so the exec check passes
  d <- tempfile(); dir.create(d)                     # empty run dir: no renf90.par
  expect_error(run_pregs(path_2_execs = execs, input_files_dir = d, output_files_dir = d),  "Parameter file not found")
  expect_error(run_postgs(path_2_execs = execs, input_files_dir = d, output_files_dir = d), "Parameter file not found")
})

test_that("run_predf validates its arguments", {
  expect_error(run_predf(path_2_execs = tempdir(), output_files_dir = tempdir()), "Provide the genotype file")
  nodir <- file.path(tempdir(), "does_not_exist_xyz123")
  expect_error(run_predf(path_2_execs = tempdir(), output_files_dir = nodir, snp_file = "x.geno"), "does not exist")
  d <- tempfile(); dir.create(d)
  expect_error(run_predf(path_2_execs = d, input_files_dir = d, output_files_dir = d, snp_file = "x.geno"), "Executable not found")
  execs <- tempfile(); dir.create(execs)
  exe <- if (.Platform$OS.type == "windows") "predf90.exe" else "predf90"
  file.create(file.path(execs, exe))                 # fake executable so the exec check passes
  d2 <- tempfile(); dir.create(d2)
  expect_error(run_predf(path_2_execs = execs, input_files_dir = d2, output_files_dir = d2, snp_file = "missing.geno"),
               "Genotype file not found")
})
