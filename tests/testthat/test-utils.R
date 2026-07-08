# Tests for internal helpers in utils.R.
# create_folds() splits a data frame into contiguous folds; the last fold
# absorbs the remainder when num_folds does not divide the row count evenly.
# (Two-column frames are used so row subsetting keeps its data.frame shape.)

test_that("create_folds splits rows into equal contiguous folds", {
  d <- data.frame(id = 1:10, x = 11:20)
  f <- create_folds(d, 5)
  expect_length(f, 5)
  expect_equal(vapply(f, nrow, integer(1)), rep(2L, 5))
  expect_equal(sort(unlist(lapply(f, `[[`, "id"))), 1:10)   # every row kept exactly once (no dups/missing)
})

test_that("create_folds puts the remainder in the last fold", {
  d <- data.frame(id = 1:11, x = 21:31)
  f <- create_folds(d, 5)
  expect_equal(vapply(f, nrow, integer(1)), c(2L, 2L, 2L, 2L, 3L))
  expect_equal(sort(unlist(lapply(f, `[[`, "id"))), 1:11)
})

test_that("create_folds with a single fold returns all rows", {
  d <- data.frame(id = 1:7, x = 8:14)
  f <- create_folds(d, 1)
  expect_length(f, 1)
  expect_equal(nrow(f[[1]]), 7)
  expect_equal(sort(f[[1]]$id), 1:7)
})
