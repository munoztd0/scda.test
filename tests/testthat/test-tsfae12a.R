test_that("tsfae12a", {
  skip_if_not_installed("envsetup")

  expect_snapshot_file(write_test_rtf_for("tsfae12a.R"), "tsfae12a.rtf")
})
