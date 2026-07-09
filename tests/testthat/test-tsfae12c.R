test_that("tsfae12c", {
  skip_if_not_installed("envsetup")

  expect_snapshot_file(write_test_rtf_for("tsfae12c.R"), "tsfae12c.rtf")
})
