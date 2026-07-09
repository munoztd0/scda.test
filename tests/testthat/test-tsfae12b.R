test_that("tsfae12b", {
  skip_if_not_installed("envsetup")

  expect_snapshot_file(write_test_rtf_for("tsfae12b.R"), "tsfae12b.rtf")
})
