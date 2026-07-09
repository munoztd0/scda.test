test_that("tsfae11b", {
  skip_if_not_installed("envsetup")

  expect_snapshot_file(write_test_rtf_for("tsfae11b.R"), "tsfae11b.rtf")
})
