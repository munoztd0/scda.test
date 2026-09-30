test_that("tsfae10d", {
  skip_if_not_installed("envsetup")

  expect_snapshot_file(write_test_rtf_for("tsfae10d.R"), "tsfae10d.rtf")
})
