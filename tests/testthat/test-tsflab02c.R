test_that("tsflab02c", {
  skip_if_not_installed("envsetup")

  expect_snapshot_file(write_test_rtf_for("tsflab02c.R"), "tsflab02c.rtf")
})
