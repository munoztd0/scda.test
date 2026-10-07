test_that("tsflab01", {
  skip_if_not_installed("envsetup")

  expect_snapshot_file(write_test_rtf_for("tsflab01.R"), "tsflab01.rtf")
})
