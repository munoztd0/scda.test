test_that("tsflab01a", {
  skip_if_not_installed("envsetup")

  expect_snapshot_file(write_test_rtf_for("tsflab01a.R"), "tsflab01a.rtf")
})
