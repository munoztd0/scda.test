test_that("tsfae11a", {
  skip_if_not_installed("envsetup")

  expect_snapshot_file(write_test_rtf_for("tsfae11a.R"), "tsfae11a.rtf")
})
