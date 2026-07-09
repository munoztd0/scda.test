test_that("tsfae20d", {
  skip_if_not_installed("envsetup")

  expect_snapshot_file(write_test_rtf_for("tsfae20d.R"), "tsfae20dpart1of3.rtf")
})
