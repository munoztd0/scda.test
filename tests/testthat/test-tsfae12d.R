test_that("tsfae12d", {
  skip_if_not_installed("envsetup")

  expect_snapshot_file(write_test_rtf_for("tsfae12d.R"), "tsfae12dpart1of2.rtf")
})
