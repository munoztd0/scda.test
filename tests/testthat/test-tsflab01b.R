test_that("tsflab01b", {
  skip_if_not_installed("envsetup")

  expect_snapshot_file(write_test_rtf_for("tsflab01b.R"), "tsflab01b.rtf")
})
