test_that("lsfae07a", {
  skip_if_not_installed("envsetup")

  expect_snapshot_file(write_test_rtf_for("lsfae07a.R"), "lsfae07a.rtf")
})
