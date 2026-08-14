test_that("lsfae07b", {
  skip_if_not_installed("envsetup")

  expect_snapshot_file(write_test_rtf_for("lsfae07b.R"), "lsfae07b.rtf")
})
