test_that("tsids03", {
  skip_if_not_installed("envsetup")

  expect_snapshot_file(write_test_rtf_for("tsids03.R"), "tsids03.rtf")
})
