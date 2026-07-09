test_that("tsfvit01apart1of4", {
  skip_if_not_installed("envsetup")

  expect_snapshot_file(write_test_rtf_for("tsfvit01a.R", part_num = 1, total_parts = 4), "tsfvit01apart1of4.rtf")
})

# test_that("tsfvit01apart2of4", {
#   skip_if_not_installed("envsetup")

#   expect_snapshot_file(write_test_rtf_for("tsfvit01a.R", part_num = 2, total_parts = 4), "tsfvit01apart2of4.rtf")
# })

# test_that("tsfvit01apart3of4", {
#   skip_if_not_installed("envsetup")

#   expect_snapshot_file(write_test_rtf_for("tsfvit01a.R", part_num = 3, total_parts = 4), "tsfvit01apart3of4.rtf")
# })

# test_that("tsfvit01apart4of4", {
#   skip_if_not_installed("envsetup")

#   expect_snapshot_file(write_test_rtf_for("tsfvit01a.R", part_num = 4, total_parts = 4), "tsfvit01apart4of4.rtf")
# })
