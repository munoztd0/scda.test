test_that("tsfvit01abpart1of4", {
  skip_if_not_installed("envsetup")

  expect_snapshot_file(write_test_rtf_for("tsfvit01ab.R", part_num = 1, total_parts = 4), "tsfvit01abpart1of4.rtf")
})

# test_that("tsfvit01abpart2of4", {
#   skip_if_not_installed("envsetup")

#   expect_snapshot_file(write_test_rtf_for("tsfvit01ab.R", part_num = 2, total_parts = 4), "tsfvit01abpart2of4.rtf")
# })

# test_that("tsfvit01abpart3of4", {
#   skip_if_not_installed("envsetup")

#   expect_snapshot_file(write_test_rtf_for("tsfvit01ab.R", part_num = 3, total_parts = 4), "tsfvit01abpart3of4.rtf")
# })

# test_that("tsfvit01abpart4of4", {
#   skip_if_not_installed("envsetup")

#   expect_snapshot_file(write_test_rtf_for("tsfvit01ab.R", part_num = 4, total_parts = 4), "tsfvit01abpart4of4.rtf")
# })
