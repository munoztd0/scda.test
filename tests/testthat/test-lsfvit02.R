test_that("lsfvit02part1of3", {
  skip_if_not_installed("envsetup")

  expect_snapshot_file(write_test_rtf_for("lsfvit02.R", part_num = 1, total_parts = 3), "lsfvit02part1of3.rtf")
})

# test_that("lsfvit02part2of3", {
#   skip_if_not_installed("envsetup")

#   expect_snapshot_file(write_test_rtf_for("lsfvit02.R", part_num = 2, total_parts = 3), "lsfvit02part2of3.rtf")
# })

# test_that("lsfvit02part3of3", {
#   skip_if_not_installed("envsetup")

#   expect_snapshot_file(write_test_rtf_for("lsfvit02.R", part_num = 3, total_parts = 3), "lsfvit02part3of3.rtf")
# })
