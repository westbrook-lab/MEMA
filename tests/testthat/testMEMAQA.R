library(MEMA)


context("QA Functions")

test_that("CV works correctly",{
  expect_equal(CV(c(5,5,5,NA)),0)
  expect_equal(CV(c(0,5,10, NA)),1)
})


context("MEMA Test Functions")

# test_that("renameSpottedWells works correctly",{
#   expect_equal(renameSpottedWells(data.table::data.table(Well=1:10),well = "A01"),data.table::data.table(Spot=1:10, Well= "A01"))
# })

# test_that("rotateMetadata works correctly",{
#   expect_equal(rotateMetadata(MEMA::spotMetadata),MEMA::spotMetadata180)
# })

# test_that("readSpotMetadata works correctly",{
#   expect_equal(smd, MEMA::readSpotMetadata(system.file("extdata", "20150403_LI8V002_16ECM_28pin.gal", package = "MEMA")))
# })

# test_that("melt8Well works correctly",{
#   expect_equal(melt8Well(MEMA::popDataRaw),MEMA::popData)
# })

