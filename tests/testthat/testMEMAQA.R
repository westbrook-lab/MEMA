library(MEMA)
context("QA Functions")

test_that("CV works correctly",{
  expect_equal(CV(c(5,5,5,NA)),0)
  expect_equal(CV(c(0,5,10, NA)),1)
})
