context("HERC")

test_that("regular HERC", {
    portfolio <- HERC(prices)
    expect_true(sum(portfolio$w) - 1 < 1e-6)
})
