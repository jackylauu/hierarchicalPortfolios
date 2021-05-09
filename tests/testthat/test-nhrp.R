context("NHRP")
data(prices)

test_that("regular NHRP", {
    portfolio <- NHRP(prices)
    expect_true(sum(portfolio$w) == 1)
    expect_true(class(portfolio$tree) == 'hclust')
})

test_that("NHRP with weight constraints", {
    portfolio <- NHRP(prices, w_min=0.2, w_max=0.4)
    expect_true(min(portfolio$w) >= 0.2 - 1e-5)
    expect_true(max(portfolio$w) <= 0.4 + 1e-5)
})
