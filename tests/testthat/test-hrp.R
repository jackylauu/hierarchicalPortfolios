context("HRP")
data(prices)

test_that("regular HRP", {
    portfolio <- HRP(prices)
    expect_true(sum(portfolio$w) == 1)
    expect_true(class(portfolio$tree) == 'hclust')
})

test_that("HRP with weight constraints", {
    portfolio <- HRP(prices, w_min=0.2, w_max=0.4)
    expect_true(min(portfolio$w) >= 0.2 - 1e-6)
    expect_true(max(portfolio$w) <= 0.4)
})
