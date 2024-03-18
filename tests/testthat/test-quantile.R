test_that("qged", {
    dist <- "ged"
    set.seed(101)
    r <- rdist(dist, 100000, mu = 0.1, sigma = 1, shape = 5)
    median_r = median(r)
    median_q <- qged(0.5, mu = 0.1, sigma = 1, shape = 5)
    p <- pged(median_q, mu = 0.1, sigma = 1, shape = 5)
    expect_equal(median_r, median_q, tolerance = 0.05)
    expect_equal(p, 0.5)
})

test_that("qstd", {
    dist <- "std"
    set.seed(101)
    r <- rdist(dist, 100000, mu = 0.1, sigma = 1, shape = 4.2)
    median_r = median(r)
    median_q <- qstd(0.5, mu = 0.1, sigma = 1, shape = 4.2)
    p <- pstd(median_q, mu = 0.1, sigma = 1, shape = 4.2)
    expect_equal(median_r, median_q, tolerance = 0.06)
    expect_equal(p, 0.5)
})


test_that("qghst", {
    dist <- "ghst"
    set.seed(101)
    r <- rdist(dist, 10000, mu = 0.01, sigma = 1, skew = -20, shape = 5)
    median_r = median(r)
    median_q <- qghst(0.5, mu = 0.01, sigma = 1, skew = -20, shape = 5)
    p <- pghst(median_q, mu = 0.01, sigma = 1, skew = -20, shape = 5)
    expect_equal(median_r, median_q, tolerance = 0.01)
    expect_equal(p, 0.5, tolerance = 0.0001)
})

test_that("qnig", {
    dist <- "nig"
    set.seed(101)
    r <- rdist(dist, 10000, mu = 0.1, sigma = 1, skew = -0.8, shape = 4)
    median_r = median(r)
    median_q <- qnig(0.5, mu = 0.1, sigma = 1, skew = -0.8, shape = 4)
    p <- pnig(median_q, mu = 0.1, sigma = 1, skew = -0.8, shape = 4)
    expect_equal(median_r, median_q, tolerance = 0.1)
    expect_equal(p, 0.5, tolerance = 0.0001)
})

test_that("qjsu", {
    dist <- "jsu"
    set.seed(101)
    r <- rdist(dist, 100000, mu = 0.1, sigma = 1, skew = 1.5, shape = 4)
    median_r = median(r)
    median_q <- qjsu(0.5, mu = 0.1, sigma = 1, skew = 1.5, shape = 4)
    p <- pjsu(median_q, mu = 0.1, sigma = 1, skew = 1.5, shape = 4)
    expect_equal(median_r, median_q, tolerance = 0.05)
    expect_equal(p, 0.5)
})

test_that("qsstd", {
    dist <- "sstd"
    set.seed(101)
    r <- rdist(dist, 100000, mu = 0.1, sigma = 2, skew = 1.5, shape = 4)
    median_r = median(r)
    median_q <- qsstd(0.5, mu = 0.1, sigma = 2, skew = 1.5, shape = 4)
    p <- psstd(median_q, mu = 0.1, sigma = 2, skew = 1.5, shape = 4)
    expect_equal(median_r, median_q, tolerance = 0.02)
    expect_equal(p, 0.5)
})

test_that("qsnorm", {
    dist <- "snorm"
    set.seed(101)
    r <- rdist(dist, 100000, mu = 0.1, sigma = 2, skew = 1.5, shape = 4)
    median_r = median(r)
    median_q <- qsnorm(0.5, mu = 0.1, sigma = 2, skew = 1.5)
    p <- psnorm(median_q, mu = 0.1, sigma = 2, skew = 1.5)
    expect_equal(median_r, median_q, tolerance = 0.02)
    expect_equal(p, 0.5)
})

test_that("qsged", {
    dist <- "sged"
    set.seed(101)
    r <- rdist(dist, 100000, mu = 0.1, sigma = 2, skew = 1.5, shape = 2)
    median_r = median(r)
    median_q <- qsged(0.5, mu = 0.1, sigma = 2, skew = 1.5, shape = 2)
    p <- psged(median_q, mu = 0.1, sigma = 2, skew = 1.5, shape = 2)
    expect_equal(median_r, median_q, tolerance = 0.02)
    expect_equal(p, 0.5)
})
