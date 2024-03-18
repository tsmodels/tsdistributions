test_that("norm", {
 dist <- "norm"
 set.seed(1)
 r <- rdist(dist, 1000, mu = 0.1, sigma = 0.3)
 spec <- distribution_modelspec(r, dist)
 mod <- estimate(spec)
 pars <- coef(mod)
 rx <- (r - pars[1])/pars[2]
 test_llh <- sum(-log(ddist(dist, rx, mu = 0, sigma = 1)/pars[2]))
 expect_equal(mod$loglik, test_llh)
})

test_that("std", {
 dist <- "std"
 set.seed(1)
 r <- rdist(dist, 1000, mu = 0.1, sigma = 0.3, shape = 5)
 spec <- distribution_modelspec(r, dist)
 mod <- estimate(spec)
 pars <- coef(mod)
 rx <- (r - pars[1])/pars[2]
 test_llh <- sum(-log(ddist(dist, rx, mu = 0, sigma = 1, shape = pars[3])/pars[2]))
 expect_equal(mod$loglik, test_llh)
})

test_that("ged", {
 dist <- "ged"
 set.seed(1)
 r <- rdist(dist, 1000, mu = 0.1, sigma = 0.3, shape = 2)
 spec <- distribution_modelspec(r, dist)
 mod <- estimate(spec)
 pars <- coef(mod)
 rx <- (r - pars[1])/pars[2]
 test_llh <- sum(-log(ddist(dist, rx, mu = 0, sigma = 1, shape = pars[3])/pars[2]))
 expect_equal(mod$loglik, test_llh)
})

test_that("ged", {
 dist <- "ged"
 set.seed(1)
 r <- rdist(dist, 1000, mu = 0.1, sigma = 0.3, shape = 2)
 spec <- distribution_modelspec(r, dist)
 mod <- estimate(spec)
 pars <- coef(mod)
 rx <- (r - pars[1])/pars[2]
 test_llh <- sum(-log(ddist(dist, rx, mu = 0, sigma = 1, shape = pars[3])/pars[2]))
 expect_equal(mod$loglik, test_llh)
})

test_that("snorm", {
 dist <- "snorm"
 set.seed(1)
 r <- rdist(dist, 1000, mu = 0.1, sigma = 0.3, skew = 1)
 spec <- distribution_modelspec(r, dist)
 mod <- estimate(spec)
 pars <- coef(mod)
 rx <- (r - pars[1])/pars[2]
 test_llh <- sum(-log(ddist(dist, rx, mu = 0, sigma = 1, skew = pars[3])/pars[2]))
 expect_equal(mod$loglik, test_llh)
})

test_that("sged", {
 dist <- "sged"
 set.seed(1)
 r <- rdist(dist, 1000, mu = 0.1, sigma = 0.3, skew = 1, shape = 1)
 spec <- distribution_modelspec(r, dist)
 mod <- estimate(spec)
 pars <- coef(mod)
 rx <- (r - pars[1])/pars[2]
 test_llh <- sum(-log(ddist(dist, rx, mu = 0, sigma = 1, skew = pars[3], shape = pars[4])/pars[2]))
 expect_equal(mod$loglik, test_llh)
})

test_that("sstd", {
 dist <- "sstd"
 set.seed(1)
 r <- rdist(dist, 1000, mu = 0.1, sigma = 0.3, skew = 1, shape = 5)
 spec <- distribution_modelspec(r, dist)
 mod <- estimate(spec)
 pars <- coef(mod)
 rx <- (r - pars[1])/pars[2]
 test_llh <- sum(-log(ddist(dist, rx, mu = 0, sigma = 1, skew = pars[3], shape = pars[4])/pars[2]))
 expect_equal(mod$loglik, test_llh)
})

test_that("nig", {
 dist <- "nig"
 set.seed(1)
 r <- rdist(dist, 1000, mu = 0.1, sigma = 0.3, skew = -0.3, shape = 2)
 spec <- distribution_modelspec(r, dist)
 mod <- estimate(spec)
 pars <- coef(mod)
 rx <- (r - pars[1])/pars[2]
 test_llh <- sum(-log(ddist(dist, rx, mu = 0, sigma = 1, skew = pars[3], shape = pars[4])/pars[2]))
 expect_equal(mod$loglik, test_llh)
})

test_that("jsu", {
 dist <- "jsu"
 set.seed(1)
 r <- rdist(dist, 1000, mu = 0.1, sigma = 0.3, skew = -10, shape = 1)
 spec <- distribution_modelspec(r, dist)
 mod <- estimate(spec)
 pars <- coef(mod)
 rx <- (r - pars[1])/pars[2]
 test_llh <- sum(-log(ddist(dist, rx, mu = 0, sigma = 1, skew = pars[3], shape = pars[4])/pars[2]))
 expect_equal(round(mod$loglik,1), round(test_llh,1))
})

test_that("ghst", {
 dist <- "ghst"
 set.seed(1)
 r <- rdist(dist, 1000, mu = 0.1, sigma = 0.3, skew = -40, shape = 5)
 spec <- distribution_modelspec(r, dist)
 mod <- estimate(spec)
 pars <- coef(mod)
 rx <- (r - pars[1])/pars[2]
 test_llh <- sum(-log(ddist(dist, rx, mu = 0, sigma = 1, skew = pars[3], shape = pars[4])/pars[2]))
 expect_equal(mod$loglik, test_llh)
})

test_that("gh", {
 dist <- "gh"
 set.seed(1)
 r <- rdist(dist, 1000, mu = 0.1, sigma = 0.3, skew = -0.5, shape = 5, lambda = 1)
 spec <- distribution_modelspec(r, dist)
 mod <- estimate(spec)
 pars <- coef(mod)
 rx <- (r - pars[1])/pars[2]
 test_llh <- sum(-log(ddist(dist, rx, mu = 0, sigma = 1, skew = pars[3], shape = pars[4], lambda = pars[5])/pars[2]))
 expect_equal(mod$loglik, test_llh)
})
