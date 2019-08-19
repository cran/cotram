
library("cotram")

set.seed(25)

yr <- rpois(100, lambda = 2)
yy <- 0:100 / 10

y <- 0:10
pr <- 1:9/10

### PREDICT
layout(matrix(c(1:4), nrow = 2, byrow = TRUE))
## cotram: log_first = FALSE
m1 <- cotram(yr ~ 1, log_first = FALSE, extrapolate = TRUE)

stopifnot(all.equal(predict(m1, newdata = data.frame(1), q = 0, type = "density"), 
                    predict(m1, newdata = data.frame(1), q = 0, type = "distribution")))

stopifnot(!isTRUE(all.equal(predict(m1, newdata = data.frame(1), q = 1, type = "density"), 
                  predict(m1, newdata = data.frame(1), q = 1, type = "distribution"))))

if (FALSE) {
main <- "predict - density: log_first = FALSE"
plot(y, predict(m1, newdata = data.frame(1), q = y , type = "density"), 
     type = "h", main = main)
points(y, predict(m1, newdata = data.frame(1), q = y, type = "density"), pch = 20)
points(y, dpois(y, lambda = 2), col = "blue")
lines(yy, predict(m1, newdata = data.frame(1), q = yy, type = "density",
                  smooth = TRUE), col = "grey")

main <- "predict - density: log_first = FALSE"
plot(y, predict(m1, newdata = data.frame(1), q = y , type = "density"), 
     type = "h", main = main)
points(y, predict(m1, newdata = data.frame(1), q = y, type = "density"), pch = 20)
points(y, dpois(y, lambda = 2), col = "blue")
lines(yy, predict(m1, newdata = data.frame(1), q = yy, type = "density",
                  smooth = TRUE), col = "grey")

main <- "predict - distribution: log_first = FALSE"
plot(y, predict(m1, newdata = data.frame(1), q = y, type = "distribution"),
            type = "s", ylim = c(0, 1), main = main)
lines(yy, predict(m1, newdata = data.frame(1), q = yy, type = "distribution"), col = "grey")
points(y, ppois(y, lambda = 2), col = "blue")
points(predict(m1, newdata = data.frame(1), prob = pr, type = "quantile", 
                smooth = "TRUE"), pr, col = "red")
legend("bottomright", legend = c("ppois", "quantiles"), pch = 1, 
       col = c("blue", "red"))
}

## cotram: log_first = TRUE
m1_lf <- cotram(yr ~ 1, log_first = TRUE, extrapolate = TRUE)

stopifnot(all.equal(predict(m1_lf, newdata = data.frame(1), q = 0, type = "density"), 
          predict(m1_lf, newdata = data.frame(1), q = 0, type = "distribution")))

stopifnot(!isTRUE(all.equal(predict(m1_lf, newdata = data.frame(1), q = 1, type = "density"), 
                  predict(m1_lf, newdata = data.frame(1), q = 1, type = "distribution"))))

if (FALSE) {
main <- "predict - density: log_first = TRUE"
plot(y, predict(m1_lf, newdata = data.frame(1), q = y, type = "density"), 
     type = "h", main = main)
points(y, predict(m1_lf, newdata = data.frame(1), q = y, type = "density"), pch = 20)
lines(yy, predict(m1_lf, newdata = data.frame(1), q = yy, type = "density",
                  smooth = TRUE))

main <- "predict - distribution: log_first = TRUE"
plot(y, predict(m1_lf, newdata = data.frame(1), q = y, type = "distribution"),
     type = "s", ylim = c(0, 1), main = main)
points(y, ppois(y, lambda = 2), col = "blue")
lines(yy, predict(m1_lf, newdata = data.frame(1), q = yy, type = "distribution"))
points(predict(m1_lf, newdata = data.frame(1), 
               smooth = TRUE, prob = pr, type = "quantile"), pr, col = "red")
legend("bottomright", legend = c("ppois", "quantiles"), pch = 1, 
       col = c("blue", "red"))


### PLOT
layout(matrix(c(1:4), nrow = 2, byrow = TRUE))

main <- "plot: log_first = FALSE"
plot(m1, newdata = data.frame(1), type = "distribution", col = "black", main = main)
plot(m1, newdata = data.frame(1), type = "distribution", col = "grey",
     smooth = TRUE, add = TRUE)

plot(m1, newdata = data.frame(1), type = "density", col = "black", main = main)
plot(m1, newdata = data.frame(1), type = "density", col = "grey",
     smooth = TRUE, add = TRUE)

main <- "plot: log_first = TRUE"
plot(m1_lf, newdata = data.frame(1), type = "distribution", col = "black", main = main)
plot(m1_lf, newdata = data.frame(1), type = "distribution", col = "grey",
     smooth = TRUE, add = TRUE)

plot(m1_lf, newdata = data.frame(1), type = "density", col = "black", main = main)
plot(m1_lf, newdata = data.frame(1), type = "density", col = "grey",
     smooth = TRUE, add = TRUE)

}
