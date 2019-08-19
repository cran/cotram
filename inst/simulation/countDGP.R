library("cotram")
library("MASS")

# settings
set.seed(29)

Nsim <- 100

## coefficients
b0 <- 1.2
b1 <- 0.8
theta <- 3

## transformation function dgp
pY <- ppois(0:100, lambda = exp(b0))
h.lo <- mlt:::.Logistic()[["q"]](pY) # logit link
h.cll <- mlt:::.MinExtrVal()[["q"]](pY) # cloglog link
h.ll <- mlt:::.MaxExtrVal()[["q"]](pY) # loglog link
h.pr <- mlt:::.Normal()[["q"]](pY) # probit link

# data-generating processes
dgp <- function(n = 1000){
  x <- runif(n, min = 0, max = 1)
  log.mu <- b0 + b1 * x
  
  ## Poisson - DGP
  y.p <- rpois(n = n, lambda = exp(log.mu))
  
  ## negbin - DGP
  y.nb <- rnbinom(n = n, mu = exp(log.mu), size = theta)
  
  ## logit - DGP
  h.lo_m <- matrix(h.lo, nrow = length(h.lo), ncol = length(x))
  p.lo <-  mlt:::.Logistic()[["p"]](t(h.lo_m) - b1 * x) # inv. logit link
  y.lo <- max.col(-(p.lo - runif(n))^2) - 1L
  
  mlo <- cotram(y.lo ~ x, data = data.frame(y.lo = y.lo, x = x), 
                method = "logit", fixed = c("x" = b1))
  
  ## cloglog - DGP
  h.cll_m <- matrix(h.cll, nrow = length(h.cll), ncol = length(x)) 
  p.cll <-  mlt:::.MinExtrVal()[["p"]](t(h.cll_m) - b1 * x) # inv. cloglog link
  y.cll <-  max.col(-(p.cll - runif(n))^2) - 1L
  
  mcll <- cotram(y.cll ~ x, data = data.frame(y.cll = y.cll, x = x), 
                 method = "cloglog", fixed = c("x" = b1))
  
  ## loglog - DGP
  h.ll_m <- matrix(h.ll, nrow = length(h.ll), ncol = length(x)) 
  p.ll <-  mlt:::.MaxExtrVal()[["p"]](t(h.ll_m) - b1 * x) # inv. loglog link
  y.ll <-  max.col(-(p.ll - runif(n))^2) - 1L
  
  mll <- cotram(y.ll ~ x, data = data.frame(y.ll = y.ll, x = x), 
                method = "loglog", fixed = c("x" = b1))
  
  ## probit - DGP
  h.pr_m <- matrix(h.pr, nrow = length(h.pr), ncol = length(x)) 
  p.pr <-  mlt:::.Normal()[["p"]](t(h.pr_m) - b1 * x) # inv. probit link
  y.pr <- max.col(-(p.pr - runif(n))^2) - 1L
  
  mpr <- cotram(y.pr ~ x, data = data.frame(y.pr = y.pr, x = x), 
                method = "probit", fixed = c("x" = b1))
  
  ## data frame
  ret <- data.frame(x = x, mu = exp(log.mu),
                    y.p = y.p, y.nb = y.nb, y.lo = y.lo,
                    y.cll = y.cll, y.ll = y.ll, y.pr = y.pr)
  
  attr(ret, "mlo") <- mlo
  attr(ret, "mcll") <- mcll
  attr(ret, "mll") <- mll
  attr(ret, "mpr") <- mpr
  return(ret)
}

# model out-of-sample log-likelihood
logLikFUN <- function(model, lhs, newdata){ 
  
  if("cotram" %in% class(model)){
    val <- logLik(model, newdata = newdata)
    
  }else if(grepl("Negative Binomial", family(model)$family)){
    val <- sum(dnbinom(newdata[,lhs],
                       mu = predict(model, newdata = newdata, type = "response"),
                       size = model$theta, log = TRUE))
    
  }else if(family(model)$family %in% "poisson"){
    val <- sum(dpois(newdata[,lhs],
                     lambda = predict(model, newdata = newdata, type = "response"),
                     log = TRUE))
  }
  return(val)}

# out-of-sample log-likelihood of generating process
logLik_refFUN <- function(dgp, newdata){
  if(dgp == "lo"){
    mlo <- attr(newdata, "mlo")
    val <- logLik(mlo, newdata = newdata, parm = coef(as.mlt(mlo), fixed = TRUE))
    
  }else if(dgp == "cll"){
    mcll <- attr(newdata, "mcll")
    val <- logLik(mcll, newdata = newdata, parm = coef(as.mlt(mcll), fixed = TRUE))
    
  }else if(dgp == "ll"){
    mll <- attr(newdata, "mll")
    val <- logLik(mll, newdata = newdata, parm = coef(as.mlt(mll), fixed = TRUE))
    
  }else if(dgp == "pr"){
    mpr <- attr(newdata, "mpr")
    val <- logLik(mpr, newdata = newdata, parm = coef(as.mlt(mpr), fixed = TRUE))
    
  }else if(dgp == "nb"){
    val <- sum(dnbinom(newdata$y.nb, mu = newdata$mu, size = theta, log = TRUE))
    
  }else if(dgp == "p"){
    val <- sum(dpois(newdata$y.p, lambda = newdata$mu, log = TRUE))
  }
  return(val)
}

# set-up
dgps <- c("p", "nb", "lo", "cll", "ll", "pr")
ret <- vector(mode = "list", length = length(dgps))
names(ret) <- dgps
for (d in dgps){
  ret[[d]] <- matrix(NA, nrow = Nsim, ncol = length(dgps),
                     dimnames = list(1:Nsim, paste0("m", dgps)))
}

# run Nsim-times
for (i in c(1:Nsim)){
  print(i)
  
  df <- dgp()
  
  # split into training and validation data-set
  trainID <- sample(1:nrow(df), size = 0.25 * nrow(df))
  df_train <- df[trainID,]
  df_test <- df[-trainID,]
  
  # for each data-generating process
  for (d in dgps){
    
    # model formula
    lhs <- paste0("y.", d)
    fm <-  as.formula(paste(lhs, "~ x"))
    
    # count regression models
    mp <- glm(fm, data = df_train, family = poisson(link = "log"))
    mnb <- glm.nb(fm, data = df_train, link = "log")
    
    mlo <- cotram(fm, data = df_train, method = "logit")
    mcll <- cotram(fm, data = df_train, method = "cloglog")
    mll <- cotram(fm, data = df_train, method = "loglog")
    mpr <- cotram(fm, data = df_train, method = "probit")
    
    fit <- list(mp = mp, mnb = mnb, mlo = mlo, mcll = mcll, mll = mll, mpr = mpr)
    
    # out-of-sample log-likelihood of models
    logLiks.d <- sapply(fit, logLikFUN, lhs = lhs, newdata = df_test)
    
    # centered out-of-sample log-likelihood
    logLiks_diff.d <- logLiks.d - logLik_refFUN(newdata = df_test, dgp = d)
    
    colnames(ret[[d]]) <- names(fit)
    ret[[d]][i, ] <- logLiks_diff.d
  }
}

for (d in names(ret)){
  ret[[d]] <- as.data.frame(ret[[d]])
  ret[[d]]$dgp <- d
}

logLiks <- ret

warnings()
