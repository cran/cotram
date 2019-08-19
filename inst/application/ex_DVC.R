library("cotram")
library("MASS")
library("multcomp")

# loading data
if (!file.exists("analysis/DVC.rda")) {
  download.file("https://zenodo.org/record/17179/files/DVC.tgz", "DVC.tgz")
  untar("DVC.tgz", file = "analysis/DVC.rda")
}
load("analysis/DVC.rda")

# setups
loc <- Sys.setlocale("LC_ALL", "en_US.UTF-8")
rm(loc)

# data frame
df <- data.frame(wild = as.vector(obs[,,"wild"]),
                 day = seq(from = start, to = end,by = "30 min")[1:prod(dim(obs)[1:2])])
df$weekday <- factor(format(df$day, "%A"),
                     levels = c("Monday", "Tuesday", "Wednesday",
                                "Thursday", "Friday", "Saturday", "Sunday"))
df$time <- as.numeric(difftime(df$day, start, unit = "days"))
df$year <- factor(format(df$day, "%Y"))
df$Datum <- factor(format(df$day, "%Y-%m-%d"))
df$daytime <- cut(as.numeric(difftime(df$day, trunc(df$day, "day"), unit = "hours")), 
                  breaks = c(-Inf, 12, Inf), labels = c("am", "pm"))

# extract sunrise, sunset
w <- weekdays[, c("Datum", "SAofficial", "SUofficial", "ArbZeitFaktor")]
w$Datum <- factor(format(w$Datum, "%Y-%m-%d"))

df <- merge(df, w, by = "Datum")
df$weekday[df$ArbZeitFaktor == 0] <- "Sunday"

a <- as.numeric(with(df, difftime(day, SAofficial, unit = "mins")))
df$sunrise <- cut(a, breaks = c(-Inf, -120, -15, 120, Inf), 
                  labels = c("night", "pre.sunrise", "post.sunrise", "day")) 

a <- as.numeric(with(df, difftime(day, SUofficial, unit = "mins")))
df$sunset <- cut(a, breaks = c(-Inf, -120, -15, 120, Inf), 
                 labels = c("day", "pre.sunset", "post.sunset", "night")) 

df$hours <- with(df, interaction(sunrise, sunset))[, drop = TRUE]
levels(df$hours) <- c("night", "pre.sunrise", "post.sunrise", "day",
                      "pre.sunset", "post.sunset", "night")
df$hours <- relevel(df$hours, "day")
df$daytime <- with(df, interaction(hours, daytime)[, drop = TRUE])

# harmonics
sfm <- function(timevar, freq = 365, S = 10) {
  S <- 1:S * 2
  paste("I(", rep(c("sin", "cos"), length(S)), "(",
        rep(S, rep(2, length(S))), "* pi *", timevar, "/", freq, "))",
        collapse = "+")
}

Xtime <- model.matrix(as.formula(paste("~", sfm("time"))), data = df)[,-1]
Xtime <- as.data.frame(Xtime)
colnames(Xtime) <- paste0("tvar", 1:ncol(Xtime))
df <- cbind(df, Xtime)

# training and validation data-set
trainID <-  as.character(2002:2009)
df_train <- subset(df, year %in% trainID)
df_train$year <- droplevels(df_train$year)
df_test <- subset(df, !year %in% trainID)

# last value carried forward
df_test$year <- "2009"
df_test$year <- factor(df_test$year, levels = levels(df_train$year))

# model formula
rhs <- paste("year + daytime * weekday +", 
             paste("daytime:", colnames(Xtime), "", collapse = "+"))

fm <- as.formula(paste("wild ~", rhs))

# count regression models
mp <- glm(fm, data = df_train, family = poisson(link = "log"))

mnb <- glm.nb(fm, data = df_train, link = "log")

mcll <- cotram(fm, data = df_train, method = "cloglog", prob = .99)

# out-of-sample log-likelihoods
logLik_mp <- sum(dpois(df_test$wild,
                       lambda = predict(mp, newdata = df_test, type = "response"),
                       log = TRUE))

logLik_mnb <- sum(dnbinom(df_test$wild,
                          mu = predict(mnb, newdata = df_test, type = "response"),
                          size = mnb$theta, log = TRUE))

logLik_mcll <- logLik(mcll, newdata = df_test)

logLiks <- c(p = logLik_mp, nb = logLik_mnb, cll = logLik_mcll)

# estimated risk and 95% simultaneous confidence bands
ci.FUN <- function(mod){
  ret <- vector(mode = "list", length = nlevels(df$daytime))
  names(ret) <- levels(df$daytime)
  
  for (d in levels(df$daytime)) {
    nd_am <- subset(df_train, day > as.POSIXlt("2009-01-01 00:00:00", tz = "UTC") &
                      day < as.POSIXlt("2009-12-31 24:00:00", tz = "UTC")) 
    nd_am <- subset(nd_am, daytime == d)
    nd_am <- subset(nd_am, !duplicated(Datum))
    nd_am$weekday <- factor("Monday",
                            levels = c("Monday", "Tuesday", "Wednesday",
                                       "Thursday", "Friday", "Saturday","Sunday"))
    
    ci_d.FUN <- function(mod, nd, FUN = exp, vcov) {
      
      K <- model.matrix(fm, data = nd[yrw <- seq(from = 1, to = 365, by = 15),])
      K[, c("year2009", colnames(K)[grep("(Intercept)", colnames(K))])] <- 0
      if("cotram" %in% class(mod)) K <- K[,-1] # remove intercept for cotram
      ci <- confint(glht(mod, linfct = K, vcov. = vcov))
      alpha <- attr(ci$confint, "calpha")
      
      K <- model.matrix(fm, data = nd)
      if("cotram" %in% class(mod)) K <- K[,-1] # remove intercept for cotram
      ci <- confint(glht(mod, linfct = K, vcov. = vcov), calpha = alpha)
      
      eJan01 <- ci$confint[1, "Estimate"]
      ci$confint[,1] <- ci$confint[,1] - eJan01
      ci$confint[,2] <- ci$confint[,2] - eJan01
      ci$confint[,3] <- ci$confint[,3] - eJan01
      data.frame(FUN(ci$confint), day = nd$day)
    }
    
    ret[[d]] <- ci_d.FUN(mod = mod, nd = nd_am, vcov = vcov(mod))
    ret[[d]]$model <- NULL
  }
  for (d in names(ret)) ret[[d]]$daytime <- d
  return(ret)
}

fit_mp <- ci.FUN(mp)
fit_mnb <- ci.FUN(mnb)
fit_mcll <- ci.FUN(mcll)

# conditional distribution functions
cdf.FUN <- function(d, q = 0:max(df_train$wild)){
  
  nd <- subset(df, day == as.POSIXlt(d, tz = "UTC"))
  nd$year <- factor(nd$year, levels = levels(df_train$year))
  nd$wild <- NULL
  
  ret <- data.frame(q = q)
  
  ret$p <- ppois(q, lambda = predict(mp, newdata = nd, type = "response"))
  ret$nb <- pnbinom(q, mu = predict(mnb, newdata = nd, type = "response"),
                    size = mnb$theta)
  ret$cll <- c(predict(mcll, newdata = nd, q = q, type = "distribution"))
  
  ret$start <- nd$day
  ret$end <- df[which(df$day == as.POSIXlt(d, tz = "UTC")) + 1, "day"]
  ret$wild <- df[which(df$day == as.POSIXlt(d, tz = "UTC")), "wild"]
  
  return(ret)
}

days <- c("2009-08-29 20:00:00", "2009-04-21 04:30:00",
          "2009-01-01 01:00:00", "2009-02-25 18:00:00")

cdf <- lapply(days, cdf.FUN)

warnings()
