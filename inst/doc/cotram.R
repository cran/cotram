## ----setup, echo = FALSE, results = "hide", message = FALSE, warning = FALSE----
set.seed(290875)

sapply(c("cotram", "tram", "mlt", "lattice"), library, char = TRUE)

trellis.par.set(list(plot.symbol = list(col=1,pch=20, cex=0.7),
                     box.rectangle = list(col=1),
                     box.umbrella = list(lty=1, col=1),
                     strip.background = list(col = "white"), 
                     axis.components = list(top = list(tck = 0),
                                            right = list(tck = 0)))
                )
ltheme <- canonical.theme(color = FALSE)     ## in-built B&W theme
ltheme$strip.background$col <- "transparent" ## change strip bg
lattice.options(default.theme = ltheme)

panel <- function(x, y, ...) {
  panel.abline(h = 1, col = "grey")
  panel.xyplot(x, y, type = "l", col = "black")
}

knitr::opts_chunk$set(echo = TRUE, results = 'markup', error = FALSE,
                      warning = FALSE, message = FALSE,
                      tidy = FALSE, cache = FALSE, size = "small",
                      fig.width = 6, fig.height = 4, fig.align = "center",
                      out.width = NULL, ###'.6\\linewidth', 
                      out.height = NULL,
                      fig.scap = NA)
knitr::render_sweave()  # use Sweave environments
knitr::set_header(highlight = '')  # do not \usepackage{Sweave}

## R settings
opt <- options(prompt = "R> ", continue = "+  ", useFancyQuotes = FALSE, width = 75)  # JSS style
library("colorspace")
col <- diverge_hcl(10, h = c(246, 40), c = 96, l = c(65, 90))
fill <- diverge_hcl(10, h = c(246, 40), c = 96, l = c(65, 90), alpha = .3)

## ----citation, echo = FALSE----------------------------------------------
year <- substr(packageDescription("cotram")$Date, 1, 4)
version <- packageDescription("cotram")$Version

## ----DVC-data, echo = FALSE, results = "hide", message = FALSE-----------
obs <- NULL
if (!file.exists("analysis/DVC.rda")) {
    op <- options(timeout = 120)
    dwnfl <- try(download.file("https://zenodo.org/record/17179/files/DVC.tgz", "DVC.tgz"))
    options(op)
    if (!inherits(dwnfl, "try-error")) {
        untar("DVC.tgz", file = "analysis/DVC.rda")
        load("analysis/DVC.rda")
    }
} else {
    load("analysis/DVC.rda")
}

## ----fail, results = "asis", echo = FALSE--------------------------------
if (is.null(obs)) {
    cat("Downloading data from zenodo failed, stop processing.", "\\end{document}\n")
    knitr::knit_exit()
}

## ----DVC-setup, echo = FALSE, results = "hide", message = FALSE----------
# set-ups
loc <- Sys.setlocale("LC_ALL", "en_US.UTF-8")
rm(loc)

# data-frame setup
df <- data.frame(margin.table(obs[,,"wild"], 2))
colnames(df) <- c("day", "DVC")
df$day <- as.Date(df$day)
df$weekday <- factor(format(df$day, "%A"),
                      levels = c("Monday", "Tuesday", "Wednesday", "Thursday",
                                 "Friday", "Saturday", "Sunday"))
df$weekday[weekdays$ArbZeitFaktor == 0] <- "Sunday"
df$year <- factor(format(df$day, "%Y"))
df$time <- as.numeric(difftime(df$day, start, unit = "days"))

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

## ----mod_HR--------------------------------------------------------------
mod_cloglog <- cotram(DVC ~ year + weekday + tvar1 + tvar2 + tvar3 +
                      tvar4 + tvar5 + tvar6 + tvar7 + tvar8 + tvar9 + 
                      tvar10 + tvar11 + tvar12 + tvar13 + tvar14 + 
                      tvar15 + tvar16 + tvar17 + tvar18 + tvar19 + tvar20,
                      data = df, method = "cloglog")
logLik(mod_cloglog)

## ----nd_day--------------------------------------------------------------
nd <- model.frame(mod_cloglog)[which(df$year == "2011"), -1]
nd$day <- df[which(df$year == "2011"), "day"]
nd$weekday <- factor("Monday", levels = levels(nd$weekday))

## ----cloglog_lp----------------------------------------------------------
fit_cloglog <- predict(mod_cloglog, type = "lp", newdata = nd) - 
  predict(mod_cloglog, type = "lp", newdata = nd)[1]
xyplot(exp(fit_cloglog) ~ day , data = cbind(nd, fit_cloglog),
       ylab = "Hazard ratio", xlab = "Day of year", panel = panel)

## ----mod_logit-----------------------------------------------------------
mod_logit <- cotram(DVC ~ year + weekday + tvar1 + tvar2 + tvar3 +
                    tvar4 + tvar5 + tvar6 + tvar7 + tvar8 + tvar9 +
                    tvar10 + tvar11 + tvar12 + tvar13 + tvar14 +
                    tvar15 + tvar16 + tvar17 + tvar18 + tvar19 + tvar20,
                    data = df, method = "logit")
logLik(mod_logit)

## ----coefs_logit---------------------------------------------------------
years <- grep("year", names(coef(mod_logit)), value = TRUE)
coef <- exp(-coef(mod_logit)[years])
ci <- exp(-confint(mod_logit)[years,])
round(cbind(coef, ci), 3)

## ----logit_cdf, echo = FALSE---------------------------------------------
nd <- model.frame(mod_logit)[unique(df$year), -1]
nd$year <- factor(levels(nd$year), levels = levels(nd$year))
nd$weekday <- factor("Monday", levels = levels(nd$weekday))
nd[, grep("tvar", colnames(nd))] <- 0
plot(mod_logit, type = "distribution", newdata = nd, q = 20:200, smooth = TRUE,
     xlab = "Number deer-vehicle collisions", ylab = "Distribution function",
     col = col, lwd = 1.5)
legend(120, .6, legend = levels(nd$year)[1:5], lty = 1, lwd = 1.5, col = col[1:5],
       bty = "n")
legend(160, .6, legend = levels(nd$year)[6:10], lty = 1, lwd = 1.5, col = col[6:10],
       bty = "n")

## ----mod_loglog----------------------------------------------------------
mod_loglog <- cotram(DVC ~ year + weekday + tvar1 + tvar2 + tvar3 +
                     tvar4 + tvar5 + tvar6 + tvar7 + tvar8 + tvar9 +
                     tvar10 + tvar11 + tvar12 + tvar13 + tvar14 +
                     tvar15 + tvar16 + tvar17 + tvar18 + tvar19 + tvar20,
                     data = df, method = "loglog")
logLik(mod_loglog)

## ----nd_density----------------------------------------------------------
nd <- model.frame(mod_loglog)[1,]

## ----loglog_density------------------------------------------------------
plot(mod_loglog, type = "density", newdata = nd, q = 0:150, col = col,
     xlab = "Number of deer-vehicle collisions", ylab = "Density function")
abline(v = nd$DVC)

## ----probit_mod----------------------------------------------------------
mod_probit <- cotram(DVC ~ year + weekday + tvar1 + tvar2 + tvar3 +
                     tvar4 + tvar5 + tvar6 + tvar7 + tvar8 + tvar9 +
                     tvar10 + tvar11 + tvar12 + tvar13 + tvar14 +
                     tvar15 + tvar16 + tvar17 + tvar18 + tvar19 + tvar20,
                     data = df, method = "probit")
logLik(mod_probit)

## ----probit_trafo--------------------------------------------------------
nd <- model.frame(mod_probit)[1, ]
trafo_probit <- predict(mod_probit, type = "trafo",
                        newdata = nd, smooth = TRUE)

## ----probit_cb-----------------------------------------------------------
cb_probit <- confband(mod_probit, type = "trafo",
                      newdata = nd, smooth = TRUE)

## ----probit_plot, echo = FALSE-------------------------------------------
layout(matrix(1:2, nrow = 1))
plot(mod_probit, type = "trafo", newdata = df[1,], smooth = TRUE, 
     xlab = "Number of deer-vehicle collisions",
     ylab = expression(paste("Transformation function ", alpha(y))),
     col = col[10], lwd = 2, confidence = "band")
plot(mod_probit, type = "distribution", newdata = df[1,], smooth = TRUE, 
     xlab = "Number of deer-vehicle collisions",
     ylab = "Distribution function",
     col = col[1], lwd = 2, confidence = "band")

## ----sessionInfo, echo = FALSE, results = "hide"------------------------------
sessionInfo()
options(opt)

