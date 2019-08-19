\name{cotram-methods}
\alias{predict.cotram}
\alias{plot.cotram}
\alias{as.mlt.cotram}
\alias{logLik.cotram}
\title{
  Methods for Count Transformation Models
}
\description{
  Methods for objects inheriting from class cotram
}
\usage{
\method{predict}{cotram}(object, newdata = model.frame(object), smooth = FALSE,
        type = c("lp", "trafo", "distribution", "survivor",
                 "density", "logdensity", "cumhazard", "quantile"), 
        q = NULL, K = 50, ...)
\method{plot}{cotram}(x, newdata, smooth = FALSE, 
     type = c("distribution", "survivor", "density", "logdensity", 
              "cumhazard", "quantile", "trafo"),
     q = NULL, prob = 1:(K - 1) / K, K = 50, col = rgb(.1, .1, .1, .1), 
     lty = 1, lwd = 1, add = FALSE, ...) 
\method{as.mlt}{cotram}(object)
\method{logLik}{cotram}(object, parm = coef(as.mlt(object), fixed = FALSE), newdata, ...)
}
\arguments{
  \item{object, x}{a fitted linear count transformation model inheriting
                from class \code{cotram}.}
  \item{newdata}{an optional data frame of observations.}
  \item{parm}{model parameters.}
  \item{type}{type of prediction, current options include
              linear predictors (\code{"lp"}, of \code{x} variables in the
              formula \code{y ~ x}), transformation functions
              (\code{"trafo"}) or distribution functions on the
              scale of the cdf (\code{"distribution"}), survivor function,
              density function, log-density function, cumulative
              hazard function or quantile function.}
  \item{q}{quantiles at which to evaluate the model.}
  \item{smooth}{logical; if \code{TRUE} a smoothed function of \code{type} is returned.}
  \item{prob}{probabilities for the evaluation of the quantile function (\code{type = "quantile"}).}
  \item{K}{number of grid points to generate (in the absence of \code{q}).}
  \item{col}{color for the lines to plot.}
  \item{lty}{line type for the lines to plot.}
  \item{lwd}{line width.}
  \item{add}{logical; indicating if a new plot shall be generated (the default)}
  \item{\dots}{additional arguments to the underlying methods for class
               \code{tram} or \code{mlt}, see \code{\link[tram]{tram-methods}}
               	or \code{\link[mlt]{mlt-methods}}.}
}
\details{
	\code{predict} and \code{plot} can be used to inspect the model on
  different scales.
}
\seealso{
  \code{\link[tram]{tram-methods}}, \code{\link[mlt]{mlt-methods}},  \code{\link[mlt]{plot.ctm}}
}
\references{

  Torsten Hothorn, Lisa Moest, Peter Buehlmann (2018), Most Likely
  Transformations, \emph{Scandinavian Journal of Statistics}, \bold{45}(1),
  110--134, \doi{10.1111/sjos.12291}.

  Sandra Siegfried, Torsten Hothorn (2019), Count Transformation Models,
  submitted manuscript (included in \code{cotram} package vignette).

}
\examples{
	
  data("birds", package = "TH.data")
  
  ### fit count transformation model with cloglog link
  m_birds <- cotram(SG5 ~ AOT + AFS + GST + DBH + DWC + LOG, data = birds,
                    method = "cloglog")
  logLik(m_birds)

  ### classical likelihood inference
  summary(m_birds)

  ### coefficients of the linear predictor (discrete hazard ratios)
  exp(-coef(m_birds))

  ### compute predicted median along with 10% and 90% quantile for the first
  ### three observations
  nd <- birds[1:3, -10] ## remove response
  predict(m_birds, newdata = nd, type = "quantile", prob = c(.1, .5, .9), 
          smooth = TRUE)

  ### plot the predicted distribution for these observations
  plot(m_birds, newdata = nd, type = "distribution",
       col = c("skyblue", "grey", "seagreen"))
    
}
