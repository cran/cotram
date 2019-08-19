
predict.cotram <- function(object, newdata = model.frame(object), smooth = FALSE,
                           type = c("lp", "trafo", "distribution", "survivor", "density", 
                                    "logdensity", "cumhazard", "quantile"),
                           q = NULL, K = 50, ...) {
    type <- match.arg(type)
    
    y <- variable.names(object, "response")
    
    if (!missing(newdata)) {
        if (object$plus_one && variable.names(object, "response") %in% names(newdata))
            newdata[,y] <- newdata[,y] + 1L
        if (!is.data.frame(newdata)) { ### newdata is list with _all_ variables
            stopifnot(is.null(q))
            if (type != "quantile") ### hm, why???
                stopifnot(variable.names(object, "response") %in% names(newdata))
        }
        if (!(y %in% names(newdata)) && is.null(q)) {
            q <- mkgrid(object, n = K)[[y]] 
            if (smooth)
                q <- seq(from = min(q), to = max(q), length.out = K)
        }
    } else {
        newdata <- NULL
        if (is.null(q)) {
            q <- mkgrid(object, n = K)[[y]] 
            if (smooth)
                q <- seq(from = min(q), to = max(q), length.out = K)
        }
    }
    q <- q + as.integer(object$plus_one)
    
    # predict.tram seems not to be exported
    if (type == "lp") {
        ret <- model.matrix(object, data = newdata) %*% 
            coef(object, with_baseline = FALSE)
        if (object$negative) return(-ret)
        return(ret)
    }
    
    if (smooth) {
        ret <- predict(as.mlt(object), newdata = newdata, type = type, q = q, K = K, ...)
        if (type == "quantile") {
            if (is.matrix(ret)) 
                zero <- matrix(0, nrow = nrow(ret), ncol = ncol(ret))
            ret <- pmax(zero, ret - as.integer(object$plus_one))
        }
    } else {
        if (type == "quantile") 
            stop("quantiles only available with smooth = TRUE")
        if (length(grep("density", type)) > 0) {
            ret <- predict(as.mlt(object), newdata = newdata, 
                           type = "distribution", q = q, K = K, ...) - 
                predict(as.mlt(object), newdata = newdata, 
                        type = "distribution", q = q - 1, K = K, ...)
            if (type == "logdensity") ret <- log(ret)
        } else {
            ret <- predict(as.mlt(object), newdata = newdata, type = type, q = q, K = K, ...)
        }
    }
    return(ret)
}
