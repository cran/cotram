
predict.cotram <- function(object, newdata = model.frame(object), 
                           type = c("lp", "trafo", "distribution", "survivor", "density", 
                                    "logdensity", "hazard", "loghazard", 
                                    "cumhazard", "logcumhazard", "odds", "logodds", 
                                    "quantile"),
                           smooth = FALSE, q = NULL, K = 20, prob = 1:(10-1)/10, ...) {
    type <- match.arg(type)
    
    y <- variable.names(object, "response")
    
    ## original newdata & q for internal calls of predict.cotram
    nd <- newdata
    nq <- q
    
    ## y + 1 for log_first
    plus_one <- as.integer(object$log_first)
    
    if (y %in% names(newdata)) {
        if (any(newdata[,y] < 0)) stop("response is non-positive")
        if (!smooth && !all(newdata[,y] %% 1 == 0)) stop("response is non-integer")
        newdata[,y] <- newdata[,y] + plus_one
    }
    if (!is.null(q)) {
        if (any(q < 0)) stop("q is non-positive")
        if (!smooth && !all(q %% 1 == 0)) stop("q is non-integer")
        q <- q + plus_one 
    }
    ## generate quantiles if response not in newdata & no q given
    if (!(y %in% names(newdata)) && is.null(q)) {
        q <- mkgrid(object, n = K)[[y]] + plus_one
        if (smooth)
            q <- seq(from = min(q), to = max(q), length.out = K)
    }
    
    
    ## linear predictor
    if (type == "lp") {
        ret <- model.matrix(object, data = newdata) %*% 
            coef(object, with_baseline = FALSE)
        if (object$negative) return(-ret)
        return(ret)
    }
    
    
    if (smooth) {
        
        ret <- predict(as.mlt(object), newdata = newdata, type = type, q = q,
                       K = K, prob = prob, ...)
        
        ### smooth quantile function
        if (type == "quantile") {
            ret <- ret - plus_one
            names <- dimnames(ret)
            zero <- array(0, dim = length(ret))
            if (is.matrix(ret)) zero <- matrix(0, nrow = nrow(ret), ncol = ncol(ret))
            ret <- pmax(zero, ret)
            dimnames(ret) <- names
            return(ret)
        }
        
    } else {
        
        if (type == "quantile")
            stop("quantiles only available with smooth = TRUE")
        
        ret <- predict(as.mlt(object), newdata = newdata, type = type, q = q,
                       K = K, prob = prob, ...)
        
        ### discrete density function
        if (length(grep("density", type)) > 0) {
            
            newdata_m1 <- newdata
            if (y %in% names(newdata_m1)) newdata_m1[[y]] <- newdata_m1[[y]] - 1L
            
            q_m1 <- q
            if (!is.null(q_m1)) q_m1 <- q_m1 - 1L
            
            ret <- predict(as.mlt(object), newdata = newdata, 
                           type = "distribution", q = q, K = K, prob = prob, ...) - 
                predict(as.mlt(object), newdata = newdata_m1, 
                        type = "distribution", q = q_m1, K = K, prob = prob, ...)
            if (type == "logdensity") ret <- log(ret)
        }

        ### discrete hazard function
        if (type %in% c("hazard", "loghazard")) {
          ## use predict.cotram
            d <- predict(object, newdata = nd, q = nq, type = "density", ...)
            p <- predict(object, newdata = nd, q = nq, type = "distribution", ...)
            if (type == "loghazard") return(d - log1p(-(p - exp(d))))
            return(d / (1 - (p - d)))
        }
    }
    
    ## correct dimnames
    if (!is.null(dimnames(ret)[[y]]))
      dimnames(ret)[[y]] <- as.numeric(dimnames(ret)[[y]]) - plus_one
    
    return(ret)
}
