function (model, cluster, parallel = FALSE, use_white = NULL,
    df_correction = TRUE, leverage = FALSE, force_posdef = FALSE,
    stata_fe_model_rank = FALSE, debug = FALSE)
{
    if (inherits(cluster, "formula")) {
        cluster_tmp <- expand.model.frame(model, cluster, na.expand = FALSE)
        cluster <- model.frame(cluster, cluster_tmp, na.action = na.pass)
    }
    else {
        cluster <- as.data.frame(cluster, stringsAsFactors = FALSE)
    }
    cluster_dims <- ncol(cluster)
    tcc <- 2^cluster_dims - 1
    acc <- list()
    for (i in 1:cluster_dims) {
        acc <- append(acc, combn(1:cluster_dims, i, simplify = FALSE))
    }
    if (debug)
        print(acc)
    vcov_sign <- sapply(acc, function(i) (-1)^(length(i) + 1))
    acc <- acc[-1:-cluster_dims]
    if (debug)
        print(acc)
    if (!is.null(model$na.action)) {
        if (class(model$na.action) == "exclude") {
            cluster <- cluster[-model$na.action, ]
            esttmp <- estfun(model)[-model$na.action, ]
        }
        else if (class(model$na.action) == "omit") {
            cluster <- cluster[-model$na.action, ]
            esttmp <- estfun(model)
        }
        cluster <- as.data.frame(cluster, stringsAsFactors = FALSE)
    }
    else {
        esttmp <- estfun(model)
    }
    if (debug)
        print(class(cluster))
    i <- !sapply(cluster, is.numeric)
    cluster[i] <- lapply(cluster[i], as.character)
    if (cluster_dims > 1) {
        for (i in acc) {
            cluster <- cbind(cluster, Reduce(paste0, cluster[,
                i]))
        }
    }
    df <- data.frame(M = integer(tcc), N = integer(tcc), K = integer(tcc))
    rank_adjustment <- 0
    if (stata_fe_model_rank == TRUE) {
        rank_adjustment <- 1
    }
    for (i in 1:tcc) {
        df[i, "M"] <- length(unique(cluster[, i]))
        df[i, "N"] <- length(cluster[, i])
        df[i, "K"] <- model$rank + rank_adjustment
    }
    if (df_correction == TRUE) {
        df$dfc <- (df$M/(df$M - 1)) * ((df$N - 1)/(df$N - df$K))
    }
    else if (is.numeric(df_correction) == TRUE) {
        df$dfc <- df_correction
    }
    else {
        df$dfc <- 1
    }
    if (is.null(use_white)) {
        if (cluster_dims > 1 && df[tcc, "M"] == prod(df[-tcc,
            "M"])) {
            use_white <- TRUE
        }
        else {
            use_white <- FALSE
        }
    }
    if (use_white) {
        df <- df[-tcc, ]
        tcc <- tcc - 1
    }
    if (debug) {
        print(acc)
        print(paste("Original Cluster Dimensions", cluster_dims))
        print(paste("Theoretical Cluster Combinations", tcc))
        print(paste("Use White VCOV for final matrix?", use_white))
    }
    if (leverage) {
        if ("x" %in% names(model)) {
            X <- model$x
        }
        else {
            X <- model.matrix(model)
        }
        ixtx <- solve(crossprod(X))
        h <- 1 - vapply(1:df[1, "N"], function(i) X[i, ] %*%
            ixtx %*% as.matrix(X[i, ]), numeric(1))
        if (leverage == 3) {
            esttmp <- esttmp/h
        }
        else if (leverage == 2) {
            esttmp <- esttmp/sqrt(h)
        }
    }
    uj <- list()
    if (length(parallel) > 1) {
        clusterExport(parallel, varlist = c("cluster", "model"),
            envir = environment())
    }
    for (i in 1:tcc) {
        if (length(parallel) > 1) {
            uj[[i]] <- crossprod(parApply(parallel, esttmp, 2,
                function(x) tapply(x, cluster[, i], sum)))
        }
        else {
            uj[[i]] <- crossprod(apply(esttmp, 2, function(x) tapply(x,
                cluster[, i], sum)))
        }
    }
    if (debug) {
        print(df)
        print(uj)
        print(vcov_sign)
    }
    vcov_matrices <- list()
    for (i in 1:tcc) {
        vcov_matrices[[i]] <- vcov_sign[i] * df[i, "dfc"] * (uj[[i]]/df[i,
            "N"])
    }
    if (use_white) {
        i <- i + 1
        vcov_matrices[[i]] <- vcov_sign[i] * meatHC(model, type = "HC0")
    }
    if (debug) {
        print(vcov_matrices)
    }
    vcov_matrix <- sandwich(model, meat. = Reduce("+", vcov_matrices))
    if (force_posdef) {
        decomp <- eigen(vcov_matrix, symmetric = TRUE)
        if (debug)
            print(decomp$values)
        pos_eigens <- pmax(decomp$values, rep.int(0, length(decomp$values)))
        vcov_matrix <- decomp$vectors %*% diag(pos_eigens) %*%
            t(decomp$vectors)
    }
    return(vcov_matrix)
}
<bytecode: 0x00000000203057c8>
<environment: namespace:multiwayvcov>
