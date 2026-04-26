#' A logistic regression model for testing differential abundance in compositional microbiome data (LOCOM2)
#'
#' This function allows you to test
#' (1). whether any OTU (or taxon) is associated with the trait of interest with FDR control, based on log ratios of relative abundances between pairs of taxa, and
#' (2). whether the whole community is associated with the trait (a global test), based on the harmonic mean method for combining individual p-values
#' The tests accommodate continuous, discrete (binary or categorical), and multivariate traits, and allow adjustment for confounders.
#'
#' This function extends LOCOM (Hu et al., 2022, PNAS) in the following ways:
#' -. accommodating both relative abundance and read count data for OTUs;
#' -. refining the weighting scheme in LOCOM to eliminate confounding by library size;
#' -. incorporating a series of adjustments to ensure stable and reliable inference, even under extreme conditions such as rare taxa and highly unbalanced case–control designs;
#' -. replacing the computationally intensive permutation procedure with a Wald-type test (using a fixed 1,000 permutation replicates).
#'
#'
#' @param otu.table The OTU table (or taxa count table), where rows correspond to samples and columns correspond to OTUs (taxa).
#' @param Y The trait of interest, which can be a vector, matrix, or data frame, must be numeric; for example, a factor should be represented by its corresponding design matrix. When specified as a matrix or data frame, all components are tested jointly for microbial association.
#' @param C The additional (confounding) covariates to be adjusted for. See the requirements for \code{Y}.
#' @param fdr.nominal The nominal FDR level, with a default of 0.1.
#' @param filter A logical value indicating whether to filter out rare taxa. The default is TRUE, using a filtering threshold of min(0.1*n.sam, 10).
#' @param permute A logical value indicating whether to perform permutation. The default is TRUE. 
#' @param n.perm.max The maximum number of permutations. The default is 1,000, used for the Wald-type test. The full permutation procedure as in LOCOM is performed when \code{n.perm.max} is set to NULL. In this case, the total number of permutations is determined as \code{n.otu} * \code{n.rej.stop} * (1/\code{fdr.nominal}),
#' where \code{n.otu} is the number of OTUs (that have non-zero counts in at least one sample). The full permutation procedure adopts a sequential stopping criterion (similar to Sandve et al. 2011),
#' which stops when all taxon-level tests have either reached the prespecified
#' number of rejections (default 100) or yielded a q-value (by the Benjamini-Hochberg [BH] procedure) below the
#' nominal FDR level (default 0.1). 
#' @param n.rej.stop The minimum number of rejections (i.e., instances where the permutation test statistic exceeds the observed test statistic) required before stopping the permutation procedure. The default is 100.
#' @param n.cores The number of cores to be used for parallel computing. The default is 1.
#' @param seed A user-supplied integer seed for the random number generator in the
#'   permutation procedure. The default is NULL, in which case an integer seed is
#'   generated internally at random. In either case, the seed is stored
#'   in the output object to enable reproducibility of the permutation replicates.
#' @param verbose A logical value indicating whether to produce verbose output during the permutation process. The default is TRUE.
#' @param Firth.thresh The threshold (between 0 and 1) of taxon prevalence for applying the Firth correction. The default is 0.4.
#' @return A list consisting of
#' \itemize{
#'   \item p.otu.Wald - Wald p-values for OTU-specific tests
#'   \item q.otu.Wald - Wald q-values (adjusted p-values by BH) for OTU-specific tests
#'   \item detected.otu.Wald - OTUs detected by the Wald test at the nominal FDR level
#'   \item p.otu.perm - permutation p-values for OTU-specific tests
#'   \item q.otu.perm - permutation q-values (adjusted p-values by BH) for OTU-specific tests
#'   \item detected.otu.perm - OTUs detected by the permutation test at the nominal FDR level
#'   \item p.otu.asymptotic - asymptotic p-values for OTU-specific tests
#'   \item q.otu.asymptotic - asymptotic q-values (adjusted p-values by BH) for OTU-specific tests
#'   \item detected.otu.asymptotic - OTUs detected by the asymptotic test at the nominal FDR level
#'   \item beta - effect size at each OTU, defined as beta_j - median (beta_j'), after Yeo–Johnson transformation if the Wald test is used
#'   \item beta.var - estimated variance for each beta
#'   \item ref.otu - reference OTU
#'   \item p.global - p-value for the global test (not available in the asymptotic version). The global test is based on the harmonic mean of individual p-values, using all permutation replicates generated up to the point when the procedure terminates.
#'   \item n.perm.completed - number of permutations completed
#'   \item seed - the seed used to generate the permutation replicates
#' }
#' @import Rcpp
#' @importFrom Rcpp evalCpp
#' @importFrom stats lm p.adjust pchisq resid setNames var
#' @importFrom BiocParallel bplapply SnowParam
#' @importFrom parallel mclapply
#' @importFrom permute shuffleSet
#' @importFrom matrixStats rowMedians colMedians colMins colSds rowRanks
#' @importFrom abind abind
#' @importFrom car powerTransform yjPower
#' @useDynLib LOCOM2
#' @export
#' @examples
#' data("throat.otu.table.filter")
#' data("throat.meta.filter")
#' data("throat.otu.taxonomy")
#' 
#' Y <- ifelse(throat.meta.filter$SmokingStatus == "NonSmoker", 0, 1)
#' C <- ifelse(throat.meta.filter$Sex == "Male", 0, 1)
#' 
#' ##################
#' # running LOCOM2
#' ##################
#' 
#' ## LOCOM2 (Wald), most recommended, better to use n.cores = 4 to speed up
#' 
#' res.Wald <- locom2(otu.table = throat.otu.table.filter, Y = Y, C = C, seed = 123)
#' res.Wald$detected.otu.Wald
#' res.Wald$p.otu.Wald[res.Wald$detected.otu.Wald] 


locom2 <- function(otu.table, Y, C = NULL, 
                   fdr.nominal = 0.1, filter = TRUE, 
                   permute = TRUE, n.perm.max = 1000, n.rej.stop = 100, n.cores = 1, seed = 123,
                   verbose = TRUE, Firth.thresh = 0.4) {
    
    tol <- 1e-10
    iter_max <- 1000

    #-----------------
    # Check input
    #-----------------
    
    otu.table.name <- deparse(substitute(otu.table))
    if (is.null(otu.table)) stop(sprintf("%s must not be NULL", otu.table.name), call. = FALSE)
    if (is.data.frame(otu.table)) otu.table <- data.matrix(otu.table)
    otu.table <- as.matrix(otu.table)
    if (nrow(otu.table) == 0L || ncol(otu.table) == 0L) stop(sprintf("%s must have positive numbers of rows and columns", otu.table.name), call. = FALSE)
    if (!is.numeric(otu.table)) stop(sprintf("%s must be a numeric matrix or data frame", otu.table.name), call. = FALSE)
    if (anyNA(otu.table)) stop(sprintf("%s contains NA values", otu.table.name), call. = FALSE)
    if (any(otu.table < 0)) stop(sprintf("%s contains negative values", otu.table.name), call. = FALSE)
    n.sam <- nrow(otu.table)

    Y.name <- deparse(substitute(Y))
    if (is.null(Y)) stop(sprintf("%s must not be NULL", Y.name), call. = FALSE)
    if (is.data.frame(Y)) Y <- data.matrix(Y)
    Y <- as.matrix(Y)
    if (!is.numeric(Y)) stop(sprintf("%s must be a numeric vector, matrix or data frame", Y.name), call. = FALSE)
    if (anyNA(Y)) stop(sprintf("%s contains NA values", Y.name), call. = FALSE)
    if (nrow(Y) != n.sam) stop(sprintf("Sample sizes in %s (%d) and %s (%d) do not match", otu.table.name, n.sam, Y.name, nrow(Y)), call. = FALSE)
    n.trait <- ncol(Y)
    
    if (!is.null(C)) {
        C.name <- deparse(substitute(C))
        if (is.data.frame(C)) C <- data.matrix(C)
        C <- as.matrix(C)
        if (!is.numeric(C)) stop(sprintf("%s must be a numeric vector, matrix or data frame", C.name), call. = FALSE)
        if (anyNA(C)) stop(sprintf("%s contains NA values", C.name), call. = FALSE)
        if (nrow(C) != n.sam) stop(sprintf("Sample sizes in %s (%d) and %s (%d) do not match", otu.table.name, n.sam, C.name, nrow(C)), call. = FALSE)
    }
    
    if (!is.logical(filter) || length(filter) != 1L || is.na(filter)) {
        stop("filter must be TRUE or FALSE", call. = FALSE)
    }
    if (!is.logical(permute) || length(permute) != 1L || is.na(permute)) {
        stop("permute must be TRUE or FALSE", call. = FALSE)
    }
    if (!is.numeric(fdr.nominal) || length(fdr.nominal) != 1L ||
        is.na(fdr.nominal) || fdr.nominal <= 0 || fdr.nominal >= 1) {
        stop("fdr.nominal must be a single number in (0, 1)", call. = FALSE)
    }
    if (!is.numeric(Firth.thresh) || length(Firth.thresh) != 1L ||
        is.na(Firth.thresh) || Firth.thresh < 0 || Firth.thresh > 1) {
        stop("Firth.thresh must be a single number in [0, 1]", call. = FALSE)
    }
    if (!is.numeric(n.cores) || length(n.cores) != 1L || is.na(n.cores) || n.cores < 1 || n.cores != as.integer(n.cores)) {
        stop("n.cores must be a positive integer", call. = FALSE)
    }
    if (!is.null(n.perm.max) && (!is.numeric(n.perm.max) || length(n.perm.max) != 1L || is.na(n.perm.max) ||
         n.perm.max < 1 || n.perm.max != as.integer(n.perm.max))) {
        stop("n.perm.max must be NULL or a positive integer", call. = FALSE)
    }
    if (!is.numeric(n.rej.stop) || length(n.rej.stop) != 1L || is.na(n.rej.stop) ||
        n.rej.stop < 1 || n.rej.stop != as.integer(n.rej.stop)) {
        stop("n.rej.stop must be a positive integer", call. = FALSE)
    }
    if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1L || is.na(seed))) {
        stop("seed must be NULL or a single numeric value", call. = FALSE)
    }
    
    # Remove rare OTUs
    
    if (filter) {
        filter.thresh <- min(ceiling(0.1*n.sam), 10)
        w <- which(colSums(otu.table > 0) >= filter.thresh)
        if (length(w) < ncol(otu.table)) {
            if (verbose) message(paste(ncol(otu.table) - length(w), ' OTU(s) present in fewer than ', filter.thresh, ' samples are filtered out', sep = ""))
            otu.table <- otu.table[, w, drop = FALSE]
        }
    }
    
    n.otu <- ncol(otu.table)
    otu.names <- colnames(otu.table)
    
    # find reference OTU
    
    lib.size <- rowSums(otu.table)
    if (any(lib.size == 0)) stop("Some samples have zero total counts in otu.table", call. = FALSE)
    ref.otu <- which.max(colMeans(otu.table/lib.size))
    
    # -----------------------
    # Observed statistic
    # -----------------------
    
    if (is.null(C)) {
        X <- cbind(Y, 1)
        Yr <- Y
    } else {
        X <- cbind(Y, 1, C)
        Yr <- stats::resid(stats::lm(Y ~ C))
        Yr <- as.matrix(Yr)
    }
    
    rel.table <- otu.table / lib.size
    rel.table.sum <- rel.table + rel.table[, ref.otu]
    n.X <- ncol(X)
    XX <- CalculateXX(X)
    
    weight <- matrix(1, nrow = n.sam, ncol = n.otu)
    
    prop.presence <- colMeans(otu.table > 0)
    
    beta_init0 <- array(0, dim = c(n.X, n.otu))
    
    res.obs <- tryCatch(
        {
            Newton_LOCOM2(rel_table = rel.table, rel_table_sum = rel.table.sum,
                      X = X, XX = XX,
                      beta_init = beta_init0, weight = weight,
                      tol = tol, iter_max = iter_max, Firth_thresh = Firth.thresh, robust_var = TRUE,
                      prop_presence = prop.presence, get_var = TRUE, step_rate = 1)
        }, 
        error = function(cond) {
            stop("Model fitting failed for the observed data: ", cond$message, call. = FALSE)
        }
    )
    
    beta <- res.obs[[1]][1:n.trait, , drop = FALSE]
    var.median <- apply(beta, 1, mckean.schrader.var)
    beta <- beta - matrixStats::rowMedians(beta)
    beta <- matrix(beta, nrow = n.trait)
    colnames(beta) <- otu.names
    
    # -----------------------
    # Asymptotic results
    # -----------------------
    
    beta.var <- array(NA, dim=c(n.trait, n.trait, n.otu))
    p.otu.asymptotic <- rep(NA, n.otu)
    for (j in 1:n.otu) {
        if (n.trait == 1) {
            beta.var[,,j] <- res.obs[[2]][1:n.trait, 1:n.trait, j] + var.median
        } else {
            beta.var[,,j] <- res.obs[[2]][1:n.trait, 1:n.trait, j] + diag(var.median)
        }
        p.otu.asymptotic[j] <- t(beta[, j, drop = FALSE]) %*% solve(beta.var[,,j]+diag(1e-8, n.trait)) %*% beta[, j, drop = FALSE]
    }
    p.otu.asymptotic <- stats::pchisq(p.otu.asymptotic, df = n.trait, lower.tail = FALSE)
    q.otu.asymptotic <- stats::p.adjust(p.otu.asymptotic, method = "BH")
    p.otu.asymptotic <- stats::setNames(p.otu.asymptotic, otu.names)
    q.otu.asymptotic <- stats::setNames(q.otu.asymptotic, otu.names)
    detected.otu.asymptotic <- otu.names[which(q.otu.asymptotic < fdr.nominal)]
    
    # -----------------------
    # Permutation
    # -----------------------
    
    p.otu.perm <- NULL
    q.otu.perm <- NULL
    detected.otu.perm <- NULL
    
    p.otu.Wald <- NULL
    q.otu.Wald <- NULL
    detected.otu.Wald <- NULL
    
    p.global <- NULL
    n.perm <- NULL
    
    if (permute) {
        
        shrinkage <- ifelse(is.null(C), 0.5, 0.75)
        beta_init <- rbind(matrix(0, nrow = n.trait, ncol = n.otu), shrinkage * res.obs[[1]][(n.trait + 1):n.X, , drop = FALSE])
        
        if (is.null(seed)) seed <- sample(1:10^6, 1)
        set.seed(seed)
        
        if (is.null(n.perm.max)) n.perm.max <- max(n.otu * n.rej.stop * (1 / fdr.nominal), 1000)
        
        beta.perm <- array(NA, dim = c(n.trait, n.otu, n.perm.max))
        n.rej.otu.left <- 0
        n.rej.otu.right <- 0
        
        tol.eq <- 10^-8
        n.perm.block <- 1000
        n.block <- ceiling(n.perm.max / n.perm.block)
        n.perm.core <- ceiling(n.perm.block/n.cores)
        
        parallel.perm <- function(i) {
            
            start <- i*n.perm.core + 1
            end   <- min((i+1)*n.perm.core, ncol(perm.mat))
            if (start > end) return(NULL)
            
            tryCatch(
                {
                    perm_Newton_LOCOM2(rel_table = rel.table, rel_table_sum = rel.table.sum, Yr = Yr, X = X, XX = XX,
                                   beta_init = beta_init, weight = weight,
                                   perm = perm.mat[, start:end, drop=FALSE],
                                   tol = tol, iter_max = iter_max, Firth_thresh = Firth.thresh, robust_var = TRUE,
                                   prop_presence = prop.presence, get_var = TRUE, step_rate = 1)
                }, 
                error = function(cond) {
                    message("Permutation replicate error: ", cond$message)
                    message("Second attempt with an adjusted step size...")
                    return(tryCatch(
                        {
                            perm_Newton_LOCOM2(rel_table = rel.table, rel_table_sum = rel.table.sum, Yr = Yr, X = X, XX = XX,
                                       beta_init = beta_init, weight = weight,
                                       perm = perm.mat[, start:end, drop=FALSE],
                                       tol = tol, iter_max = iter_max * 10, Firth_thresh = Firth.thresh, robust_var = TRUE,
                                       prop_presence = prop.presence, get_var = TRUE, step_rate = 0.1)
                        }, 
                        error = function(cond2) {
                            message("Second attempt failed: ", cond2$message)
                            NULL
                        }
                    ))
                }
            )
        }
        
        n.perm <- 0
        
        for (i.block in 1:n.block) {
            
            n_this <- min(n.perm.block, n.perm.max - n.perm)
            perm.mat <- t(permute::shuffleSet(n.sam, n_this, quietly = TRUE)) - 1
            
            if (n.cores > 1) {
                
                worker.id <- 0:(min(n.cores, ceiling(ncol(perm.mat) / n.perm.core)) - 1)
                
                if (Sys.info()[['sysname']] == 'Windows') {
                    parallel.stat <- BiocParallel::bplapply(worker.id, parallel.perm, BPPARAM = BiocParallel::SnowParam(workers = length(worker.id))) # replaced MulticoreParam
                } else {
                    parallel.stat <- parallel::mclapply(worker.id, parallel.perm, mc.cores = length(worker.id))
                }
                
                # check whether there is any error in any core
                if (any(vapply(parallel.stat, is.null, logical(1)))) next
                
                res.perm <- do.call(abind::abind, parallel.stat)
                
            } else {
                
                res.perm <- tryCatch(
                    {
                        perm_Newton_LOCOM2(rel_table = rel.table, rel_table_sum = rel.table.sum, Yr = Yr, X = X, XX = XX,
                                       beta_init = beta_init, weight = weight,
                                       perm = perm.mat,
                                       tol = tol, iter_max = iter_max, Firth_thresh = Firth.thresh, robust_var = TRUE,
                                       prop_presence = prop.presence, get_var = FALSE, step_rate = 1)
                    }, 
                    error = function(cond) {
                        message("Permutation replicate error: ", cond$message)
                        message("Second attempt with an adjusted step size...")
                        
                        return(tryCatch(
                            {
                                perm_Newton_LOCOM2(rel_table = rel.table, rel_table_sum = rel.table.sum, Yr = Yr, X = X, XX = XX,
                                           beta_init = beta_init, weight = weight,
                                           perm = perm.mat,
                                           tol = tol, iter_max = iter_max*10, Firth_thresh = Firth.thresh, robust_var = TRUE,
                                           prop_presence = prop.presence, get_var = FALSE, step_rate = 0.1)
                            }, 
                            error = function(cond2) {
                                message("Second attempt failed: ", cond2$message)
                                NULL
                            }
                        ))
                    }
                )
                
            } # End of parallel
            
            
            # Permutation batch skipped due to fitting failure
            if (is.null(res.perm) || !is.array(res.perm) || anyNA(res.perm)) next
            
            n.valid <- dim(res.perm)[3]
            w <- (n.perm + 1):(n.perm + n.valid)
            n.perm <- n.perm + n.valid
            n.perm.inv <- 1 / (n.perm + 1)
            
            for (i.trait in 1:n.trait) {
                beta.perm[i.trait, , w] <- t(t(res.perm[i.trait, , ]) - matrixStats::colMedians(res.perm[i.trait, , ]))
            }
            diff <- sweep(beta.perm[, , w, drop = FALSE], c(1, 2), beta, FUN = "-")
            n.rej.otu.equal <- 0.5 * (abs(diff) < tol.eq)
            n.rej.otu.left <- n.rej.otu.left + rowSums((diff < -tol.eq) + n.rej.otu.equal, dims = 2)
            n.rej.otu.right <- n.rej.otu.right + rowSums((diff > tol.eq) + n.rej.otu.equal, dims = 2)
            n.rej.otu <- 2 * pmin(n.rej.otu.left, n.rej.otu.right) + 1
            
            if (n.trait >= 2) {
                pmin.otu <- matrixStats::colMins(n.rej.otu)
                pnull.otu <- n.perm + 0.5 - array(matrixStats::rowRanks(abs(beta.perm[, , 1:n.perm, drop = FALSE]), ties.method = "average", dim. = c(n.trait * n.otu, n.perm)), dim = c(n.trait, n.otu, n.perm))
                pnullmin.otu <- array(matrixStats::colMins(pnull.otu, dim. = c(n.trait, n.otu * n.perm)), dim = c(n.otu, n.perm))
                diff.1 <- pnullmin.otu - c(pmin.otu)
                n.rej.otu <- rowSums((diff.1 < -tol.eq) + 0.5 * (abs(diff.1) < tol.eq))
            } 
            p.otu.perm <- n.rej.otu * n.perm.inv
            q.otu.perm <- fdr.Sandve(p.otu.perm)
            
            if (verbose) message(paste("permutations:", n.perm))
                                 
            if (all(q.otu.perm <= fdr.nominal | n.rej.otu >= n.rej.stop)) break

        } # End of permutation loop
        
        q.otu.perm <- stats::p.adjust(p.otu.perm, method = "BH")
        p.otu.perm <- stats::setNames(p.otu.perm, otu.names)
        q.otu.perm <- stats::setNames(q.otu.perm, otu.names)
        detected.otu.perm <- otu.names[which(q.otu.perm < fdr.nominal)] # previously <=
        
        # ----------------------------------------
        # Global p-value
        # ----------------------------------------
        
        if (n.trait == 1) {
            beta.all <- cbind(beta[1, ], beta.perm[1, , 1:n.perm])
        } else {
            beta.all <- cbind(pmin.otu, pnullmin.otu)
        }
        r <- matrixStats::rowRanks(beta.all, ties.method = "average")
        p.otu1 <- (2 * pmin(r, n.perm + 2 - r) - 1.5)
        
        stat.global <- colSums(1 / p.otu1)
        p.global <- (sum(stat.global[1] <= stat.global[-1]) + 1) * n.perm.inv
        
        # ----------------------------------------
        # Wald test (allowing multiple traits)
        # ----------------------------------------
        
        # Yeo-Johnson transformation
        
        if (n.perm < 100) stop("Not enough permutations to run Wald transformation-based test", call. = FALSE)
        n.Wald <- min(1000, n.perm)
        z2.value.yj <- rep(NA, n.otu)
        beta.Wald <- beta
        beta.var.Wald <- beta.var

        for (j in 1:n.otu) {
            beta.perm.j <- matrix(beta.perm[,j,1:n.Wald], nrow=n.trait, ncol=n.Wald)
            beta.var.Wald[,,j] <- stats::var(t(beta.perm.j), na.rm = TRUE)
            
            z2.value.yj[j] <- tryCatch(
                {
                    for (t in 1:n.trait) {
                        pt <- car::powerTransform(beta.perm.j[t,] ~ 1, family = "yjPower")

                        beta.perm.j[t,] <- car::yjPower(beta.perm.j[t,], pt$lambda)
                        beta.Wald[t,j] <- car::yjPower(beta[t,j], pt$lambda) - mean(beta.perm.j[t,], na.rm = TRUE)
                    }
                    beta.var.Wald[,,j] <- stats::var(t(beta.perm.j), na.rm = TRUE)
                    t(beta.Wald[,j,drop=FALSE]) %*% solve(beta.var.Wald[,,j]+diag(1e-8, n.trait)) %*% beta.Wald[,j,drop=FALSE]
                },
                 warning = function(w) {
                    return(t(beta.Wald[,j,drop=FALSE]) %*% solve(beta.var.Wald[,,j]+diag(1e-8, n.trait)) %*% beta.Wald[,j,drop=FALSE])
                }
            )
        }

        p.otu.Wald <- stats::pchisq(z2.value.yj, df = n.trait, lower.tail = FALSE)
        q.otu.Wald <- stats::p.adjust(p.otu.Wald, method = "BH")
        p.otu.Wald <- stats::setNames(p.otu.Wald, otu.names)
        q.otu.Wald <- stats::setNames(q.otu.Wald, otu.names)
        detected.otu.Wald <- otu.names[which(q.otu.Wald < fdr.nominal)]
        
    } # if (permute)
    
    if (n.trait==1) {
        beta <- beta[1,]
        beta.var <- beta.var[1,1,]
        if (permute) {
            beta.Wald <- beta.Wald[1,]
            beta.var.Wald <- beta.var.Wald[1,1,]
        }
    }
    
    if (permute) {
        beta.final <- beta.Wald
        beta.var.final <- beta.var.Wald
    } else {
        beta.final <- beta
        beta.var.final <- beta.var
    }
        
    results <- list(
            p.otu.Wald = p.otu.Wald,
            q.otu.Wald = q.otu.Wald,
            detected.otu.Wald = detected.otu.Wald,
            
            p.otu.perm = p.otu.perm,
            q.otu.perm = q.otu.perm,
            detected.otu.perm = detected.otu.perm,
            
            p.otu.asymptotic = p.otu.asymptotic,
            q.otu.asymptotic = q.otu.asymptotic,
            detected.otu.asymptotic = detected.otu.asymptotic,
            
            beta = beta.final,
            beta.var = beta.var.final, 
            ref.otu = ref.otu,
            
            p.global = p.global,
            n.perm.completed = n.perm,
            seed = seed
    )
    return(results)
    
} # LOCOM2


fdr.Sandve = function(p.otu) {
    m = length(p.otu)
    p.otu.sort = sort(p.otu)
    n.otu.detected = seq(1, m)
    pi0 = min(1, 2/m*sum(p.otu))
    
    qval.sort = m * pi0 * p.otu.sort / n.otu.detected
    j.min.q = 1
    while (j.min.q < m) {
        min.q = min( qval.sort[j.min.q:m] )
        new.j.min.q = (j.min.q-1) + max( which(qval.sort[j.min.q:m]==min.q) )
        qval.sort[j.min.q:new.j.min.q] = qval.sort[new.j.min.q]
        j.min.q = new.j.min.q+1
    }
    mat = match(p.otu, p.otu.sort)
    qval.orig = qval.sort[mat]
    results = qval.orig
    
    return(results)
    
} # fdr.Sandve


mckean.schrader.var = function( x ) {
    
    # The mckean-schrader standard error and variance of a sample median
    
    n.data = length(x)
    c = max(1, round( (n.data+1)/2 - 0.98*sqrt(n.data) ) )    #  0.98=1.96/sqrt(4)
    x.sort = sort(x)
    med.se = ( x.sort[n.data-c+1] - x.sort[c] )/3.92     #    3.92=2*1.96
    med.var = med.se^2
    
    return(med.var)
    
} # mckean.schrader.var
