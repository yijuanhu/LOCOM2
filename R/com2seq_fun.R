#' com2seq: combining two sequencing datasets (e.g., from 16S marker-gene and shotgun metagenomic sequencing) for integrative analysis of differential abundance at both the genus and community levels
#'
#' This function allows you to test
#' (1). whether any genus is associated with the trait of interest with FDR control, based on the log ratio of relative abundances between pairs of taxa, and
#' (2). whehter the whole community is associated with the trait (a global test).
#' The tests accommodate continuous, discrete (binary, categorical), and multivariate traits, and allows adjustment of confounders.
#'
#' This function uses a sequential stopping criterion (similar to that of Sandve et al. 2011) for the permutation procedure,
#' which stops when all OTU-level tests have either reached the pre-specified
#' number (default 100) of rejections or yielded a q-value (by the Benjamini-Hochberg [BH] procedure) that is below the
#' nominal FDR level (default 0.2). The permutation procedure is always terminated if a pre-specified maximum number (see description of \code{n.perm.max}) of
#' permutations have been generated. The global test is based on all permutation replicates when the procedure stops/terminates.
#'
#'
#' @param formula.1 A formula in the form of \code{otu.table.1 | C1 ~ Y1}, where \code{otu.table.1} is the OTU table (or taxa count table) for the first sequencing dataset, \code{C1} contains confounding covariates to be adjusted for, and \code{Y1} is the trait of interest. In the OTU table, rows correspond to samples and columns correspond to OTUs; both row and column names are required for aligning two OTU tables. \code{C1} can be replaced by multiple (say \code{K}) covariates in the form of \code{C11 + C12 + ... + C1K}. 
#' @param data.1 A data frame for metadata of the first sequencing dataset, with samples ordered in the same way as in \code{otu.table.1}.
#' @param formula.2 A formula in the form of \code{otu.table.2 | C2 ~ Y2}, where \code{otu.table.2} is the OTU table (or taxa count table) for the second sequencing dataset, \code{C2} contains confounding covariates to be adjusted for, and \code{Y2} is the trait of interest. Same requirements for each component as for \code{formula.1}.
#' @param data.2 A data frame for metadata of the second sequencing dataset, with samples ordered in the same way as in \code{otu.table.2}.
#' @param fdr.nominal the nominal FDR. The default is 0.2.
#' @param seed a user-supplied integer seed for the random number generator in the
#'   permutation procedure. The default is NULL, which means that an integer seed will be
#'   generated internally and randomly. In either case, the integer seed will be stored
#'   in the output object in case the user wants to reproduce the permutation replicates.
#' @param n.perm.max the maximum number of permutations. The default is NULL, in which case \code{n.taxa} * \code{n.rej.stop} * (1/\code{fdr.nominal})
#'   is used where \code{n.taxa} is the total number of taxa (that have non-zero counts in at least one sample).
#' @param n.rej.stop the minimum number of rejections (i.e., the permutation test
#'   statistic exceeds the observed test statistic) to obtain before stopping the permutation procedure. The
#'   default is 100.
#' @param n.cores the number of cores to be used for parallel computing. The default is 1.
#' @param filter.thresh a real value between 0 and 1; OTUs present in fewer than \code{filter.thresh} * 100\% samples are filtered out. The default is 0.2.
#' @param verbose a logical value indicating whether to generate verbose output during the permutation process. Default is TRUE.
#' @return A list consisting of
#' \itemize{
#'   \item effect.size.wc - effect size at each OTU, i.e., beta_j,1 - median_j'=1,...J (beta_j',1), based on the count weight
#'   \item p.taxa.wc - p-values by the OTU-specific tests based on the count weight
#'   \item q.taxa.wc - q-values (adjusted p-values by the BH procedure) by the OTU-specific tests based on the count weight
#'   \item detected.taxa.wc - detected taxa (using the column names of the OTU table) at the nominal FDR by the OTU-specific tests based on the count weight
#'   \item p.global.wc - p-value by the global test based on the count weight
#'   \item effect.size.we - effect size at each OTU, i.e., beta_j,1 - median_j'=1,...J(beta_j',1), based on the equal weight
#'   \item p.taxa.we - p-values by the OTU-specific tests based on the equal weight
#'   \item q.taxa.we - q-values (adjusted p-values by the BH procedure) by the OTU-specific tests based on the equal weight
#'   \item detected.taxa.we - detected taxa (using the column names of the OTU table) at the nominal FDR by the OTU-specific tests based on the equal weight
#'   \item p.global.we - p-value by the global test based on the equal weight
#'   \item p.taxa.omni - p-values by the OTU-specific omnibus tests 
#'   \item q.taxa.omni - q-values (adjusted p-values by the BH procedure) by the OTU-specific omnibus tests 
#'   \item detected.taxa.omni - detected taxa (using the column names of the OTU table) at the nominal FDR by the OTU-specific omnibus tests 
#'   \item p.global.omni - p-value by the global omnibus test 
#'   \item n.perm.completed - number of permutations completed
#'   \item ref.otu - the selected reference OTU
#'   \item seed - the seed used to generate the permutation replicates
#' }
#' @export
#' @examples
#' library(LOCOM2)
#' data("otu.table.1")
#' data("otu.table.2")
#' data("meta.1")
#' data("meta.2")
#' 
#' # restricting to a subset of samples and OTUs for speed
#' sub.sam.1 <- c(31:60)
#' sub.sam.2 <- c(1:20, 51:60)
#' otu.table.sub.1 <- otu.table.1[sub.sam.1, 1:70]
#' otu.table.sub.2 <- otu.table.2[sub.sam.2, 1:70]
#' meta.sub.1 <- meta.1[sub.sam.1,]
#' meta.sub.2 <- meta.2[sub.sam.2,]
#' 
#' # running com2seq
#' res.com2seq <- com2seq(formula.1 = otu.table.sub.1 | C ~ as.factor(Y), data.1 = meta.sub.1, 
#'                        formula.2 = otu.table.sub.2 | C ~ as.factor(Y), data.2 = meta.sub.2, 
#'                        seed = 123, n.cores = 1, n.perm.max = 1000)
#' res.com2seq$p.global.omni 
#' res.com2seq$detected.taxa.omni 


com2seq <- function(formula.1, formula.2, data.1 = NULL, data.2 = NULL,
                    fdr.nominal = 0.2, seed = NULL,
                    n.perm.max = NULL, n.rej.stop = 100, n.cores = 1,
                    filter.thresh = 0.2, verbose = TRUE) {
    
    ## ================================================================
    ## Process formula.1
    ## ================================================================
    
    if (!inherits(formula.1, "formula")) {
        stop("'formula' must be a formula.")
    }
    
    env.1 <- environment(formula.1)
    
    ## Helper function for evaluating expressions
    eval_expr <- function(expr) {
        if (is.null(data.1)) {
            eval(expr, envir = env.1)
        } else {
            eval(expr, envir = data.1, enclos = env.1)
        }
    }
    
    ## ------------------------------------------------
    ## 1. Parse formula.1
    ## ------------------------------------------------
    
    lhs.1 <- formula.1[[2]]
    rhs.1 <- formula.1[[3]]
    
    ## Formula can be:
    ##     otu.table ~ Y
    ##     otu.table | C ~ Y
    ##     otu.table | age + sex + BMI ~ Y
    
    if (is.call(lhs.1) && identical(lhs.1[[1]], as.name("|"))) {
        
        otu.expr.1 <- lhs.1[[2]]
        C.expr.1   <- lhs.1[[3]]
        
    } else {
        
        otu.expr.1 <- lhs.1
        C.expr.1   <- NULL
    }
    
    
    ## ------------------------------------------------
    ## 2. Extract OTU table1
    ## ------------------------------------------------
    
    table1 <- eval_expr(otu.expr.1)
    
    if (!is.matrix(table1) && !is.data.frame(table1)) {
        stop("'otu.table' must be a matrix or data frame.")
    }
    
    n.1 <- nrow(table1)
    
    
    ## ------------------------------------------------
    ## 3. Process covariates C1
    ## ------------------------------------------------
    
    if (is.null(C.expr.1)) {
        
        C1 <- NULL
        
    } else {
        
        ## First check whether C.expr is the name of an existing
        ## matrix/data.frame, e.g. otu.table | C ~ Y
        
        C.raw.1 <- try(eval_expr(C.expr.1), silent = TRUE)
        
        if (is.symbol(C.expr.1) &&
            !inherits(C.raw.1, "try-error") &&
            (is.matrix(C.raw.1) || is.data.frame(C.raw.1))) {
            
            C.raw.1 <- as.data.frame(C.raw.1)
            
            ## Convert factors to dummy variables
            C1 <- model.matrix(~ ., data = C.raw.1)
            
        } else {
            
            ## General case, e.g., otu.table | age + sex + BMI ~ Y
            
            C.formula.1 <- as.formula(call("~", C.expr.1), env = env.1)
            
            C1 <- model.matrix(C.formula.1, data = data.1)
        }
        
        ## Remove intercept
        C1 <- C1[, colnames(C1) != "(Intercept)", drop = FALSE]
        
        if (nrow(C1) != n.1) {
            stop("The number of observations in 'C' does not match ",
                 "the number of rows of 'otu.table'.")
        }
    }
    
    
    ## ------------------------------------------------
    ## 4. Process trait Y1
    ## ------------------------------------------------
    
    Y.raw.1 <- try(eval_expr(rhs.1), silent = TRUE)
    
    if (is.symbol(rhs.1) &&
        !inherits(Y.raw.1, "try-error") &&
        (is.matrix(Y.raw.1) || is.data.frame(Y.raw.1))) {
        
        Y.raw.1 <- as.data.frame(Y.raw.1)
        Y1 <- model.matrix(~ ., data = Y.raw.1)
        
    } else {
        
        Y.formula.1 <- as.formula(call("~", rhs.1), env = env.1)
        Y1 <- model.matrix(Y.formula.1, data = data.1)
    }
    
    ## Remove intercept
    Y1 <- Y1[, colnames(Y1) != "(Intercept)", drop = FALSE]
    
    if (nrow(Y1) != n.1) {
        stop("The number of observations in 'Y' does not match ",
             "the number of rows of 'otu.table'.")
    }
    
    if (ncol(Y1) == 0L) {
        stop("'Y' contains no estimable trait columns.")
    }
    
    
    ## ================================================================
    ## Process formula.2
    ## ================================================================
    
    if (!inherits(formula.2, "formula")) {
        stop("'formula' must be a formula.")
    }
    
    env.2 <- environment(formula.2)
    
    ## Helper function for evaluating expressions
    eval_expr <- function(expr) {
        if (is.null(data.2)) {
            eval(expr, envir = env.2)
        } else {
            eval(expr, envir = data.2, enclos = env.2)
        }
    }
    
    ## ------------------------------------------------
    ## 1. Parse formula.2
    ## ------------------------------------------------
    
    lhs.2 <- formula.2[[2]]
    rhs.2 <- formula.2[[3]]
    
    ## Formula can be:
    ##     otu.table ~ Y
    ##     otu.table | C ~ Y
    ##     otu.table | age + sex + BMI ~ Y
    
    if (is.call(lhs.2) && identical(lhs.2[[1]], as.name("|"))) {
        
        otu.expr.2 <- lhs.2[[2]]
        C.expr.2   <- lhs.2[[3]]
        
    } else {
        
        otu.expr.2 <- lhs.2
        C.expr.2   <- NULL
    }
    
    
    ## ------------------------------------------------
    ## 2. Extract OTU table2
    ## ------------------------------------------------
    
    table2 <- eval_expr(otu.expr.2)
    
    if (!is.matrix(table2) && !is.data.frame(table2)) {
        stop("'otu.table' must be a matrix or data frame.")
    }
    
    n.2 <- nrow(table2)
    
    
    ## ------------------------------------------------
    ## 3. Process covariates C2
    ## ------------------------------------------------
    
    if (is.null(C.expr.2)) {
        
        C2 <- NULL
        
    } else {
        
        ## First check whether C.expr is the name of an existing
        ## matrix/data.frame, e.g. otu.table | C ~ Y
        
        C.raw.2 <- try(eval_expr(C.expr.2), silent = TRUE)
        
        if (is.symbol(C.expr.2) &&
            !inherits(C.raw.2, "try-error") &&
            (is.matrix(C.raw.2) || is.data.frame(C.raw.2))) {
            
            C.raw.2 <- as.data.frame(C.raw.2)
            
            ## Convert factors to dummy variables
            C2 <- model.matrix(~ ., data = C.raw.2)
            
        } else {
            
            ## General case, e.g., otu.table | age + sex + BMI ~ Y
            
            C.formula.2 <- as.formula(call("~", C.expr.2), env = env.2)
            
            C2 <- model.matrix(C.formula.2, data = data.2)
        }
        
        ## Remove intercept
        C2 <- C2[, colnames(C2) != "(Intercept)", drop = FALSE]
        
        if (nrow(C2) != n.2) {
            stop("The number of observations in 'C' does not match ",
                 "the number of rows of 'otu.table'.")
        }
    }
    
    
    ## ------------------------------------------------
    ## 4. Process trait Y2
    ## ------------------------------------------------
    
    Y.raw.2 <- try(eval_expr(rhs.2), silent = TRUE)
    
    if (is.symbol(rhs.2) &&
        !inherits(Y.raw.2, "try-error") &&
        (is.matrix(Y.raw.2) || is.data.frame(Y.raw.2))) {
        
        Y.raw.2 <- as.data.frame(Y.raw.2)
        Y2 <- model.matrix(~ ., data = Y.raw.2)
        
    } else {
        
        Y.formula.2 <- as.formula(call("~", rhs.2), env = env.2)
        Y2 <- model.matrix(Y.formula.2, data = data.2)
    }
    
    ## Remove intercept
    Y2 <- Y2[, colnames(Y2) != "(Intercept)", drop = FALSE]
    
    if (nrow(Y2) != n.2) {
        stop("The number of observations in 'Y' does not match ",
             "the number of rows of 'otu.table'.")
    }
    
    if (ncol(Y2) == 0L) {
        stop("'Y' contains no estimable trait columns.")
    }
    
    
    #########################
    ## Process input datasets
    #########################
    
    if (is.data.frame(table1)) table1 <- data.matrix(table1)
    if (is.data.frame(table2)) table2 <- data.matrix(table2)
    if (is.data.frame(Y1)) Y1 <- data.matrix(Y1)
    if (is.data.frame(Y2)) Y2 <- data.matrix(Y2)
    if (is.data.frame(C1)) C1 <- data.matrix(C1)
    if (is.data.frame(C2)) C2 <- data.matrix(C2)
    
    n.trait <- ncol(Y1)
    if (ncol(Y2) != n.trait) {
        stop("The trait designs in the two datasets do not match.")
    }
    if (xor(is.null(C1), is.null(C2))) {
        stop("Confounders must be specified for both datasets or neither.")
    }
    if (!is.null(C1)) {
        n.conf <- ncol(C1)
        if (ncol(C2) != n.conf) {
            stop("The confounder designs in the two datasets do not match.")
        }
    }
    
    # detect whether the data are relative abundances (which can be expressed as percentages have row sums around 100)
    
    if (mean(rowSums(table1)) < 200 | mean(rowSums(table2)) < 200) {
        we.only = TRUE
    } else {
        we.only = FALSE
    }

    # sort samples in the two tables
    
    if (is.null(rownames(table1)) | is.null(rownames(table2)) | is.null(colnames(table1)) | is.null(colnames(table2))) {
        stop("Column and row names are required for both OTU tables.", call. = FALSE)
    }
    
    n.sam1 <- nrow(table1)
    n.sam2 <- nrow(table2)
    sam.common <- intersect(rownames(table1), rownames(table2))
    n.sam.common <- length(sam.common)
    if (n.sam.common==0) {
        stop("No common samples between the two OTU tables. Analyze each OTU table using LOCOM2 separately.", call. = FALSE)
    }
    
    n.sam1.only <- n.sam1 - n.sam.common
    n.sam2.only <- n.sam2 - n.sam.common
    
    sam1.sort <- sam.common
    sam2.sort <- sam.common
    if (n.sam1.only > 0) sam1.sort <- c(setdiff(rownames(table1), sam.common), sam.common)
    if (n.sam2.only > 0) sam2.sort <- c(sam.common, setdiff(rownames(table2), sam.common))
    
    mat1 <- match(sam1.sort, rownames(table1))
    mat2 <- match(sam2.sort, rownames(table2))
    
    table1 <- table1[mat1, ] 
    table2 <- table2[mat2, ]
    
    Y1 <- Y1[mat1,,drop=FALSE]
    Y2 <- Y2[mat2,,drop=FALSE]
    Y <- rbind(Y1, Y2) 
    
    if (is.null(C1) | is.null(C2)) {
        C <- NULL
    } else {
            C1 <- C1[mat1,,drop=FALSE]
            C2 <- C2[mat2,,drop=FALSE]
            C <- rbind(C1, C2)
    }
    
    # create three tables (table.comTaxa.filter, table1.only.filter, table2.only.filter)
    
    taxa.common <- intersect(colnames(table1), colnames(table2))
 
    table1.filter <- table1[ ,which(colMeans(table1 > 0) >= filter.thresh), drop = FALSE] # filter 1
    table2.filter <- table2[ ,which(colMeans(table2 > 0) >= filter.thresh), drop = FALSE] # filter 2
    
    taxa1.only.filter <- setdiff(colnames(table1.filter), taxa.common)
    taxa2.only.filter <- setdiff(colnames(table2.filter), taxa.common)
    n.taxa1.only.filter <- length(taxa1.only.filter)
    n.taxa2.only.filter <- length(taxa2.only.filter)
    if (n.taxa1.only.filter > 0) {
        table1.only.filter <- table1.filter[, taxa1.only.filter, drop=FALSE]
    }
    if (n.taxa2.only.filter > 0) {
        table2.only.filter <- table2.filter[, taxa2.only.filter, drop=FALSE]
    }
    rm(taxa1.only.filter, taxa2.only.filter)
    
    sam.union <- union(rownames(table1), rownames(table2))
    table1.uniSam.comTaxa <- matrix(0, nrow = length(sam.union), ncol = length(taxa.common))
    table2.uniSam.comTaxa <- matrix(0, nrow = length(sam.union), ncol = length(taxa.common))
    table1.uniSam.comTaxa[match(rownames(table1), sam.union), ] <- table1[, match(taxa.common, colnames(table1))]
    table2.uniSam.comTaxa[match(rownames(table2), sam.union), ] <- table2[, match(taxa.common, colnames(table2))]
    table.uniSam.comTaxa <- table1.uniSam.comTaxa + table2.uniSam.comTaxa
    
    taxa.common.filter3 <- taxa.common[which(colMeans(table.uniSam.comTaxa > 0) >= filter.thresh)] # filter 3
    
    rm(table1.uniSam.comTaxa, table2.uniSam.comTaxa, table.uniSam.comTaxa, sam.union)
    
    taxa.3filters <- unique(c(colnames(table1.filter), colnames(table2.filter), taxa.common.filter3)) # passing three filters
    taxa.common.filter <- intersect(taxa.common, taxa.3filters)
    
    table.comTaxa.filter <- rbind(table1[, match(taxa.common.filter, colnames(table1)), drop = FALSE], 
                                  table2[, match(taxa.common.filter, colnames(table2)), drop = FALSE]) 
    
    n.sam.com <- nrow(table.comTaxa.filter)
    n.taxa.common.filter <- ncol(table.comTaxa.filter)
    if (n.taxa.common.filter == 0L) {
        stop("No common OTUs remain after filtering.", call. = FALSE)
    }
    
    # find reference OTU
    
    ref.otu <- NULL
    if (is.null(ref.otu)) {
        freq.table.comTaxa.filter <- table.comTaxa.filter/(rowSums(table.comTaxa.filter)+1)
        mean.freq <- colMeans(freq.table.comTaxa.filter)
        sd.freq <- colSds(freq.table.comTaxa.filter)
        ref.otu <- which.max(mean.freq - sd.freq)
        rm(freq.table.comTaxa.filter)
    }
    impute <- min(table.comTaxa.filter[table.comTaxa.filter[, ref.otu]>0, ref.otu])/2
    table.comTaxa.filter[table.comTaxa.filter[, ref.otu]==0, ref.otu] <- impute
    
    # -----------------------
    # Observed statistic
    # -----------------------
    
    int1 <- c(rep(1, n.sam1), rep(0, n.sam2))
    int2 <- c(rep(0, n.sam1), rep(1, n.sam2))
    
    if (is.null(C) == TRUE) {
        X <- cbind(Y, int1, int2)
        Yr <- resid(lm(Y ~ -1 + int1 + int2))
        Yr <- matrix(Yr, ncol=length(Yr)/nrow(Y))
        
        X1 <- cbind(Y1, 1)
        Y1r <- Y1
        
        X2 <- cbind(Y2, 1)
        Y2r <- Y2
    } else {
        X <- cbind(Y, int1, int2, C)
        Yr <- resid(lm(Y ~ -1 + int1 + int2 + C))
        Yr <- matrix(Yr, ncol=length(Yr)/nrow(Y))
        
        X1 <- cbind(Y1, 1, C1)
        Y1r <- resid(lm(Y1 ~ C1))
        Y1r <- matrix(Y1r, ncol=length(Y1r)/nrow(Y1))
        
        X2 <- cbind(Y2, 1, C2)
        Y2r <- resid(lm(Y2 ~ C2))
        Y2r <- matrix(Y2r, ncol=length(Y2r)/nrow(Y2))
    }
    
    freq.table.ref <- table.comTaxa.filter/(table.comTaxa.filter + table.comTaxa.filter[, ref.otu])
    n.X <- ncol(X)
    XX <- CalculateXX0(X)
    if (!we.only) weight.count = (table.comTaxa.filter + table.comTaxa.filter[, ref.otu])
    weight.equal = matrix(rep(colMeans(table.comTaxa.filter + table.comTaxa.filter[, ref.otu]), n.sam.com), ncol=n.taxa.common.filter, byrow=TRUE) # temporary
    prop.presence <- colMeans(table.comTaxa.filter > 0)
    
    if (n.taxa1.only.filter > 0) {
        freq.table1.ref <- table1.only.filter/(table1.only.filter + table.comTaxa.filter[seq_len(n.sam1), ref.otu])
        n.X1 <- ncol(X1)
        XX1 <- CalculateXX0(X1)
        weight1 = (table1.only.filter + table.comTaxa.filter[seq_len(n.sam1), ref.otu])
        prop.presence1 <- colMeans(table1.only.filter>0)
    }
    if (n.taxa2.only.filter > 0) {
        freq.table2.ref <- table2.only.filter/(table2.only.filter + table.comTaxa.filter[-(seq_len(n.sam1)), ref.otu])
        n.X2 <- ncol(X2)
        XX2 <- CalculateXX0(X2)
        weight2 = (table2.only.filter + table.comTaxa.filter[-(seq_len(n.sam1)), ref.otu])
        prop.presence2 <- colMeans(table2.only.filter>0)
    }
    beta_init0 <- array(0, dim = c(n.X, n.taxa.common.filter))
    
    tol = 1e-6
    iter_max = 50
    Firth.thresh = 0.4
    
    if (!we.only) {
        res.obs.wc <- Newton(freq_table = freq.table.ref, X = X, XX = XX, 
                        beta_init = beta_init0, weight = weight.count,
                        tol = tol, iter_max = iter_max, Firth_thresh = Firth.thresh, robust_var = FALSE, 
                        prop_presence = prop.presence, get_var = FALSE)
        beta.wc <- res.obs.wc[[1]][seq_len(n.trait),,drop=FALSE] 
    }
    
    res.obs.we <- Newton(freq_table = freq.table.ref, X = X, XX = XX, 
                             beta_init = beta_init0, weight = weight.equal,
                             tol = tol, iter_max = iter_max, Firth_thresh = Firth.thresh, robust_var = FALSE, 
                             prop_presence = prop.presence, get_var = FALSE)
    beta.we <- res.obs.we[[1]][seq_len(n.trait),,drop=FALSE] 
    
    if (n.taxa1.only.filter > 0) {
        res1.obs <- Newton(freq_table = freq.table1.ref, X = X1, XX = XX1, 
                                 beta_init = array(0, dim = c(n.X1, n.taxa1.only.filter)), weight = weight1,
                                 tol = tol, iter_max = iter_max, Firth_thresh = Firth.thresh, robust_var = FALSE,
                                 prop_presence = prop.presence1, get_var = FALSE)
        if (!we.only) beta.wc <- cbind(beta.wc, res1.obs[[1]][seq_len(n.trait),,drop=FALSE] )
        beta.we <- cbind(beta.we, res1.obs[[1]][seq_len(n.trait),,drop=FALSE] )
    }
    
    if (n.taxa2.only.filter > 0) {
        res2.obs <- Newton(freq_table = freq.table2.ref, X = X2, XX = XX2, 
                                 beta_init = array(0, dim = c(n.X2, n.taxa2.only.filter)), weight = weight2,
                                 tol = tol, iter_max = iter_max, Firth_thresh = Firth.thresh, robust_var = FALSE,
                                 prop_presence = prop.presence2, get_var = FALSE)
        if (!we.only) beta.wc <- cbind(beta.wc, res2.obs[[1]][seq_len(n.trait),,drop=FALSE] )
        beta.we <- cbind(beta.we, res2.obs[[1]][seq_len(n.trait),,drop=FALSE] )
    }
    
    if (!we.only) beta.wc <- beta.wc - rowMedians(beta.wc)
    beta.we <- beta.we - rowMedians(beta.we)
    
    shrinkage <- ifelse(is.null(C), 0.5, 0.75)
    if (!we.only) beta_init_wc = rbind(matrix(0, nrow=n.trait, ncol=ncol(res.obs.wc[[1]])), shrinkage*res.obs.wc[[1]][(n.trait+1):n.X,,drop=FALSE])
    beta_init_we = rbind(matrix(0, nrow=n.trait, ncol=ncol(res.obs.we[[1]])), shrinkage*res.obs.we[[1]][(n.trait+1):n.X,,drop=FALSE])
    if (n.taxa1.only.filter > 0) beta_init1 = rbind(matrix(0, nrow=n.trait, ncol=ncol(res1.obs[[1]])), shrinkage*res1.obs[[1]][(n.trait+1):n.X1,,drop=FALSE])
    if (n.taxa2.only.filter > 0) beta_init2 = rbind(matrix(0, nrow=n.trait, ncol=ncol(res2.obs[[1]])), shrinkage*res2.obs[[1]][(n.trait+1):n.X2,,drop=FALSE])
    
    
    # -----------------------
    # Permutation
    # -----------------------
    
    if (is.null(seed)) seed = sample(1:10^6, 1)
    set.seed(seed)
    
    n.taxa <- n.taxa.common.filter + n.taxa1.only.filter + n.taxa2.only.filter # used to be n.taxa.union.filter
    
    if (is.null(n.perm.max)) n.perm.max = ceiling(n.taxa * n.rej.stop * (1/fdr.nominal))
    
    if (!we.only) {
        beta.perm.wc <- array(NA, dim = c(n.trait, n.taxa, n.perm.max))
        n.rej.taxa.left.wc <- rep(0, n.taxa)
        n.rej.taxa.right.wc <- n.rej.taxa.left.wc
    }
    beta.perm.we <- array(NA, dim = c(n.trait, n.taxa, n.perm.max))
    n.rej.taxa.left.we <- rep(0, n.taxa)
    n.rej.taxa.right.we <- n.rej.taxa.left.we
    
    pnull.taxa.both <- array(NA, dim=c(n.trait*2, n.taxa, n.perm.max))
    
    tol.eq = 10^-8
    n.perm.block <- 1000
    n.block <- ceiling(n.perm.max / n.perm.block)
    

    if (!we.only) {
        parallel.perm.wc <- function(i) {
            start <- i * n.perm.core + 1L
            end <- min((i + 1L) * n.perm.core, ncol(perm.mat))
            if (start > end) return(NULL)
            perm_Newton(freq_table = freq.table.ref, Yr = Yr, X = X, XX = XX, 
                                    beta_init = beta_init_wc, weight = weight.count,
                                    perm = perm.mat[,  start:end, drop=FALSE],
                                    tol = tol, iter_max = iter_max, Firth_thresh = Firth.thresh, robust_var = FALSE,
                                    prop_presence = prop.presence, get_var = FALSE)
        }
    }
    parallel.perm.we <- function(i) {
        start <- i * n.perm.core + 1L
        end <- min((i + 1L) * n.perm.core, ncol(perm.mat))
        if (start > end) return(NULL)
        perm_Newton(freq_table = freq.table.ref, Yr = Yr, X = X, XX = XX, 
                                beta_init = beta_init_we, weight = weight.equal,
                                perm = perm.mat[,  start:end, drop=FALSE],
                                tol = tol, iter_max = iter_max, Firth_thresh = Firth.thresh, robust_var = FALSE,
                                prop_presence = prop.presence, get_var = FALSE)
    }
    parallel.perm1 <- function(i) {
        start <- i * n.perm.core + 1L
        end <- min((i + 1L) * n.perm.core, ncol(perm.mat1))
        if (start > end) return(NULL)
        perm_Newton(freq_table = freq.table1.ref, Yr = Y1r, X = X1, XX = XX1, 
                    beta_init = beta_init1, weight = weight1,
                    perm = perm.mat1[,  start:end, drop=FALSE],
                    tol = tol, iter_max = iter_max, Firth_thresh = Firth.thresh, robust_var = FALSE,
                    prop_presence = prop.presence1, get_var = FALSE)
    } 
    parallel.perm2 <- function(i) {
        start <- i * n.perm.core + 1L
        end <- min((i + 1L) * n.perm.core, ncol(perm.mat2))
        if (start > end) return(NULL)
        perm_Newton(freq_table = freq.table2.ref, Yr = Y2r, X = X2, XX = XX2, 
                    beta_init = beta_init2, weight = weight2,
                    perm = perm.mat2[,  start:end, drop=FALSE],
                    tol = tol, iter_max = iter_max, Firth_thresh = Firth.thresh, robust_var = FALSE,
                    prop_presence = prop.presence2, get_var = FALSE)
    } 
    
    
    n.perm <- 0 # used to be n.perm.completed
    
    for (i.block in seq_len(n.block)) {
        
        # Permutation with three strata
        n_this <- min(n.perm.block, n.perm.max - n.perm)
        n.perm.core <- ceiling(n_this/n.cores)
        
        perm.mat.common <- t(permute::shuffleSet(n.sam.common, n_this, quietly = TRUE, check = FALSE)) - 1
        perm.mat <- rbind(perm.mat.common, perm.mat.common + n.sam.common)
        if (n.sam1.only > 0) {
            perm.mat1.only <- t(shuffleSet(n.sam1.only, n_this, quietly = TRUE, check = FALSE)) - 1
            perm.mat <- rbind(perm.mat1.only, (perm.mat + n.sam1.only))
        }
        if (n.sam2.only > 0) {
            perm.mat2.only <- t(shuffleSet(n.sam2.only, n_this, quietly = TRUE, check = FALSE)) - 1
            perm.mat <- rbind(perm.mat, perm.mat2.only + n.sam1.only + 2*n.sam.common)
        }
        if (n.taxa1.only.filter > 0) {
            perm.mat1 <- perm.mat.common
            if (n.sam1.only > 0) {
                perm.mat1 <- rbind(perm.mat1.only, perm.mat.common + n.sam1.only)
            }
        }
        if (n.taxa2.only.filter > 0) {
            perm.mat2 <- perm.mat.common
            if (n.sam2.only > 0) {
                perm.mat2 <- rbind(perm.mat.common, perm.mat2.only + n.sam.common)
            }
        }
        
        if (n.cores > 1) {
            
            worker.id <- 0:(min(n.cores, ceiling(n_this / n.perm.core)) - 1)
            n.worker <- length(worker.id)
            
            if (Sys.info()[['sysname']] == 'Windows') {
                if (!we.only) parallel.stat.wc = bplapply(worker.id, parallel.perm.wc, BPPARAM = SnowParam(workers=n.worker))
                parallel.stat.we = bplapply(worker.id, parallel.perm.we, BPPARAM = SnowParam(workers=n.worker))
                if (n.taxa1.only.filter > 0) parallel.stat.only.1 = bplapply(worker.id, parallel.perm1, BPPARAM = SnowParam(workers=n.worker))
                if (n.taxa2.only.filter > 0) parallel.stat.only.2 = bplapply(worker.id, parallel.perm2, BPPARAM = SnowParam(workers=n.worker))
            } else {
                if (!we.only) parallel.stat.wc = mclapply(worker.id, parallel.perm.wc, mc.cores = n.worker)
                parallel.stat.we = mclapply(worker.id, parallel.perm.we, mc.cores = n.worker)
                if (n.taxa1.only.filter > 0) parallel.stat.only.1 = mclapply(worker.id, parallel.perm1, mc.cores = n.worker)
                if (n.taxa2.only.filter > 0) parallel.stat.only.2 = mclapply(worker.id, parallel.perm2, mc.cores = n.worker)
            }
            if (!we.only) res.perm.wc <- abind(parallel.stat.wc) 
            res.perm.we <- abind(parallel.stat.we)
            if (n.taxa1.only.filter > 0) res.perm1 <- abind(parallel.stat.only.1) 
            if (n.taxa2.only.filter > 0) res.perm2 <- abind(parallel.stat.only.2) 
            
        } else {
            if (!we.only) {
                res.perm.wc <- perm_Newton(freq_table = freq.table.ref, Yr = Yr, X = X, XX = XX,
                                                        beta_init = beta_init_wc, weight = weight.count,
                                                        perm = perm.mat,
                                                        tol = tol, iter_max = iter_max, Firth_thresh = Firth.thresh, robust_var = FALSE,
                                                        prop_presence = prop.presence, get_var = FALSE)
            }
            res.perm.we <- perm_Newton(freq_table = freq.table.ref, Yr = Yr, X = X, XX = XX,
                                                   beta_init = beta_init_we, weight = weight.equal,
                                                   perm = perm.mat,
                                                   tol = tol, iter_max = iter_max, Firth_thresh = Firth.thresh, robust_var = FALSE,
                                                   prop_presence = prop.presence, get_var = FALSE)
            if (n.taxa1.only.filter > 0) res.perm1 <- perm_Newton(freq_table = freq.table1.ref, Yr = Y1r, X = X1, XX = XX1,
                                               beta_init = beta_init1, weight = weight1,
                                               perm = perm.mat1,
                                               tol = tol, iter_max = iter_max, Firth_thresh = Firth.thresh, robust_var = FALSE,
                                               prop_presence = prop.presence1, get_var = FALSE)
            if (n.taxa2.only.filter > 0) res.perm2 <- perm_Newton(freq_table = freq.table2.ref, Yr = Y2r, X = X2, XX = XX2,
                                               beta_init = beta_init2, weight = weight2,
                                               perm = perm.mat2,
                                               tol = tol, iter_max = iter_max, Firth_thresh = Firth.thresh, robust_var = FALSE,
                                               prop_presence = prop.presence2, get_var = FALSE)
        }
        
        w <- (n.perm + 1):(n.perm + n_this)
        n.perm <- n.perm + n_this
        n.perm.inv <- 1 / (n.perm+1)
        if (verbose) message("permutations completed: ", n.perm)
        
        
        if (!we.only) if (n.taxa1.only.filter > 0) res.perm.wc <- abind(res.perm.wc, res.perm1, along=2)
        if (!we.only) if (n.taxa2.only.filter > 0) res.perm.wc <- abind(res.perm.wc, res.perm2, along=2)
        if (n.taxa1.only.filter > 0) res.perm.we <- abind(res.perm.we, res.perm1, along=2)
        if (n.taxa2.only.filter > 0) res.perm.we <- abind(res.perm.we, res.perm2, along=2)
        
        for (i.trait in seq_len(n.trait)) {
            beta.perm.we[i.trait,,w] <- t(t(res.perm.we[i.trait,,]) - colMedians(res.perm.we[i.trait,,]))
        }
        diff <- sweep(beta.perm.we[,,w,drop=FALSE], c(1,2), beta.we, FUN="-")
        n.rej.taxa.equal.we <- 0.5*(abs(diff) < tol.eq) # 2 148 1000
        n.rej.taxa.left.we <- n.rej.taxa.left.we + rowSums( (diff < - tol.eq) + n.rej.taxa.equal.we, dims=2) # 2 148
        n.rej.taxa.right.we <- n.rej.taxa.right.we + rowSums( (diff > tol.eq) + n.rej.taxa.equal.we, dims=2) # 2 148
        n.rej.taxa.we <- 2*pmin(n.rej.taxa.left.we, n.rej.taxa.right.we) # 2 148
        if (n.trait >= 2) {
            pmin.taxa.we <- colMins(n.rej.taxa.we) # 148
            pnull.taxa.we <- n.perm + 0.5 - array(rowRanks(abs(beta.perm.we[,,seq_len(n.perm),drop=FALSE]), ties.method = "average", dim.=c(n.trait*n.taxa, n.perm)), dim = c(n.trait, n.taxa, n.perm)) # 2 148 1000
            pnullmin.taxa.we <- array(colMins(pnull.taxa.we, dim.=c(n.trait, n.taxa*n.perm)), dim = c(n.taxa, n.perm)) # 148 1000 
            diff <- pnullmin.taxa.we - c(pmin.taxa.we) # 148 1000 
            n.rej.taxa.we2 <- rowSums( (diff < - tol.eq) + 0.5 * ( abs(diff) < tol.eq) )
            p.taxa.we <- (n.rej.taxa.we2+1) * n.perm.inv # 1 148
        } else {
            p.taxa.we <- (n.rej.taxa.we+1) * n.perm.inv # 1 148
        }
        q.taxa.we <- fdr.Sandve(p.taxa.we) # 148
        
        if (!we.only) {
            
            for (i.trait in seq_len(n.trait)) {
                beta.perm.wc[i.trait,,w] <- t(t(res.perm.wc[i.trait,,]) - colMedians(res.perm.wc[i.trait,,]))
            }
            diff <- sweep(beta.perm.wc[,,w,drop=FALSE], c(1,2), beta.wc, FUN="-")
            n.rej.taxa.equal.wc <- 0.5*(abs(diff) < tol.eq) # 2 148 1000
            n.rej.taxa.left.wc <- n.rej.taxa.left.wc + rowSums( (diff < - tol.eq) + n.rej.taxa.equal.wc, dims=2) # 2 148
            n.rej.taxa.right.wc <- n.rej.taxa.right.wc + rowSums( (diff > tol.eq) + n.rej.taxa.equal.wc, dims=2) # 2 148
            n.rej.taxa.wc <- 2*pmin(n.rej.taxa.left.wc, n.rej.taxa.right.wc)  # 2 148
            if (n.trait >= 2) {
                pmin.taxa.wc <- colMins(n.rej.taxa.wc) # 148
                pnull.taxa.wc <- n.perm + 0.5 - array(rowRanks(abs(beta.perm.wc[,,seq_len(n.perm),drop=FALSE]), ties.method="average", dim.=c(n.trait*n.taxa, n.perm)), dim = c(n.trait, n.taxa, n.perm)) # 2 148 1000
                pnullmin.taxa.wc <- array(colMins(pnull.taxa.wc, dim.=c(n.trait, n.taxa*n.perm)), dim = c(n.taxa, n.perm))   # 148 1000 
                diff <- pnullmin.taxa.wc - c(pmin.taxa.wc) 
                n.rej.taxa.wc2 <- rowSums( (diff < - tol.eq) + 0.5 * ( abs(diff) < tol.eq) )
                p.taxa.wc <- (n.rej.taxa.wc2+1) * n.perm.inv 
            } else {
                p.taxa.wc <- (n.rej.taxa.wc+1) * n.perm.inv 
            }
            q.taxa.wc <- fdr.Sandve(p.taxa.wc) 
            
            pmin.taxa.omni <- colMins(rbind(n.rej.taxa.wc, n.rej.taxa.we)) # 148
            if (n.trait == 1) {
                pnull.taxa.we <- n.perm + 0.5 - array(rowRanks(abs(beta.perm.we[,,seq_len(n.perm),drop=FALSE]), ties.method="average", dim.=c(n.trait*n.taxa, n.perm)), dim = c(n.trait, n.taxa, n.perm)) # 2 148 1000
                pnull.taxa.wc <- n.perm + 0.5 - array(rowRanks(abs(beta.perm.wc[,,seq_len(n.perm),drop=FALSE]), ties.method="average", dim.=c(n.trait*n.taxa, n.perm)), dim = c(n.trait, n.taxa, n.perm)) # 2 148 1000
            }
            pnull.taxa.both[seq_len(n.trait),,seq_len(n.perm)] <- pnull.taxa.wc
            pnull.taxa.both[(n.trait+1):(2*n.trait),,seq_len(n.perm)] <- pnull.taxa.we # 4 148 1000
            pnullmin.taxa.omni <- array(colMins(pnull.taxa.both[,,seq_len(n.perm)], dim.=c(n.trait*2, n.taxa*n.perm)), dim = c(n.taxa, n.perm)) # 148 1000 
            diff <- pnullmin.taxa.omni- c(pmin.taxa.omni) 
            n.rej.taxa.omni <- rowSums( (diff < - tol.eq) + 0.5 * ( abs(diff) < tol.eq) )
            p.taxa.omni <- (n.rej.taxa.omni+1) * n.perm.inv
            q.taxa.omni <- fdr.Sandve(p.taxa.omni)
        } 
        
        if (n.trait >= 2) {
            n.rej.stop.we <- n.rej.taxa.we2
            if (!we.only) n.rej.stop.wc <- n.rej.taxa.wc2
        } else {
            n.rej.stop.we <- c(n.rej.taxa.we)
            if (!we.only) n.rej.stop.wc <- c(n.rej.taxa.wc)
        }
        
        if (!we.only) {
            if (n.perm >= 10000 
                & all(q.taxa.wc <= fdr.nominal | n.rej.stop.wc >= n.rej.stop)
                & all(q.taxa.we <= fdr.nominal | n.rej.stop.we >= n.rej.stop)
                & all(q.taxa.omni <= fdr.nominal | n.rej.taxa.omni >= n.rej.stop)
            ) break
        } else {
            if (n.perm >= 10000 
                & all(q.taxa.we <= fdr.nominal | n.rej.stop.we >= n.rej.stop)
            ) break
        }
        
    } # permutation
    
    taxa.union <- colnames(table.comTaxa.filter)
    if (n.taxa1.only.filter > 0) taxa.union <- c(taxa.union, colnames(table1.only.filter))
    if (n.taxa2.only.filter > 0) taxa.union <- c(taxa.union, colnames(table2.only.filter))
    
    # Multiple test correction by BH
    q.taxa.we <- p.adjust(p.taxa.we, method = "BH")
    detected.taxa.we <- taxa.union[which(q.taxa.we < fdr.nominal)]
    if (!we.only) {
        q.taxa.wc <- p.adjust(p.taxa.wc, method = "BH")
        q.taxa.omni <- p.adjust(p.taxa.omni, method = "BH")
        detected.taxa.wc <- taxa.union[which(q.taxa.wc < fdr.nominal)]
        detected.taxa.omni <- taxa.union[which(q.taxa.omni < fdr.nominal)]
    }
    
    
    # ------------------------
    # Global p-value
    # ------------------------
    
    if (n.trait == 1) {
        beta.all.we <- cbind(beta.we[1,], beta.perm.we[1,,seq_len(n.perm)])
        r <- matrixStats::rowRanks(beta.all.we, ties.method = "average") # 148 1001
        p.taxa1.we <- (2*pmin(r, n.perm+2-r) - 1.5) # omit n.perm.inv # 148 1001
    } else {
        beta.all.we <- cbind(pmin.taxa.we, pnullmin.taxa.we)  # 148 1001 
        r <- matrixStats::rowRanks(beta.all.we, ties.method = "average") # 148 1001
        p.taxa1.we <- r / (n.perm + 1)
    }

    stat.global.we <- colSums(1/p.taxa1.we)
    p.global.we <- (sum(stat.global.we[1] <= stat.global.we[-1]) + 1) * n.perm.inv
    
    if (!we.only) {
        if (n.trait == 1) {
            beta.all.wc <- cbind(beta.wc[1,], beta.perm.wc[1,,seq_len(n.perm)])
            r <- rowRanks(beta.all.wc, ties.method = "average") # 148 1001
            p.taxa1.wc <- (2*pmin(r, n.perm+2-r) - 1.5)
        } else {
            beta.all.wc <- cbind(pmin.taxa.wc, pnullmin.taxa.wc)
            r <- rowRanks(beta.all.wc, ties.method = "average") # 148 1001
            p.taxa1.wc <- r / (n.perm + 1)
        }
        stat.global.wc <- colSums(1/p.taxa1.wc)
        p.global.wc <- (sum(stat.global.wc[1] <= stat.global.wc[-1]) + 1) * n.perm.inv
        
        beta.all.omni <- cbind(pmin.taxa.omni, pnullmin.taxa.omni) # 148 1001
        p.taxa1.omni <- (rowRanks(beta.all.omni, ties.method = "average") - 0.5) # 148 1001
        stat.global.omni <- colSums(1/p.taxa1.omni) # 1001
        p.global.omni<- (sum(stat.global.omni[1] <= stat.global.omni[-1]) + 1) * n.perm.inv
    } 
    
    beta.we <- matrix(beta.we, nrow=n.trait)
    colnames(beta.we) <- taxa.union
    p.taxa.we <- matrix(c(p.taxa.we), nrow = 1)
    q.taxa.we <- matrix(c(q.taxa.we), nrow = 1)
    colnames(p.taxa.we) <- taxa.union
    colnames(q.taxa.we) <- taxa.union
    
    if (!we.only) {
        beta.wc <- matrix(beta.wc, nrow=n.trait)
        colnames(beta.wc) <- taxa.union
        p.taxa.wc <- matrix(c(p.taxa.wc), nrow = 1)
        q.taxa.wc <- matrix(c(q.taxa.wc), nrow = 1)
        colnames(p.taxa.wc) <- taxa.union
        colnames(q.taxa.wc) <- taxa.union
    
        p.taxa.omni <- matrix(c(p.taxa.omni), nrow = 1)
        q.taxa.omni <- matrix(c(q.taxa.omni), nrow = 1)
        colnames(p.taxa.omni) <- taxa.union
        colnames(q.taxa.omni) <- taxa.union
    } else {
        beta.wc <- NULL
        p.taxa.wc <- NULL
        q.taxa.wc <- NULL
        detected.taxa.wc <- NULL
        p.global.wc <- NULL
        
        p.taxa.omni <- NULL
        q.taxa.omni <- NULL
        detected.taxa.omni = NULL
        p.global.omni = NULL
    }
    
    return(list(effect.size.wc = beta.wc,
                p.taxa.wc = p.taxa.wc,
                q.taxa.wc = q.taxa.wc,
                detected.taxa.wc = detected.taxa.wc,
                p.global.wc = p.global.wc,
                
                effect.size.we = beta.we,
                p.taxa.we = p.taxa.we,
                q.taxa.we = q.taxa.we,
                detected.taxa.we = detected.taxa.we,
                p.global.we = p.global.we,
                
                p.taxa.omni = p.taxa.omni,
                q.taxa.omni = q.taxa.omni,
                detected.taxa.omni = detected.taxa.omni,
                p.global.omni = p.global.omni,
                
                n.perm.completed = n.perm,
                ref.otu = colnames(table.comTaxa.filter)[ref.otu],
                seed = seed))
} # com2seq

