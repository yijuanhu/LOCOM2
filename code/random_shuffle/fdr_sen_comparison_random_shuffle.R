# Random-shuffle (null) simulation of the LOCOM2 paper: the real count tables of the throat, vaginal
# and IBD templates are kept unchanged and the binary trait Y is assigned at random, so every taxon
# is null. The raw taxon-level p-values of LOCOM2, LOCOM, MaAsLin3, the MaAsLin2-style log-linear
# model, ANCOM-BC2 and LinDA are saved for every shuffle; the empirical type I error is computed
# from them afterwards.

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 7) {
  stop(
    "Usage: Rscript fdr_sen_comparison_random_shuffle.R ",
    "<start.seed> <n.rep> <sim.method> <n.perm.max> ",
    "<prop.case> <Firth.thresh> <n.wald>"
  )
}

start.seed   <- as.integer(args[1])
n.rep        <- as.integer(args[2])
sim.method   <- as.integer(args[3])
n.perm.max   <- as.integer(args[4])
prop.case    <- as.numeric(args[5])
Firth.thresh <- as.numeric(args[6])
n.wald       <- as.integer(args[7])

if (!is.finite(start.seed)) stop("start.seed must be finite.")
if (!is.finite(n.rep) || n.rep < 1) stop("n.rep must be >= 1.")
if (!(sim.method %in% 1:4)) stop("sim.method must be 1, 2, 3, or 4.")
if (!is.finite(prop.case) || prop.case <= 0 || prop.case >= 1) {
  stop("prop.case must be strictly between 0 and 1.")
}

library(RcppArmadillo)
library(Rcpp)
library(MASS)
library(permute)
library(parallel)
library(abind)
library(matrixStats)
library(maaslin3)
library(S4Vectors)
library(TreeSummarizedExperiment)
library(ANCOMBC)
library(LinDA)

# folder with the template count tables and the LOCOM2 / LOCOM code
data.dir <- "../data_templates"

source(file.path(data.dir, "LOCOM2_func.R"))
sourceCpp(file.path(data.dir, "Newton_LOCOM2.cpp"))

# original LOCOM
source(file.path(data.dir, "LOCOM_fun.R"))
sourceCpp(file.path(data.dir, "Newton.cpp"))

# the four template count tables (throat.otu.tab, count.vaginal, count.female, count.female.filtered);
# only Y is randomized across replicates
load(file.path(data.dir, "all.templates.RData"))

# method-internal settings (the script saves raw p-values; the nominal level for the empirical
# type I error is chosen later). LOCOM / LOCOM2 use fdr.nominal and n.rej.stop in their adaptive
# permutation procedure, so the values of the main simulations are kept.
internal.fdr.nominal <- 0.2
internal.n.rej.stop <- 100
internal.linda.alpha <- 0.2

# -----------------------------------------------------------------------------
# Template selection
# -----------------------------------------------------------------------------
get_template <- function(sim.method) {
  if (sim.method == 1) {
    list(
      otu.table = throat.otu.tab,
      ref.otu = "4414",
      name = "throat"
    )
  } else if (sim.method == 2) {
    list(
      otu.table = count.vaginal,
      ref.otu = NULL,
      name = "vaginal"
    )
  } else if (sim.method == 3) {
    list(
      otu.table = count.female,
      ref.otu = NULL,
      name = "ibd_female"
    )
  } else if (sim.method == 4) {
    list(
      otu.table = count.female.filtered,
      ref.otu = NULL,
      name = "ibd_female_filtered"
    )
  } else {
    stop("sim.method must be 1, 2, 3, or 4.")
  }
}

# -----------------------------------------------------------------------------
# Null simulation: the template count table is kept unchanged and Y is assigned at random
# -----------------------------------------------------------------------------
simulate_data <- function(otu.template, ref.otu, seed, prop.case,
                          filter.thresh) {
  n.sam <- nrow(otu.template)
  
  if (n.sam < 2) stop("Template must contain at least two samples.")
  if (is.null(colnames(otu.template))) {
    colnames(otu.template) <- paste0("OTU", seq_len(ncol(otu.template)))
  }
  
  lib.size <- rowSums(otu.template)
  if (any(!is.finite(lib.size)) || any(lib.size <= 0)) {
    stop("Template contains at least one sample with non-positive library size.")
  }
  
  n1.sam <- round(n.sam * prop.case)
  n0.sam <- n.sam - n1.sam
  
  if (n1.sam < 1 || n0.sam < 1) {
    stop("prop.case produces an empty case or control group for this template.")
  }
  
  set.seed(seed)
  Y <- sample(c(rep(0L, n0.sam), rep(1L, n1.sam)))
  
  otu.table.sim <- otu.template
  
  # keep taxa present in at least filter.thresh * n.sam samples (independent of the shuffled labels)
  prop.count <- colSums(otu.table.sim > 0)
  otus.keep <- which(prop.count >= filter.thresh * n.sam)
  
  if (length(otus.keep) == 0) {
    stop("No taxa remain after prevalence filtering.")
  }
  
  otu.table.sim.filter <- otu.table.sim[, otus.keep, drop = FALSE]
  
  if (is.null(ref.otu)) {
    ref.otu.filter <- NULL
  } else {
    ref.otu.filter <- match(ref.otu, colnames(otu.table.sim.filter))
    if (is.na(ref.otu.filter)) ref.otu.filter <- NULL
  }
  
  list(
    ref.otu.filter = ref.otu.filter,
    otu.table.sim.filter = otu.table.sim.filter,
    Y = Y
  )
}

# -----------------------------------------------------------------------------
# Utility functions for raw taxon-level p-values
# -----------------------------------------------------------------------------
make_na_p <- function(target_taxa) {
  out <- rep(NA_real_, length(target_taxa))
  names(out) <- target_taxa
  out
}

align_values <- function(value, taxa, target_taxa) {
  out <- make_na_p(target_taxa)
  
  if (is.null(value) || is.null(taxa)) return(out)
  
  idx <- match(target_taxa, as.character(taxa))
  keep <- !is.na(idx)
  out[keep] <- as.numeric(value[idx[keep]])
  out
}

align_p_by_length <- function(p, target_taxa, method.name) {
  if (is.null(p)) return(make_na_p(target_taxa))
  
  # keep the taxon labels of a named vector or of a 1 x J / J x 1 matrix (LOCOM / LOCOM2 return 1 x J matrices)
  if (!is.null(names(p))) {
    return(align_values(p, names(p), target_taxa))
  }
  
  if (is.matrix(p) || is.data.frame(p)) {
    if (nrow(p) == 1 && !is.null(colnames(p))) {
      return(align_values(as.numeric(p[1, ]), colnames(p), target_taxa))
    }
    if (ncol(p) == 1 && !is.null(rownames(p))) {
      return(align_values(as.numeric(p[, 1]), rownames(p), target_taxa))
    }
  }
  
  p <- as.numeric(p)
  if (length(p) != length(target_taxa)) {
    warning(
      method.name,
      " returned ", length(p),
      " p-values for ", length(target_taxa),
      " retained taxa; treating this replicate as unavailable for that method."
    )
    return(make_na_p(target_taxa))
  }
  
  names(p) <- target_taxa
  p
}

# extract a p-value vector from a result list by trying a short list of component names
extract_p_component <- function(res, candidates, target_taxa, method.name) {
  if (inherits(res, "try-error") || is.null(res)) {
    return(make_na_p(target_taxa))
  }
  
  hit <- candidates[vapply(
    candidates,
    function(nm) !is.null(res[[nm]]),
    logical(1)
  )]
  
  if (length(hit) == 0) {
    warning(
      method.name,
      " did not return any expected raw p-value component: ",
      paste(candidates, collapse = ", "),
      ". Available components are: ",
      paste(names(res), collapse = ", ")
    )
    return(make_na_p(target_taxa))
  }
  
  align_p_by_length(res[[hit[1]]], target_taxa, method.name)
}

extract_locom_pvalue <- function(res, target_taxa) {
  extract_p_component(
    res = res,
    candidates = c("p.otu", "p.otu.wr"),
    target_taxa = target_taxa,
    method.name = "LOCOM"
  )
}

extract_ancombc2_pvalue <- function(res, target_taxa) {
  if (inherits(res, "try-error") || is.null(res$res)) {
    return(make_na_p(target_taxa))
  }
  
  tab <- res$res
  taxa <- if ("taxon" %in% colnames(tab)) {
    as.character(tab$taxon)
  } else {
    rownames(tab)
  }
  
  if (!"p_Y" %in% colnames(tab)) {
    warning("ANCOM-BC2 result does not contain raw p-value column p_Y.")
    return(make_na_p(target_taxa))
  }
  
  # the raw p-value itself is used; passed_ss_Y belongs to ANCOM-BC2's discovery rule, not to the
  # taxon-wise type I error
  align_values(tab$p_Y, taxa, target_taxa)
}

extract_linda_pvalue <- function(res, target_taxa) {
  if (inherits(res, "try-error") ||
      is.null(res$output) || length(res$output) == 0) {
    return(make_na_p(target_taxa))
  }
  
  if ("Y" %in% names(res$output)) {
    tab <- res$output[["Y"]]
  } else {
    y.idx <- grep("Y", names(res$output))
    tab <- if (length(y.idx) > 0) res$output[[y.idx[1]]] else res$output[[1]]
  }
  
  if (is.null(tab) || !"pvalue" %in% colnames(tab)) {
    return(make_na_p(target_taxa))
  }
  
  align_values(tab$pvalue, rownames(tab), target_taxa)
}

# -----------------------------------------------------------------------------
# Prepare fixed template information
# -----------------------------------------------------------------------------
template.info <- get_template(sim.method)
raw.template <- template.info$otu.table
n.sam <- nrow(raw.template)
filter.thresh <- min(0.1, 10 / n.sam)
method.name <- template.info$name

message(
  "Null simulation template: ", method.name,
  "; n = ", n.sam,
  "; prop.case = ", prop.case,
  "; prevalence threshold = ", signif(filter.thresh, 4)
)

# -----------------------------------------------------------------------------
# Output file and simulation loop
# -----------------------------------------------------------------------------
# one row per (shuffle seed, taxon) with one raw p-value column per method, so that the taxon-specific
# empirical type I error can be computed later at any nominal level
raw.pvalue.file <- paste0(
  "null_raw_pvalues",
  "_thresh", signif(filter.thresh, 4),
  "_sim_", start.seed,
  "_nrep_", n.rep,
  "_n_", n.sam,
  "_template_", method.name,
  "_pcase_", prop.case,
  "_nperm_", n.perm.max,
  "_firth_", Firth.thresh,
  "_nwald_", n.wald,
  ".csv"
)

if (file.exists(raw.pvalue.file)) {
  file.remove(raw.pvalue.file)
}

file.initialized <- FALSE
n.completed <- 0L

for (seed.current in seq.int(start.seed, length.out = n.rep)) {
  
  simulated <- try(
    simulate_data(
      otu.template = raw.template,
      ref.otu = template.info$ref.otu,
      seed = seed.current,
      prop.case = prop.case,
      filter.thresh = filter.thresh
    ),
    silent = TRUE
  )
  
  if (inherits(simulated, "try-error")) {
    message("Seed ", seed.current, ": null shuffling failed: ", simulated)
    next
  }
  
  n.sam.current <- nrow(simulated$otu.table.sim.filter)
  feature.names <- colnames(simulated$otu.table.sim.filter)
  
  if (is.null(feature.names)) {
    feature.names <- paste0(
      "OTU", seq_len(ncol(simulated$otu.table.sim.filter))
    )
    colnames(simulated$otu.table.sim.filter) <- feature.names
  }
  
  if (anyDuplicated(feature.names)) {
    stop("Filtered OTU table contains duplicated taxon names; safe p-value alignment is impossible.")
  }
  
  # ---------------------------------------------------------------------------
  # LOCOM2: raw taxon-level p-values
  # ---------------------------------------------------------------------------
  res.locom2 <- try(
    locom2(
      otu.table = simulated$otu.table.sim.filter,
      Y = simulated$Y,
      C = NULL,
      fdr.nominal = internal.fdr.nominal,
      ref.otu = simulated$ref.otu.filter,
      n.perm.max = n.perm.max,
      n.rej.stop = internal.n.rej.stop,
      n_cores = 1,
      filter.thresh = 0,
      seed = seed.current,
      n.wald = n.wald,
      Firth.thresh = Firth.thresh
    ),
    silent = TRUE
  )
  
  locom2.p.components <- list(
    "LOCOM-count-wr" = c("p.otu.wr"),
    "LOCOM-Wald-count-wr" = c("p.otu.wald.wr"),
    "LOCOM-Wald-count-wr-yj" = c("p.otu.wald.wr.yj"),
    "LOCOM-Wald-count-wr-sig" = c("p.otu.wald.wr.sig"),
    "LOCOM-Wald-count-wr-pval" = c("p.otu.wald.wr.pval")
  )
  
  if (inherits(res.locom2, "try-error")) {
    pvalue.list <- setNames(
      replicate(
        length(locom2.p.components),
        make_na_p(feature.names),
        simplify = FALSE
      ),
      names(locom2.p.components)
    )
    message("Seed ", seed.current, ": LOCOM2 failed")
  } else {
    pvalue.list <- lapply(
      names(locom2.p.components),
      function(method.current) {
        extract_p_component(
          res = res.locom2,
          candidates = locom2.p.components[[method.current]],
          target_taxa = feature.names,
          method.name = method.current
        )
      }
    )
    names(pvalue.list) <- names(locom2.p.components)
    message("Seed ", seed.current, ": LOCOM2 finished")
  }
  
  # ---------------------------------------------------------------------------
  # Comparison methods (metadata: Y only)
  # ---------------------------------------------------------------------------
  sample.names <- paste0("sam", seq_len(n.sam.current))
  rownames(simulated$otu.table.sim.filter) <- sample.names
  
  comparison.meta <- data.frame(
    Y = simulated$Y,
    row.names = sample.names
  )
  
  maaslin.formula <- "~ Y"
  ancom.formula <- "Y"
  linda.formula <- "~ Y"
  
  # ---------------------------------------------------------------------------
  # LOCOM
  # ---------------------------------------------------------------------------
  res.locom <- try(
    locom(
      otu.table = simulated$otu.table.sim.filter,
      Y = simulated$Y,
      C = NULL,
      fdr.nominal = internal.fdr.nominal,
      n.cores = 1,
      n.perm.max = n.perm.max,
      filter.thresh = 0,
      seed = seed.current
    ),
    silent = TRUE
  )
  
  p.locom <- extract_locom_pvalue(res.locom, feature.names)
  message(
    "Seed ", seed.current, ": LOCOM ",
    if (inherits(res.locom, "try-error")) "failed" else "finished"
  )
  
  # ---------------------------------------------------------------------------
  # MaAsLin3
  # ---------------------------------------------------------------------------
  res.maaslin3 <- try({
    read_data_list <- maaslin_read_data(
      input_data = simulated$otu.table.sim.filter,
      input_metadata = comparison.meta,
      feature_specific_covariate = NULL,
      unscaled_abundance = NULL
    )
    
    read_data_list <- maaslin_reorder_data(
      read_data_list$data,
      read_data_list$metadata,
      read_data_list$feature_specific_covariate,
      read_data_list$unscaled_abundance
    )
    
    data <- read_data_list$data
    metadata <- read_data_list$metadata
    
    formulas <- maaslin_check_formula(
      data,
      metadata,
      maaslin.formula,
      feature_specific_covariate_name = NULL
    )
    
    norm <- maaslin_normalize(data, output = NULL)
    filt <- maaslin_filter(norm, output = NULL)
    trans <- maaslin_transform(filt, output = NULL)
    
    metadata <- maaslin_process_metadata(
      metadata = metadata,
      formula = formulas$formula,
      standardize = FALSE,
      fixed_effects = NULL,
      reference = "Y,0",
      feature_specific_covariate_name = NULL
    )
    
    invisible(
      capture.output(
        tmp.res <- suppressMessages(
          maaslin_fit(
            filtered_data = filt,
            transformed_data = trans,
            metadata = metadata,
            formula = formulas$formula,
            random_effects_formula = NULL,
            data = data
          )
        )
      )
    )
    
    tmp.res
  }, silent = TRUE)
  
  if (inherits(res.maaslin3, "try-error")) {
    p.maaslin3 <- make_na_p(feature.names)
    message("Seed ", seed.current, ": MaAsLin3 failed")
  } else {
    maaslin3.tab <- res.maaslin3$fit_data_abundance$results
    
    if (is.null(maaslin3.tab) ||
        !("metadata" %in% colnames(maaslin3.tab)) ||
        !("pval_joint" %in% colnames(maaslin3.tab)) ||
        !("feature" %in% colnames(maaslin3.tab))) {
      p.maaslin3 <- make_na_p(feature.names)
      message("Seed ", seed.current, ": MaAsLin3 returned an unexpected result format")
    } else {
      maaslin3.tab <- maaslin3.tab[
        maaslin3.tab$metadata == "Y",
        ,
        drop = FALSE
      ]
      
      p.maaslin3 <- align_values(
        value = maaslin3.tab$pval_joint,
        taxa = maaslin3.tab$feature,
        target_taxa = feature.names
      )
      message("Seed ", seed.current, ": MaAsLin3 finished")
    }
  }
  
  # ---------------------------------------------------------------------------
  # MaAsLin2-style log-linear model: relative abundances (TSS) + half of the smallest non-zero
  # value of each taxon as pseudo-count, log-transformed and analyzed with a linear model
  # ---------------------------------------------------------------------------
  p.maaslin2 <- try({
    lib.size <- rowSums(simulated$otu.table.sim.filter)
    if (any(lib.size <= 0)) stop("Zero library size encountered.")
    
    rel.tab <- simulated$otu.table.sim.filter / lib.size
    rel.tab.pseudo <- matrix(
      NA_real_,
      nrow = nrow(rel.tab),
      ncol = ncol(rel.tab),
      dimnames = dimnames(rel.tab)
    )
    
    for (jj in seq_len(ncol(rel.tab))) {
      x <- rel.tab[, jj]
      positive <- x[x > 0 & is.finite(x)]
      
      if (length(positive) == 0) next
      rel.tab.pseudo[, jj] <- x + min(positive) / 2
    }
    
    pseudo.rowsum <- rowSums(rel.tab.pseudo)
    if (any(!is.finite(pseudo.rowsum)) || any(pseudo.rowsum <= 0)) {
      stop("Invalid pseudo-count relative abundance table.")
    }
    
    rel.tab.pseudo <- rel.tab.pseudo / pseudo.rowsum
    p.tmp <- rep(NA_real_, ncol(rel.tab.pseudo))
    
    for (jj in seq_len(ncol(rel.tab.pseudo))) {
      x <- rel.tab.pseudo[, jj]
      if (any(!is.finite(x)) || any(x <= 0)) next
      
      fit <- lm(log(x) ~ Y, data = comparison.meta)
      coef.tab <- summary(fit)$coefficients
      
      if ("Y" %in% rownames(coef.tab)) {
        p.tmp[jj] <- coef.tab["Y", "Pr(>|t|)"]
      }
    }
    
    names(p.tmp) <- feature.names
    p.tmp
  }, silent = TRUE)
  
  if (inherits(p.maaslin2, "try-error")) {
    p.maaslin2 <- make_na_p(feature.names)
    message("Seed ", seed.current, ": MaAsLin2-style analysis failed")
  } else {
    message("Seed ", seed.current, ": MaAsLin2-style analysis finished")
  }
  
  # ---------------------------------------------------------------------------
  # ANCOM-BC2
  # ---------------------------------------------------------------------------
  dat.assays <- S4Vectors::SimpleList(
    counts = t(simulated$otu.table.sim.filter)
  )
  
  smd <- S4Vectors::DataFrame(comparison.meta)
  
  tse <- TreeSummarizedExperiment::TreeSummarizedExperiment(
    assays = dat.assays,
    colData = smd
  )
  
  res.ancombc2 <- try(
    ANCOMBC::ancombc2(
      data = tse,
      assay_name = "counts",
      group = NULL,
      tax_level = NULL,
      fix_formula = ancom.formula,
      prv_cut = 0,
      verbose = FALSE
    ),
    silent = TRUE
  )
  
  p.ancombc2 <- extract_ancombc2_pvalue(
    res.ancombc2,
    feature.names
  )
  
  message(
    "Seed ", seed.current, ": ANCOM-BC2 ",
    if (inherits(res.ancombc2, "try-error")) "failed" else "finished"
  )
  
  # ---------------------------------------------------------------------------
  # LinDA
  # ---------------------------------------------------------------------------
  res.linda <- try(
    LinDA::linda(
      otu.tab = t(simulated$otu.table.sim.filter),
      meta = comparison.meta,
      formula = linda.formula,
      alpha = internal.linda.alpha,
      prev.cut = 0,
      lib.cut = 0,
      winsor.quan = NULL
    ),
    silent = TRUE
  )
  
  p.linda <- extract_linda_pvalue(
    res.linda,
    feature.names
  )
  
  message(
    "Seed ", seed.current, ": LinDA ",
    if (inherits(res.linda, "try-error")) "failed" else "finished"
  )
  
  pvalue.list <- c(
    pvalue.list,
    list(
      "LOCOM" = p.locom,
      "MaAsLin3" = p.maaslin3,
      "MaAsLin2" = p.maaslin2,
      "ANCOMBC2" = p.ancombc2,
      "LinDA" = p.linda
    )
  )
  
  # ---------------------------------------------------------------------------
  # Save raw taxon-level p-values for this shuffle
  # ---------------------------------------------------------------------------
  if (!identical(names(pvalue.list), c(
    "LOCOM-count-wr",
    "LOCOM-Wald-count-wr",
    "LOCOM-Wald-count-wr-yj",
    "LOCOM-Wald-count-wr-sig",
    "LOCOM-Wald-count-wr-pval",
    "LOCOM",
    "MaAsLin3",
    "MaAsLin2",
    "ANCOMBC2",
    "LinDA"
  ))) {
    stop("Unexpected method names/order in pvalue.list.")
  }
  
  bad.length <- vapply(
    pvalue.list,
    function(x) length(x) != length(feature.names),
    logical(1)
  )
  
  if (any(bad.length)) {
    stop(
      "Raw p-value length mismatch for method(s): ",
      paste(names(pvalue.list)[bad.length], collapse = ", ")
    )
  }
  
  raw.current <- data.frame(
    seed = rep(seed.current, length(feature.names)),
    taxon = feature.names,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  for (method.current in names(pvalue.list)) {
    p.current <- as.numeric(pvalue.list[[method.current]])
    
    # failed or unavailable tests stay NA ("not tested" is kept apart from "not rejected")
    p.current[!is.finite(p.current)] <- NA_real_
    
    invalid <- !is.na(p.current) & (p.current < 0 | p.current > 1)
    if (any(invalid)) {
      warning(
        "Seed ", seed.current, ", method ", method.current,
        ": ", sum(invalid), " p-value(s) outside [0,1]; saved as NA."
      )
      p.current[invalid] <- NA_real_
    }
    
    raw.current[[method.current]] <- p.current
  }
  
  # each shuffle is appended to the file immediately, so partial results survive a failure
  write.table(
    raw.current,
    file = raw.pvalue.file,
    sep = ",",
    quote = TRUE,
    row.names = FALSE,
    col.names = !file.initialized,
    append = file.initialized,
    na = "NA"
  )
  
  file.initialized <- TRUE
  n.completed <- n.completed + 1L
}

if (n.completed == 0L) {
  stop("No simulation replicate completed successfully.")
}

message(
  "Null shuffling completed for ", n.completed,
  " of ", n.rep, " requested shuffle(s)."
)
message("Raw taxon-level p-values saved to: ", raw.pvalue.file)
message(
  "Each row is one (seed, taxon) pair; empirical taxon-specific Type I error ",
  "can be calculated later as mean(p <= alpha, na.rm = TRUE) for each method."
)
