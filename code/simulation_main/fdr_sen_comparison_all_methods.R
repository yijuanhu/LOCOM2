# Main simulation of the LOCOM2 paper: FDR and sensitivity of LOCOM2 and of the competing methods
# (LOCOM, MaAsLin3, MaAsLin2-style log-linear model, ANCOM-BC2, LinDA) on MIDASim data.

args=(commandArgs(TRUE))

if(length(args) < 12) stop("Not enough arguments provided.")

n.sam <- as.numeric(args[1])
causal.type <- as.numeric(args[2])
exp.beta <- as.numeric(args[3])
start.seed <- as.numeric(args[4])
n.rep <- as.numeric(args[5])
diff.lib <- as.numeric(args[6])
beta.c <- as.numeric(args[7])
sim.method <- as.numeric(args[8])
n.perm.max <- as.numeric(args[9])
prop.case <- as.numeric(args[10])
Firth.thresh <- as.numeric(args[11])
n.wald <- as.numeric(args[12])

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

# folder with the MIDASim templates and the LOCOM2 / LOCOM code
data.dir <- "../data_templates"

source(file.path(data.dir, "LOCOM2_func.R"))
sourceCpp(file.path(data.dir, "Newton_LOCOM2.cpp"))
source(file.path(data.dir, "LOCOM_fun.R"))
sourceCpp(file.path(data.dir, "Newton.cpp"))
source("simulation_help_funcs.R")

sim.template <- switch(as.character(sim.method),
                       "1" = "MIDASim_throat",
                       "2" = "MIDASim_momspi",
                       "3" = "MIDASim_ibd_female_recover",
                       stop("sim.method must be 1 (throat), 2 (momspi) or 3 (IBD)."))
beta = log(exp.beta)
filter.thresh = min(0.1, 10 / n.sam)

fdr.target <- 0.2

# ============================================================
# Helper functions
# ============================================================

# realized FDR, sensitivity and number of rejections of one method in one replicate
calc_fdr_sensitivity <- function(qvalue, causal.otus, noncausal.otus,
                                 fdr.target = 0.2) {
  qvalue <- as.numeric(qvalue)

  if (length(qvalue) == 0 || all(!is.finite(qvalue))) {
    return(c(fdr = NA_real_, sen = NA_real_, n.rej = NA_real_))
  }

  # taxa without a q-value (not tested by the method) are treated as non-significant
  qvalue[!is.finite(qvalue)] <- 1
  qvalue <- pmin(pmax(qvalue, 0), 1)

  rejected <- which(qvalue <= fdr.target)
  n.rej <- length(rejected)

  fdr <- if (n.rej == 0) {
    0
  } else {
    sum(rejected %in% noncausal.otus) / n.rej
  }

  # sensitivity is undefined under the global null (no causal taxa)
  sen <- if (length(causal.otus) == 0) {
    NA_real_
  } else {
    sum(rejected %in% causal.otus) / length(causal.otus)
  }

  c(fdr = fdr, sen = sen, n.rej = n.rej)
}

# place the values reported for taxa in the order of target_taxa (NA for taxa not reported)
align_values <- function(value, taxa, target_taxa) {
  out <- rep(NA_real_, length(target_taxa))
  names(out) <- target_taxa

  idx <- match(target_taxa, as.character(taxa))
  keep <- !is.na(idx)
  out[keep] <- as.numeric(value[idx[keep]])

  out
}

extract_locom_qvalue <- function(res, target_taxa) {
  if (!is.null(res$q.otu)) {
    q <- res$q.otu
  } else if (!is.null(res$q.otu.wr)) {
    q <- res$q.otu.wr
  } else {
    return(rep(NA_real_, length(target_taxa)))
  }

  q <- as.numeric(q)
  if (length(q) != length(target_taxa)) {
    stop("LOCOM q-value vector length does not match the filtered OTU table.")
  }
  names(q) <- target_taxa
  q
}

extract_ancombc2_qvalue <- function(res, target_taxa) {

  if (inherits(res, "try-error") || is.null(res$res)) {
    return(rep(NA_real_, length(target_taxa)))
  }

  tab <- res$res
  taxa <- if ("taxon" %in% colnames(tab)) {
    as.character(tab$taxon)
  } else {
    rownames(tab)
  }

  # adjusted p-values of ANCOM-BC2; taxa that fail the sensitivity-score criterion are not declared significant
  if ("q_Y" %in% colnames(tab)) {

    q <- tab$q_Y

    if ("passed_ss_Y" %in% colnames(tab)) {
      q[!tab$passed_ss_Y | is.na(tab$passed_ss_Y)] <- 1
    }

    return(
      align_values(
        q,
        taxa,
        target_taxa
      )
    )
  }

  # fallback: BH-adjust the raw p-values if q_Y is unavailable
  if ("p_Y" %in% colnames(tab)) {

    p <- tab$p_Y

    if ("passed_ss_Y" %in% colnames(tab)) {
      p[!tab$passed_ss_Y | is.na(tab$passed_ss_Y)] <- 1
    }

    q <- p.adjust(p, method = "BH")

    return(
      align_values(
        q,
        taxa,
        target_taxa
      )
    )
  }

  rep(NA_real_, length(target_taxa))
}

extract_linda_qvalue <- function(res, target_taxa) {
  if (inherits(res, "try-error") || is.null(res$output) || length(res$output) == 0) {
    return(rep(NA_real_, length(target_taxa)))
  }

  if ("Y" %in% names(res$output)) {
    tab <- res$output[["Y"]]
  } else {
    y_idx <- grep("Y", names(res$output))
    tab <- if (length(y_idx) > 0) res$output[[y_idx[1]]] else res$output[[1]]
  }

  if (is.null(tab)) {
    return(rep(NA_real_, length(target_taxa)))
  }

  taxa <- rownames(tab)

  # the adjusted-p-value column is named differently across LinDA versions
  q_candidates <- c("padj", "qvalue", "qval", "adj.pvalue", "adj_pvalue")
  q_col <- q_candidates[q_candidates %in% colnames(tab)]

  if (length(q_col) > 0) {
    return(align_values(tab[[q_col[1]]], taxa, target_taxa))
  }

  if (!"pvalue" %in% colnames(tab)) {
    return(rep(NA_real_, length(target_taxa)))
  }

  p <- align_values(tab$pvalue, taxa, target_taxa)
  p.adjust(p, method = "BH")
}

# ============================================================
# Simulation replicates
# ============================================================

success <- start.seed
tmp.performance <- NULL
fail <- 0

while (success < (n.rep + start.seed)) {

  # ---------- simulate data ----------

  if (sim.method == 1) {
    simulated <- try(
      simulate_midasim(
        beta = beta,
        seed = success,
        n.sam = n.sam,
        causal.type = causal.type,
        filter.thresh = filter.thresh,
        diff.lib = diff.lib,
        template = "throat",
        beta.c = beta.c,
        prop.case = prop.case,
        template.dir = data.dir
      ),
      silent = TRUE
    )
  } else if (sim.method == 2) {
    simulated <- try(
      simulate_midasim(
        beta = beta,
        seed = success,
        n.sam = n.sam,
        causal.type = causal.type,
        filter.thresh = filter.thresh,
        diff.lib = diff.lib,
        template = "momspi",
        beta.c = beta.c,
        prop.case = prop.case,
        template.dir = data.dir
      ),
      silent = TRUE
    )
  } else {
    simulated <- try(
      simulate_ibd_data(
        seed = success,
        n.sam = n.sam,
        prop.case = prop.case,
        diff.lib = diff.lib,
        causal.type = causal.type,
        exp.beta = exp.beta,
        beta.c = beta.c,
        filter.thresh = filter.thresh,
        template.dir = data.dir
      ),
      silent = TRUE
    )
  }

  if (inherits(simulated, "try-error")) {
    message("Seed ", success, ": simulation failed")
    success <- success + 1
    fail <- fail + 1
    next
  }

  feature.names <- colnames(simulated$otu.table.sim.filter)
  if (is.null(feature.names)) {
    feature.names <- paste0("OTU", seq_len(ncol(simulated$otu.table.sim.filter)))
    colnames(simulated$otu.table.sim.filter) <- feature.names
  }

  # under the global null all retained taxa are non-causal and sensitivity is undefined
  if (beta == 0) {
    simulated$causal.otus.filter <- integer(0)
    simulated$noncausal.otus.filter <- seq_len(ncol(simulated$otu.table.sim.filter))
  }

  # ---------- LOCOM2 ----------

  res.locom2 <- locom2(
    otu.table = simulated$otu.table.sim.filter,
    Y = simulated$Y,
    C = simulated$C,
    fdr.nominal = fdr.target,
    ref.otu = simulated$ref.otu.filter,
    n.perm.max = n.perm.max,
    n.rej.stop = 100,
    n_cores = 1,
    filter.thresh = 0,
    seed = success,
    n.wald = n.wald,
    Firth.thresh = Firth.thresh
  )

  qvalue.list <- list(
    "LOCOM-count-wr" = res.locom2$q.otu.wr,
    "LOCOM-Wald-count-wr" = res.locom2$q.otu.wald.wr,
    "LOCOM-Wald-count-wr-yj" = res.locom2$q.otu.wald.wr.yj,
    "LOCOM-Wald-count-wr-sig" = res.locom2$q.otu.wald.wr.sig,
    "LOCOM-Wald-count-wr-pval" = res.locom2$q.otu.wald.wr.pval
  )

  # ---------- comparison methods ----------

  sample.names <- paste0("sam", seq_len(n.sam))
  rownames(simulated$otu.table.sim.filter) <- sample.names

  # the comparison methods are adjusted for C whenever a covariate was simulated
  adjust.for.C <- (length(simulated$C) > 0)

  if (adjust.for.C) {
    comparison.meta <- data.frame(
      Y = simulated$Y,
      C = as.numeric(simulated$C),
      row.names = sample.names
    )
    maaslin.formula <- "~ Y + C"
    ancom.formula <- "Y + C"
    linda.formula <- "~ Y + C"
  } else {
    comparison.meta <- data.frame(
      Y = simulated$Y,
      row.names = sample.names
    )
    maaslin.formula <- "~ Y"
    ancom.formula <- "Y"
    linda.formula <- "~ Y"
  }

  # ---------- LOCOM ----------

  res.locom <- try(
    locom(
      otu.table = simulated$otu.table.sim.filter,
      Y = simulated$Y,
      C = simulated$C,
      fdr.nominal = fdr.target,
      n.cores = 1,
      n.perm.max = n.perm.max,
      filter.thresh = 0,
      seed = success
    ),
    silent = TRUE
  )

  if (inherits(res.locom, "try-error")) {
    q.locom <- rep(NA_real_, length(feature.names))
    names(q.locom) <- feature.names
    message("Seed ", success, ": LOCOM failed")
  } else {
    q.locom <- extract_locom_qvalue(res.locom, feature.names)
    message("Seed ", success, ": LOCOM finished")
  }

  # ---------- MaAsLin3 ----------

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
        tmp_res <- suppressMessages(
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

    tmp_res
  }, silent = TRUE)

  if (inherits(res.maaslin3, "try-error")) {
    q.maaslin3 <- rep(NA_real_, length(feature.names))
    names(q.maaslin3) <- feature.names
    message("Seed ", success, ": MaAsLin3 failed")
  } else {
    maaslin3.tab <- res.maaslin3$fit_data_abundance$results
    maaslin3.tab <- maaslin3.tab[
      maaslin3.tab$metadata == "Y",
      ,
      drop = FALSE
    ]

    q.maaslin3 <- align_values(
      value = maaslin3.tab$qval_joint,
      taxa = maaslin3.tab$feature,
      target_taxa = feature.names
    )
    message("Seed ", success, ": MaAsLin3 finished")
  }

  # ---------- MaAsLin2-style log-linear model ----------
  # relative abundances (TSS) + half of the smallest non-zero value of each taxon as pseudo-count,
  # log-transformed and analyzed with a linear model; BH-adjusted p-values for Y

  q.maaslin2 <- try({
    rel.tab <- simulated$otu.table.sim.filter /
      rowSums(simulated$otu.table.sim.filter)

    rel.tab.pseudo <- apply(
      rel.tab,
      2,
      function(x) {
        positive <- x[x > 0]
        if (length(positive) == 0) {
          return(rep(NA_real_, length(x)))
        }
        x + min(positive) / 2
      }
    )

    rel.tab.pseudo <- rel.tab.pseudo / rowSums(rel.tab.pseudo)

    p.tmp <- rep(NA_real_, ncol(rel.tab.pseudo))

    for (jj in seq_len(ncol(rel.tab.pseudo))) {
      if (all(!is.finite(rel.tab.pseudo[, jj]))) next

      if (adjust.for.C) {
        fit <- lm(
          log(rel.tab.pseudo[, jj]) ~ Y + C,
          data = comparison.meta
        )
      } else {
        fit <- lm(
          log(rel.tab.pseudo[, jj]) ~ Y,
          data = comparison.meta
        )
      }

      coef.tab <- summary(fit)$coefficients
      if ("Y" %in% rownames(coef.tab)) {
        p.tmp[jj] <- coef.tab["Y", "Pr(>|t|)"]
      }
    }

    q.tmp <- p.adjust(p.tmp, method = "BH")
    names(q.tmp) <- feature.names
    q.tmp
  }, silent = TRUE)

  if (inherits(q.maaslin2, "try-error")) {
    q.maaslin2 <- rep(NA_real_, length(feature.names))
    names(q.maaslin2) <- feature.names
    message("Seed ", success, ": MaAsLin2-style analysis failed")
  } else {
    message("Seed ", success, ": MaAsLin2-style analysis finished")
  }

  # ---------- ANCOM-BC2 ----------

  ancom.meta <- comparison.meta
  ancom.meta$Sample <- sample.names

  dat.assays <- S4Vectors::SimpleList(
    counts = t(simulated$otu.table.sim.filter)
  )
  smd <- S4Vectors::DataFrame(ancom.meta)
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

  q.ancombc2 <- extract_ancombc2_qvalue(
    res.ancombc2,
    feature.names
  )

  if (inherits(res.ancombc2, "try-error")) {
    message("Seed ", success, ": ANCOM-BC2 failed")
  } else {
    message("Seed ", success, ": ANCOM-BC2 finished")
  }

  # ---------- LinDA ----------

  linda.meta <- comparison.meta

  res.linda <- try(
    LinDA::linda(
      otu.tab = t(simulated$otu.table.sim.filter),
      meta = linda.meta,
      formula = linda.formula,
      alpha = fdr.target,
      prev.cut = 0,
      lib.cut = 0,
      winsor.quan = NULL
    ),
    silent = TRUE
  )

  q.linda <- extract_linda_qvalue(
    res.linda,
    feature.names
  )

  if (inherits(res.linda, "try-error")) {
    message("Seed ", success, ": LinDA failed")
  } else {
    message("Seed ", success, ": LinDA finished")
  }

  qvalue.list <- c(
    qvalue.list,
    list(
      "LOCOM" = q.locom,
      "MaAsLin3" = q.maaslin3,
      "MaAsLin2" = q.maaslin2,
      "ANCOMBC2" = q.ancombc2,
      "LinDA" = q.linda
    )
  )

  # ---------- realized FDR and sensitivity of every method ----------

  bad.length <- vapply(
    qvalue.list,
    function(x) length(x) != length(feature.names),
    logical(1)
  )
  if (any(bad.length)) {
    stop(
      "Q-value length mismatch for method(s): ",
      paste(names(qvalue.list)[bad.length], collapse = ", ")
    )
  }

  performance.current <- do.call(
    rbind,
    lapply(names(qvalue.list), function(method.name.current) {
      stat <- calc_fdr_sensitivity(
        qvalue = qvalue.list[[method.name.current]],
        causal.otus = simulated$causal.otus.filter,
        noncausal.otus = simulated$noncausal.otus.filter,
        fdr.target = fdr.target
      )

      data.frame(
        method = method.name.current,
        fdr = unname(stat["fdr"]),
        sen = unname(stat["sen"]),
        n.rej = unname(stat["n.rej"]),
        stringsAsFactors = FALSE
      )
    })
  )

  performance.current$seed <- success
  performance.current$fdr.target <- fdr.target
  performance.current$beta <- beta
  performance.current$exp.beta <- exp.beta
  performance.current$n.sam <- n.sam
  performance.current$causal.type <- causal.type
  performance.current$n.perm.max <- n.perm.max
  performance.current$sim.template <- sim.template
  performance.current$diff.lib <- diff.lib
  performance.current$beta.c <- mean(beta.c)
  performance.current$sim.method <- sim.method
  performance.current$prop.case <- prop.case
  performance.current$Firth.thresh <- Firth.thresh
  performance.current$n.wald <- n.wald

  tmp.performance <- rbind(tmp.performance, performance.current)

  success <- success + 1
}

# ============================================================
# Save results
# ============================================================

n.success <- n.rep - fail

file.name <- paste0(
  "fdr_all_methods_target", fdr.target,
  "_thresh", filter.thresh,
  "_sim_", start.seed,
  "_n.rep_", n.rep,
  "_cau", causal.type,
  "_n", n.sam,
  "_beta", beta,
  "_difflib", diff.lib,
  "_simmethod", sim.method,
  "_betac", mean(beta.c),
  "_pcase", prop.case,
  "_firth", Firth.thresh,
  "_nwald", n.wald,
  "_otu.txt"
)

write.table(
  tmp.performance,
  file = file.name,
  quote = FALSE,
  sep = ",",
  col.names = TRUE,
  row.names = FALSE
)

message(
  "Simulation completed with ", n.success,
  " successful replicate(s). Results saved to: ",
  file.name
)
