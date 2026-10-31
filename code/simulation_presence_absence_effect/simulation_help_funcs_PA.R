# Help functions for the presence-absence-effect simulations of the LOCOM2 paper

library(MIDASim)

simulate_midasim_PA = function(beta, seed, n.sam, causal.type, filter.thresh,
                               ref.otu = NULL, depth.mu = 10^4, depth.fix = NULL,
                               diff.lib = 0, template = "throat", beta.c = 0,
                               prop.case = 0.5, template.dir = ".") {

  set.seed(seed)
  n0.sam = as.integer(round(n.sam * (1-prop.case)))
  n1.sam = as.integer(round(n.sam * prop.case))
  C = NULL

  if (template == "throat") {
    load(file.path(template.dir, "fitted.throat.856.RData"), verbose = T)
    fitted = fitted.throat.856
    pi_taxa <- read.table(file.path(template.dir, "fit_dirmult_pi.txt"), header=FALSE, as.is=TRUE)[,1]
    ref.otu = "4414"

    if (beta.c == 0) {

      # taxon 468 is the reference taxon 4414 (the most abundant taxon): causal.type 2 skips it
      # (the 14th taxon with pi >= 0.01), causal.type 3 and 4 include it
      if (causal.type == 1){
        causal.otus <- c(which(pi_taxa >= 0.005))[1:20]
      }
      if (causal.type == 2){
        causal.otus <- c(which(pi_taxa >= 0.01))[c(1:13, 15:21)]
      }
      if (causal.type == 3) {
        causal.otus <- c(which(pi_taxa >= 0.02))[1:20]
      }
      if (causal.type == 4) {
        causal.otus <- c(c(which(pi_taxa >= 0.005))[1:9], 468)
      }

    } else {

      # Y affects the taxa in Y.index and C the taxa in C.index; the two sets partly overlap
      if (causal.type == 1){
        w <- which(pi_taxa >= 0.005)
        Y.index <- c(w[1:10], w[11:20])
        C.index <- c(w[21:30], w[11:20])
      }
      if (causal.type == 2){
        w <- which(pi_taxa >= 0.002)
        Y.index <- c(w[1:10], w[11:20])
        C.index <- c(w[21:30], w[11:20])
      }
      if (causal.type == 3){
        w <- which(pi_taxa >= 0.001)
        Y.index <- c(w[1:10], w[11:20])
        C.index <- c(w[21:30], w[11:20])
      }
      if (causal.type == 4){
        o <- order(pi_taxa, decreasing = TRUE)
        Y.index <- c(o[1:3], o[4:5])
        C.index <- c(o[31:33], o[4:5])
      }
      causal.otus <- Y.index
      # C is associated with Y: P(C = 1) = rho in controls and 1 - rho in cases
      rho <- 0.1
      C <- c(rbinom(n0.sam, 1, rho), rbinom(n1.sam, 1, 1 - rho))

    } # beta.c == 0?

  } # throat

  if (template == "momspi") {
    load(file.path(template.dir, "fitted.vaginal.RData"), verbose = T)
    fitted = fitted.vaginal
    pi_taxa = pi.taxa.vaginal
    ref.otu = ref.otu.vaginal

    if (beta.c == 0) {
      if (causal.type == 1){
        causal.otus <- c(which(pi.taxa.vaginal >= 0.0005 & pi.taxa.vaginal < 0.2))[1:20]
      }
      if (causal.type == 2){
        causal.otus <- c(which(pi.taxa.vaginal >= 0.001 & pi.taxa.vaginal < 0.3))[1:20]
      }
      if (causal.type == 3){
        causal.otus <- order(pi.taxa.vaginal, decreasing = TRUE)[2:6]
      }
      if (causal.type == 4){
        causal.otus <- order(pi.taxa.vaginal, decreasing = TRUE)[3:7]
      }
    } else {
      if (causal.type == 1){
        w = c(which(pi.taxa.vaginal >= 0.0005 & pi.taxa.vaginal < 0.2))
        Y.index <- c(w[1:10], w[11:20])
        C.index <- c(w[21:30], w[11:20])
      }
      if (causal.type == 2){
        causal.otus <- c(which(pi.taxa.vaginal >= 0.001 & pi.taxa.vaginal < 0.3))[1:20]
      }
      if (causal.type == 3){
        causal.otus <- order(pi.taxa.vaginal, decreasing = TRUE)[2:6]
      }
      if (causal.type == 4){
        causal.otus <- order(pi.taxa.vaginal, decreasing = TRUE)[3:7]
      }

      causal.otus <- Y.index
      rho <- 0.2
      C <- c(rbinom(n0.sam, 1, rho), rbinom(n1.sam, 1, 1 - rho))
    }

  } # momspi

  n.otus = fitted$n.taxa
  exp.beta = exp(beta)
  noncausal.otus <- setdiff(1:n.otus, causal.otus)
  Y <- c(rep(0, n0.sam), rep(1, n1.sam))

  # Library sizes (as in simulate_midasim()); MIDASim is always run at the library sizes depth.sim and
  # the counts are obtained by scaling the final relative abundances by depth.final
  #   diff.lib = 0: normal with mean depth.mu (or fixed at depth.fix)
  #   diff.lib > 0: common library size 2*10^4 for MIDASim; final library sizes with case:control
  #                 mean ratio 1.5:1 (diff.lib = 1), 3:1 (2), 5:1 (3) or 1:1 (4)
  if (diff.lib == 0) {

    if (is.null(depth.fix)) {
      depth.lower = 2000
      depth.sd = depth.mu/3

      depth.sim <- rnorm(n.sam, depth.mu, depth.sd)
      depth.sim[depth.sim < depth.lower] <- depth.lower
      depth.sim <- round(depth.sim)
    } else {
      depth.sim = rep(depth.fix, n.sam)
    }

    depth.final <- depth.sim

  } else {

    mean.lib.size.control = 2*10^4

    if (diff.lib == 1) { # 1.5:1
      mean.lib.size.case = 1.5 * mean.lib.size.control
    }
    if (diff.lib == 2) { # 3:1
      mean.lib.size.case = 3 * mean.lib.size.control
    }
    if (diff.lib == 3) { # 5:1
      mean.lib.size.case = 5 * mean.lib.size.control
    }
    if (diff.lib == 4) { # 1:1
      mean.lib.size.case = 1 * mean.lib.size.control
    }

    depth.diff.lib <- c(
      round(rnorm(n0.sam,
                  mean.lib.size.control,
                  mean.lib.size.control/3)),
      round(rnorm(n1.sam,
                  mean.lib.size.case,
                  mean.lib.size.case/3))
    )

    depth.diff.lib <- pmax(depth.diff.lib, 2000)

    depth.sim <- rep(mean.lib.size.control, n.sam)
    depth.final <- depth.diff.lib
  }

  otu.table.sim <- matrix(nrow = n.sam, ncol = n.otus)

  # relative abundances from MIDASim: no Y effect; C acts multiplicatively on the taxa in C.index
  if (beta.c == 0) {

    fitted.modified <- MIDASim.modify(
      fitted,
      lib.size = depth.sim,
      mean.rel.abund = fitted$mean.rel.abund
    )

    temp <- MIDASim(fitted.modified)

    otu.ra.sim <- temp$sim_rel

  } else {

    individual.rel.abund <- matrix(
      fitted$mean.rel.abund,
      nrow = n.sam,
      ncol = length(fitted$mean.rel.abund),
      byrow = TRUE
    )

    individual.rel.abund[, C.index] <-
      individual.rel.abund[, C.index, drop = FALSE] *
      exp(
        outer(
          C,
          rep(beta.c, length(C.index))
        )
      )

    individual.rel.abund <-
      individual.rel.abund /
      rowSums(individual.rel.abund)

    fitted.modified <- MIDASim.modify(
      fitted,
      lib.size = depth.sim,
      individual.rel.abund = individual.rel.abund
    )

    temp <- MIDASim(fitted.modified)

    otu.ra.sim <- temp$sim_rel
  }

  # Y effect: dropout of the causal taxa in the case samples with probability beta
  otu.ra.sim <- add_dropout(
    ra.mat = otu.ra.sim,
    sample.index = which(Y == 1),
    taxa.index = causal.otus,
    dropout.prob = beta
  )

  otu.table.sim <- round(
    otu.ra.sim * depth.final
  )

  colnames(otu.table.sim) <- colnames(fitted$mat01)

  # keep taxa present in at least filter.thresh * n.sam samples
  prop.count <- colSums(otu.table.sim > 0)
  otus.keep <- which(prop.count >= filter.thresh * n.sam)

  otu.table.sim.filter <- otu.table.sim[, otus.keep]
  if (is.null(ref.otu)) {
    ref.otu.filter = NULL
  } else {
    ref.otu.filter <- which(colnames(otu.table.sim.filter) == ref.otu)
  }
  causal.otus.filter = which(otus.keep %in% causal.otus)
  noncausal.otus.filter = setdiff(1:length(otus.keep), causal.otus.filter)


  simulated = list(causal.otus.filter = causal.otus.filter,
                   noncausal.otus.filter = noncausal.otus.filter,
                   ref.otu.filter = ref.otu.filter,
                   otu.table.sim.filter = otu.table.sim.filter,
                   Y = Y,
                   C = C,
                   exp.beta = exp.beta,
                   beta = beta)

  return(simulated)

}


# IBD (female) template with the confounder age; see simulate_ibd_data() in Simulation_Main/simulation_help_funcs.R
# for the contents of fitted.ibd.female.RData and the library-size designs. Here beta and beta.c are
# dropout probabilities; the Y effect is dropout of the causal taxa (id.causal for causal.type 1,
# id.causal.20 restricted by abundance for causal.type 2-4) in the case samples.

simulate_ibd_data_PA <- function(seed, n.sam, prop.case, diff.lib,
                                 causal.type, beta, beta.c, filter.thresh, template.dir = ".") {

  if (beta < 0 || beta > 1) {
    stop("beta must be between 0 and 1.")
  }

  if (beta.c < 0 || beta.c > 1) {
    stop("beta.c must be between 0 and 1.")
  }

  # fitted.ibd.female.RData also contains an object named beta.c (the age effects used by
  # simulate_ibd_data()); keep the dropout probability given as argument
  beta.c.dropout <- beta.c

  load(file.path(template.dir, "fitted.ibd.female.RData"), verbose = FALSE)
  fitted <- fit.mix

  set.seed(seed)

  pi_taxa <- fitted$mean.rel.abund

  if (causal.type == 1) {

    causal.otus <- id.causal

  } else {

    if (causal.type == 2) {

      target_idx <- seq_along(id.causal.20)

    } else if (causal.type == 3) {

      target_idx <- which(
        id.causal.20 %in% which(pi_taxa >= 0.001)
      )

    } else if (causal.type == 4) {

      target_idx <- which(
        id.causal.20 %in% which(pi_taxa >= 0.005)
      )
    }

    causal.otus <- id.causal.20[target_idx]
  }

  n1.sam <- as.integer(round(n.sam * prop.case))
  n0.sam <- n.sam - n1.sam

  Y <- c(
    rep(0, n0.sam),
    rep(1, n1.sam)
  )

  # Library sizes
  #   diff.lib = 0: resampled from the observed library sizes of the template
  #   diff.lib > 0: normal around mean(lib.size.case) times a control:case multiplier of
  #                 1:1 (diff.lib = 1, 5), 1.5:1 (2), 3:1 (3), 5:1 (4), 1:3 (6) or 1:5 (7)
  if (diff.lib == 0) {

    depth.diff.lib <- c(
      sample(lib.size.control, n0.sam, replace = TRUE),
      sample(lib.size.case, n1.sam, replace = TRUE)
    )

    depth.sim <- rep(
      round(mean(lib.size.control)),
      n.sam
    )

  } else {

    base_mean <- mean(lib.size.case)

    mult_ctrl <- 1
    mult_case <- 1

    if (diff.lib == 2) mult_ctrl <- 1.5
    if (diff.lib == 3) mult_ctrl <- 3.0
    if (diff.lib == 4) mult_ctrl <- 5.0

    if (diff.lib == 6) mult_case <- 3.0
    if (diff.lib == 7) mult_case <- 5.0

    mean.ctrl <- base_mean * mult_ctrl
    mean.case <- base_mean * mult_case

    depth.diff.lib <- c(
      round(
        rnorm(
          n0.sam,
          mean.ctrl,
          mean.ctrl / 2
        )
      ),
      round(
        rnorm(
          n1.sam,
          mean.case,
          mean.case / 2
        )
      )
    )

    depth.diff.lib <- pmax(
      depth.diff.lib,
      2000
    )

    depth.sim <- rep(
      round(mean(lib.size.control)),
      n.sam
    )

    if (diff.lib == 1) {
      depth.sim <- rep(
        round(base_mean),
        n.sam
      )
    }
  }

  # confounder C (age), resampled from the observed ages of controls and cases
  C.control <- sample(
    age.control,
    n0.sam,
    replace = TRUE
  )

  C.case <- sample(
    age.case,
    n1.sam,
    replace = TRUE
  )

  C <- c(
    C.control,
    C.case
  )

  # relative abundances from MIDASim under the null (no Y and no C effect)
  fitted.modified <- MIDASim.modify(
    fitted,
    lib.size = depth.sim,
    mean.rel.abund = fitted$mean.rel.abund
  )

  temp <- MIDASim(fitted.modified)

  otu.ra.sim <- temp$sim_rel

  # Y effect: dropout of the causal taxa in the case samples with probability beta
  otu.ra.sim <- add_dropout(
    ra.mat = otu.ra.sim,
    sample.index = which(Y == 1),
    taxa.index = causal.otus,
    dropout.prob = beta
  )

  # C effect: dropout of the taxa in id.C, with probability beta.c, in the samples with age above the median
  if (beta.c.dropout > 0) {

    C.binary <- as.integer(
      C > median(C)
    )

    otu.ra.sim <- add_dropout(
      ra.mat = otu.ra.sim,
      sample.index = which(C.binary == 1),
      taxa.index = id.C,
      dropout.prob = beta.c.dropout
    )

  } else {

    C.binary <- NULL
  }

  otu.table.sim <- round(
    otu.ra.sim * depth.diff.lib
  )

  colnames(otu.table.sim) <- colnames(fitted$mat01)

  # keep taxa present in at least filter.thresh * n.sam samples
  prop.count <- colSums(
    otu.table.sim > 0
  )

  otus.keep <- which(
    prop.count >= filter.thresh * n.sam
  )

  otu.table.sim.filter <- otu.table.sim[
    ,
    otus.keep,
    drop = FALSE
  ]

  causal.otus.filter <- which(
    otus.keep %in% causal.otus
  )

  noncausal.otus.filter <- setdiff(
    seq_along(otus.keep),
    causal.otus.filter
  )

  return(
    list(
      causal.otus.filter = causal.otus.filter,
      noncausal.otus.filter = noncausal.otus.filter,
      ref.otu.filter = NULL,
      otu.table.sim.filter = otu.table.sim.filter,
      Y = Y,
      C = matrix(C),
      C.binary = C.binary,
      beta = beta
    )
  )
}


# Dropout: set the taxa in taxa.index to zero in the samples in sample.index, independently with
# probability dropout.prob, and renormalize those samples

add_dropout <- function(ra.mat,
                        sample.index,
                        taxa.index,
                        dropout.prob) {

  if (dropout.prob == 0 || length(sample.index) == 0 || length(taxa.index) == 0) {
    return(ra.mat)
  }

  drop.mat <- matrix(
    rbinom(length(sample.index) * length(taxa.index), size = 1, prob = dropout.prob),
    nrow = length(sample.index),
    ncol = length(taxa.index)
  )

  ra.sub <- ra.mat[sample.index, taxa.index, drop = FALSE]
  ra.sub[drop.mat == 1] <- 0
  ra.mat[sample.index, taxa.index] <- ra.sub

  row.total <- rowSums(ra.mat[sample.index, , drop = FALSE])
  ra.mat[sample.index,] <- ra.mat[sample.index,,drop = FALSE] / row.total

  return(ra.mat)
}
