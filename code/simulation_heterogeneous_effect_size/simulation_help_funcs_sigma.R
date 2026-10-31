# Help function for the heterogeneous-effect-size simulations of the LOCOM2 paper

library(MIDASim)

simulate_midasim_sigma = function(
    sigma = 0,
    effect.scenario = c("balanced", "biased"),
    seed, n.sam, causal.type, filter.thresh,
    ref.otu = NULL, depth.mu = 10^4, depth.fix = NULL,
    diff.lib = 0, template = "throat", beta.c = 0,
    prop.case = 0.5, template.dir = ".") {

  effect.scenario <- match.arg(effect.scenario)

  set.seed(seed)

  n0.sam = as.integer(round(n.sam * (1-prop.case)))
  n1.sam = as.integer(round(n.sam * prop.case))
  C = NULL
  load(file.path(template.dir, "u.RData"), verbose = T)

  if (template == "throat") {

    load(file.path(template.dir, "fitted.throat.856.RData"), verbose = T)
    fitted = fitted.throat.856
    pi_taxa <- read.table(
      file.path(template.dir, "fit_dirmult_pi.txt"),
      header = FALSE,
      as.is = TRUE
    )[,1]

    ref.otu = "4414"

    if (beta.c == 0) {

      # taxon 468 is the reference taxon 4414 (the most abundant taxon), so the reference
      # taxon is itself causal under causal.type 2-4
      if (causal.type == 1){
        causal.otus <- c(which(pi_taxa >= 0.005))[1:20]
      }

      if (causal.type == 2){
        causal.otus <- order(pi_taxa, decreasing = TRUE)[1:5]
      }

      if (causal.type == 3) {
        causal.otus <- c(
          c(which(pi_taxa >= 0.001))[1:20],
          468
        )
      }

      if (causal.type == 4) {
        causal.otus <- c(
          c(which(pi_taxa >= 0.005))[1:9],
          468
        )
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

      rho <- 0.5
      C <- c(
        rbinom(n0.sam, 1, rho),
        rbinom(n1.sam, 1, 1-rho)
      )
    }
  }


  if (template == "momspi") {

    load(file.path(template.dir, "fitted.vaginal.RData"), verbose = T)

    fitted = fitted.vaginal
    pi_taxa = pi.taxa.vaginal
    ref.otu = ref.otu.vaginal

    if (beta.c == 0) {

      if (causal.type == 1){
        causal.otus <- c(
          which(
            pi.taxa.vaginal >= 0.0005 &
              pi.taxa.vaginal < 0.2
          )
        )[1:20]
      }

      if (causal.type == 2){
        causal.otus <- c(
          which(
            pi.taxa.vaginal >= 0.001 &
              pi.taxa.vaginal < 0.3
          )
        )[1:20]
      }

      if (causal.type == 3){
        causal.otus <-
          order(pi.taxa.vaginal, decreasing = TRUE)[2:6]
      }

      if (causal.type == 4){
        causal.otus <-
          order(pi.taxa.vaginal, decreasing = TRUE)[3:7]
      }

    } else {

      if (causal.type == 1){

        w <- which(
          pi.taxa.vaginal >= 0.0005 &
            pi.taxa.vaginal < 0.2
        )

        Y.index <- c(w[1:10], w[11:20])
        C.index <- c(w[21:30], w[11:20])
      }

      causal.otus <- Y.index

      rho <- 0.5

      C <- c(
        rbinom(n0.sam, 1, rho),
        rbinom(n1.sam, 1, 1-rho)
      )
    }
  }

  n.otus = fitted$n.taxa

  # taxon-specific effects of the causal taxa
  beta0.taxa <- u

  if (effect.scenario == "balanced") {

    beta0.taxa[neg.weak] <- -beta0.taxa[neg.weak]

  } else if (effect.scenario == "biased") {

    beta0.taxa[neg.strong] <- -beta0.taxa[neg.strong]
  }

  beta.taxa <- sigma * beta0.taxa
  exp.beta.taxa <- exp(beta.taxa)

  names(beta0.taxa) <- causal.otus
  names(beta.taxa) <- causal.otus
  names(exp.beta.taxa) <- causal.otus


  noncausal.otus <- setdiff(
    1:n.otus,
    causal.otus
  )

  Y <- c(
    rep(0, n0.sam),
    rep(1, n1.sam)
  )


  # Library sizes (as in simulate_midasim())
  #   diff.lib = 0: normal with mean depth.mu (or fixed at depth.fix); counts are taken from MIDASim
  #   diff.lib > 0: MIDASim is run at a common library size of 2*10^4 and the simulated relative
  #                 abundances are scaled by sample-specific library sizes whose case:control mean
  #                 ratio is 1.5:1 (diff.lib = 1), 3:1 (2), 5:1 (3) or 1:1 (4)

  if (diff.lib == 0) {

    if (is.null(depth.fix)) {

      depth.lower = 2000
      depth.sd = depth.mu/3

      depth.sim <- rnorm(
        n.sam,
        depth.mu,
        depth.sd
      )

      depth.sim[depth.sim < depth.lower] <- depth.lower
      depth.sim <- round(depth.sim)

    } else {

      depth.sim = rep(depth.fix, n.sam)
    }

    depth.diff.lib <- depth.sim

  } else {

    mean.lib.size.control = 2*10^4

    if (diff.lib == 1) {
      mean.lib.size.case = 1.5 * mean.lib.size.control
    }

    if (diff.lib == 2) {
      mean.lib.size.case = 3 * mean.lib.size.control
    }

    if (diff.lib == 3) {
      mean.lib.size.case = 5 * mean.lib.size.control
    }

    if (diff.lib == 4) {
      mean.lib.size.case = mean.lib.size.control
    }

    depth.diff.lib <- c(
      round(
        rnorm(
          n0.sam,
          mean.lib.size.control,
          mean.lib.size.control/3
        )
      ),
      round(
        rnorm(
          n1.sam,
          mean.lib.size.case,
          mean.lib.size.case/3
        )
      )
    )

    depth.sim <- rep(
      mean.lib.size.control,
      n.sam
    )

    depth.diff.lib <- pmax(
      depth.diff.lib,
      2000
    )
  }


  otu.table.sim <- matrix(
    nrow = n.sam,
    ncol = n.otus
  )


  if (beta.c == 0) {

    # controls: template mean relative abundances
    fitted.modified <- MIDASim.modify(
      fitted,
      lib.size = depth.sim[Y == 0],
      mean.rel.abund = fitted$mean.rel.abund
    )

    temp <- MIDASim(fitted.modified)

    if (diff.lib == 0) {

      otu.table.sim[Y == 0, ] <- temp$sim_count

    } else {

      otu.table.sim[Y == 0, ] <-
        round(
          temp$sim_rel *
            depth.diff.lib[Y == 0]
        )
    }


    # cases: mean relative abundances of the causal taxa multiplied by exp(beta.taxa)
    new.mean.rel.abund <- fitted$mean.rel.abund

    new.mean.rel.abund[causal.otus] <-
      exp.beta.taxa *
      new.mean.rel.abund[causal.otus]

    new.mean.rel.abund <-
      new.mean.rel.abund /
      sum(new.mean.rel.abund)


    fitted.modified <- MIDASim.modify(
      fitted,
      lib.size = depth.sim[Y == 1],
      mean.rel.abund = new.mean.rel.abund
    )

    temp <- MIDASim(fitted.modified)

    if (diff.lib == 0) {

      otu.table.sim[Y == 1, ] <- temp$sim_count

    } else {

      otu.table.sim[Y == 1, ] <-
        round(
          temp$sim_rel *
            depth.diff.lib[Y == 1]
        )
    }


  } else {

    # four groups defined by (Y, C); samples are ordered so that each group is a block of n.sam/4
    for (group in 1:4) {

      Y.group <-
        Y[(group-1)*n.sam/4 + 1]

      C.group <-
        C[(group-1)*n.sam/4 + 1]

      effect1 <- Y.group * beta.taxa

      effect2 <- C.group * beta.c

      group.ind <-
        ((group-1)*n.sam/4 + 1):
        (group*n.sam/4)

      lib.size.group <-
        depth.sim[group.ind]


      new.mean.rel.abund <-
        fitted$mean.rel.abund

      new.mean.rel.abund[Y.index] <-
        exp(effect1) *
        new.mean.rel.abund[Y.index]

      new.mean.rel.abund[C.index] <-
        exp(effect2) *
        new.mean.rel.abund[C.index]

      new.mean.rel.abund <-
        new.mean.rel.abund /
        sum(new.mean.rel.abund)


      fitted.modified <- MIDASim.modify(
        fitted,
        lib.size = lib.size.group,
        mean.rel.abund = new.mean.rel.abund
      )

      temp <- MIDASim(fitted.modified)


      if (diff.lib == 0) {

        otu.table.sim[group.ind, ] <-
          temp$sim_count

      } else {

        otu.table.sim[group.ind, ] <-
          round(
            temp$sim_rel *
              depth.diff.lib[group.ind]
          )
      }
    }
  }


  colnames(otu.table.sim) <-
    colnames(fitted$mat01)


  # keep taxa present in at least filter.thresh * n.sam samples
  prop.count <-
    colSums(otu.table.sim > 0)

  otus.keep <-
    which(
      prop.count >=
        filter.thresh * n.sam
    )


  otu.table.sim.filter <-
    otu.table.sim[, otus.keep]


  if (is.null(ref.otu)) {

    ref.otu.filter <- NULL

  } else {

    ref.otu.filter <-
      which(
        colnames(otu.table.sim.filter) ==
          ref.otu
      )
  }


  causal.otus.filter <-
    which(
      otus.keep %in% causal.otus
    )

  noncausal.otus.filter <-
    setdiff(
      seq_along(otus.keep),
      causal.otus.filter
    )


  simulated <- list(

    causal.otus.filter =
      causal.otus.filter,

    noncausal.otus.filter =
      noncausal.otus.filter,

    ref.otu.filter =
      ref.otu.filter,

    otu.table.sim.filter =
      otu.table.sim.filter,

    Y = Y,
    C = C,

    beta0.taxa = beta0.taxa,
    beta.taxa = beta.taxa,
    exp.beta.taxa = exp.beta.taxa,

    sigma = sigma,
    effect.scenario = effect.scenario
  )

  return(simulated)
}
