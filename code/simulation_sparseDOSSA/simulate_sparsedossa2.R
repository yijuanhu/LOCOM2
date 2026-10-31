# Help functions for the SparseDOSSA2 simulations of the LOCOM2 paper

library(SparseDOSSA2)

#############################################
################# simulator #################
#############################################

simulate_sparsedossa2 = function(beta, seed, n.sam, causal.type, filter.thresh,
                                 ref.otu = NULL, depth.mu = 10^4, depth.fix = NULL,
                                 diff.lib = 0, template = "Stool", beta.c = 0,
                                 prop.case = 0.5,
                                 spike.property = c("abundance", "prevalence", "both", "multiplicative"),
                                 spike.property.c = c("abundance", "prevalence", "both", "multiplicative"),
                                 pi.max = 0.05) {

  spike.property = match.arg(spike.property)
  spike.property.c = match.arg(spike.property.c)
  if (spike.property == "multiplicative") {
    spike.property.c = "multiplicative"
  }
  if (spike.property.c == "multiplicative" && spike.property != "multiplicative") {
    stop("spike.property.c = \"multiplicative\" requires spike.property = \"multiplicative\"")
  }

  if (template == "momspi") {
    template = "Vaginal"
  }
  if (!(template %in% c("Stool", "Vaginal", "IBD"))) {
    stop("template must be \"Stool\", \"Vaginal\" (or \"momspi\") or \"IBD\"")
  }

  # template mean relative abundances; the taxa are renamed f1, ..., fJ
  pi_taxa = sd2_pi_taxa(template)
  n.otus = length(pi_taxa)
  features = names(pi_taxa)
  taxa.names = paste0("f", 1:n.otus)

  set.seed(seed)
  n0.sam = as.integer(round(n.sam * (1-prop.case)))
  n1.sam = as.integer(round(n.sam * prop.case))
  Y <- c(rep(0, n0.sam), rep(1, n1.sam))
  C = NULL
  C.index = NULL

  # reference taxon (index): the most abundant taxon by default (the role of taxon 4414 in the throat template)
  if (is.null(ref.otu)) {
    ref.otu = which.max(pi_taxa)
  } else if (is.character(ref.otu)) {
    ref.otu = ifelse(ref.otu %in% taxa.names, match(ref.otu, taxa.names), match(ref.otu, features))
  }

  if (template %in% c("Stool", "IBD")) {

    # throat-type rules: causal.type 1 = 20 abundant non-reference taxa; 2 = the five most abundant taxa
    # (includes the reference taxon); 3 and 4 = 20 (9) abundant non-reference taxa plus the reference
    # taxon (taxon 468 in the throat rules)
    if (beta.c == 0) {

      if (causal.type == 1){
        w = setdiff(which(pi_taxa >= 0.005 & pi_taxa < pi.max), ref.otu)
        causal.otus <- pick_causal(w, 20, pi_taxa, pi.max, ref.otu, "causal.type 1")
      }
      if (causal.type == 2){
        causal.otus <- order(pi_taxa, decreasing = TRUE)[1:5]
      }
      if (causal.type == 3) {
        w = setdiff(which(pi_taxa >= 0.001 & pi_taxa < pi.max), ref.otu)
        causal.otus <- c(pick_causal(w, 20, pi_taxa, pi.max, ref.otu, "causal.type 3"), ref.otu)
      }
      if (causal.type == 4) {
        w = setdiff(which(pi_taxa >= 0.005 & pi_taxa < pi.max), ref.otu)
        causal.otus <- c(pick_causal(w, 9, pi_taxa, pi.max, ref.otu, "causal.type 4"), ref.otu)
      }

    } else {

      # Y affects the taxa in Y.index and C the taxa in C.index; the two sets partly overlap
      if (causal.type == 1){
        w = setdiff(which(pi_taxa >= 0.005 & pi_taxa < pi.max), ref.otu)
        w <- pick_causal(w, 30, pi_taxa, pi.max, ref.otu, "causal.type 1 with C")
      }
      if (causal.type == 2){
        w = setdiff(which(pi_taxa >= 0.002 & pi_taxa < pi.max), ref.otu)
        w <- pick_causal(w, 30, pi_taxa, pi.max, ref.otu, "causal.type 2 with C")
      }
      if (causal.type == 3){
        w = setdiff(which(pi_taxa >= 0.001 & pi_taxa < pi.max), ref.otu)
        w <- pick_causal(w, 30, pi_taxa, pi.max, ref.otu, "causal.type 3 with C")
      }
      if (causal.type %in% 1:3) {
        Y.index <- c(w[1:10], w[11:20])
        C.index <- c(w[21:30], w[11:20])
      }
      if (causal.type == 4){
        o <- order(pi_taxa, decreasing = TRUE)
        Y.index <- c(o[1:3], o[4:5])
        C.index <- c(o[31:33], o[4:5])
      }
      causal.otus <- Y.index

    } # beta.c == 0?

  } # Stool, IBD

  if (template == "Vaginal") {

    # momspi-type rules
    if (beta.c == 0) {
      if (causal.type == 1){
        causal.otus <- c(which(pi_taxa >= 0.002 & pi_taxa < 0.2))[c(1:17, 19:21)]
      }
      if (causal.type == 2){
        w = which(pi_taxa >= 0.002 & pi_taxa < 0.2)
        causal.otus <- pick_causal(w, 20, pi_taxa, pi.max, ref.otu, "causal.type 2")
      }
      if (causal.type == 3){
        causal.otus <- order(pi_taxa, decreasing = TRUE)[2:6]
      }
      if (causal.type == 4){
        causal.otus <- order(pi_taxa, decreasing = TRUE)[3:7]
      }
    } else {
      w = which(pi_taxa >= 0.0005 & pi_taxa < 0.2)
      w <- pick_causal(w, 30, pi_taxa, pi.max, ref.otu, "causal.type with C")
      Y.index <- c(w[1:10], w[11:20])
      C.index <- c(w[21:30], w[11:20])
      causal.otus <- Y.index
    }

  } # Vaginal

  if (beta.c != 0) {
    # C is independent of Y
    rho <- 0.5
    C <- c(rbinom(n0.sam, 1, rho), rbinom(n1.sam, 1, 1 - rho))
  }

  exp.beta = exp(beta)
  noncausal.otus <- setdiff(1:n.otus, causal.otus)

  # Library sizes (as in simulate_midasim())
  #   diff.lib = 0: normal with mean depth.mu (or fixed at depth.fix)
  #   diff.lib > 0: case:control mean ratio 1.5:1 (diff.lib = 1), 3:1 (2), 5:1 (3) or 1:1 (4)
  if (diff.lib == 0) {

    if (is.null(depth.fix)) {
      depth.lower = 2000    # lower bound for library size
      depth.sd = depth.mu/3 # SD of library size

      depth.final <- rnorm(n.sam, depth.mu, depth.sd)
      depth.final[depth.final < depth.lower] <- depth.lower
      depth.final <- round(depth.final)
    } else {
      depth.final = rep(depth.fix, n.sam)
    }

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
    depth.final <- c(round(rnorm(n0.sam, mean.lib.size.control, mean.lib.size.control/3)),
                     round(rnorm(n1.sam, mean.lib.size.case, mean.lib.size.case/3)) )
    depth.final <- pmax(depth.final, 2000)

  } # diff.lib == 0?

  # relative abundances from SparseDOSSA2 with the Y (and C) effects; metadata column 1 = Y, 2 = C
  metadata = cbind(Y = Y)
  if (!is.null(C)) {
    metadata = cbind(metadata, C = C)
  }
  rownames(metadata) = paste0("Sample", 1:n.sam)

  if (spike.property == "multiplicative") {

    # MIDASim analogue: null draw, then the latent absolute abundances of the causal taxa are multiplied
    # by exp(beta * Y) (and those of the C taxa by exp(beta.c * C)), then TSS
    sim = SparseDOSSA2(template = template, n_sample = n.sam, new_features = FALSE,
                       spike_metadata = "none", median_read_depth = 20000, verbose = FALSE)
    a = sim$simulated_matrices$a_null[features, , drop = FALSE]   # taxa x samples
    if (beta != 0) {
      a[causal.otus, ] = a[causal.otus, , drop = FALSE] *
        matrix(exp(beta * Y), nrow = length(causal.otus), ncol = n.sam, byrow = TRUE)
    }
    if (beta.c != 0) {
      a[C.index, ] = a[C.index, , drop = FALSE] *
        matrix(exp(beta.c * C), nrow = length(C.index), ncol = n.sam, byrow = TRUE)
    }
    otu.ra.sim = t(apply(a, 2, function(x) if (sum(x) == 0) x else x / sum(x)))   # samples x taxa

  } else {

    # SparseDOSSA2's own spike-in of the effects into the latent absolute abundances
    spike = NULL
    if (beta != 0) {
      spike = rbind(spike, sd2_spike_table(1, causal.otus, beta, spike.property, features))
    }
    if (beta.c != 0) {
      spike = rbind(spike, sd2_spike_table(2, C.index, beta.c, spike.property.c, features))
    }
    if (is.null(spike)) {
      spike = "none"
    }
    sim = SparseDOSSA2(template = template, n_sample = n.sam, new_features = FALSE,
                       spike_metadata = spike, metadata_matrix = metadata,
                       median_read_depth = 20000, verbose = FALSE)
    otu.ra.sim = t(sim$simulated_matrices$rel)[, features, drop = FALSE]   # samples x taxa

  } # spike.property == "multiplicative"?

  # counts: multinomial draws at the library sizes depth.final
  otu.table.sim <- matrix(0, nrow = n.sam, ncol = n.otus)
  for (i in 1:n.sam) {
    if (sum(otu.ra.sim[i, ]) > 0) {
      otu.table.sim[i, ] <- as.vector(rmultinom(n = 1, size = depth.final[i], prob = otu.ra.sim[i, ]))
    }
  }
  colnames(otu.table.sim) <- taxa.names

  # keep taxa present in at least filter.thresh * n.sam samples
  prop.count <- colSums(otu.table.sim > 0)
  otus.keep <- which(prop.count >= filter.thresh * n.sam)

  otu.table.sim.filter <- otu.table.sim[, otus.keep, drop = FALSE]

  # a very shallow sample can lose all of its reads to the filtered taxa; such samples are removed
  ok = rowSums(otu.table.sim.filter) > 0
  if (!all(ok)) {
    otu.table.sim.filter = otu.table.sim.filter[ok, , drop = FALSE]
    Y = Y[ok]
    if (!is.null(C)) {
      C = C[ok]
    }
  }

  ref.otu.filter <- which(otus.keep == ref.otu)
  if (length(ref.otu.filter) == 0) {
    ref.otu.filter = NULL
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


# mean relative abundance of each taxon of a template (the role of pi_taxa in the MIDASim simulations):
# empirical means in a null draw of n.null samples with a fixed seed, so every replicate uses the same values

sd2_pi_taxa = function(template, n.null = 2000, seed = 20240101) {

  set.seed(seed)
  sim0 = SparseDOSSA2(template = template, n_sample = n.null, new_features = FALSE,
                      spike_metadata = "none", median_read_depth = 20000, verbose = FALSE)
  pi_taxa = rowMeans(sim0$simulated_matrices$rel)
  names(pi_taxa) = rownames(sim0$simulated_matrices$rel)

  return(pi_taxa)
}


# the first k taxa of a candidate set; if the set has fewer than k taxa, it is topped up with the most
# abundant remaining taxa below pi.max (excluding the reference taxon), with a message

pick_causal = function(candidates, k, pi_taxa, pi.max, ref.otu, label) {

  if (length(candidates) < k) {
    extra = setdiff(order(pi_taxa, decreasing = TRUE), c(candidates, ref.otu))
    extra = extra[pi_taxa[extra] < pi.max]
    if (length(candidates) + length(extra) < k) {
      stop("Not enough taxa for ", label)
    }
    message(label, ": only ", length(candidates), " taxa in the threshold set; topped up to ", k,
            " with the most abundant remaining taxa (pi < pi.max)")
    candidates = c(candidates, extra[1:(k - length(candidates))])
  }

  return(candidates[1:k])
}


# spike-in table for SparseDOSSA2: one row per spiked taxon and property
# (metadata_datum = column of the metadata matrix, 1 = Y and 2 = C)

sd2_spike_table = function(datum, taxa.index, effect, property, features) {

  if (property == "both") {
    properties = c("prevalence", "abundance")
  } else {
    properties = property
  }

  spike = NULL
  for (p in properties) {
    spike = rbind(spike, data.frame(metadata_datum = as.integer(datum),
                                    feature_spiked = features[taxa.index],
                                    associated_property = p,
                                    effect_size = effect,
                                    stringsAsFactors = FALSE))
  }

  return(spike)
}
