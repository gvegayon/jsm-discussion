#' Generates the conditional probability of observing e given (t, theta.)
#' @noRd
d_ergm <- function(a., b., S., theta., edges. = 0, theta_edges. = 0) {
  
  if (length(a.) > 1)
    return(sapply(
      a., d_ergm, b. = b., theta. = theta.,
      edges. = edges.,
      theta_edges. = theta_edges., S. = S.)
      )
  
  # Which matches
  idx <- which(colSums(t(S.$statmat) == c(a., b.)) == 2)
  if (!length(idx))
    return(0)
  
  idx0 <- which(S.$statmat[,2] == b.)

  # if ("ttriple" %in% colnames(S.$statmat) & any(theta. != 0)) {
  #   idx0 <- idx0
  # }
  
  # Computing the probability
  ans <- tryCatch({
    exp(
      c(theta., theta_edges.) %*% t(cbind(a., b., edges.[idx]))
      ) %*% S.$weights[idx] /
      (rbind(S.$weights[idx0]) %*% exp(
        cbind(S.$statmat, edges.)[idx0, , drop=FALSE] %*%
        c(theta., theta_edges.)
        ))[1]
      }, error = function(w) w
      )

  if (inherits(ans, "error")) {
    stop("There was a warning")
  }

  ans
  
}

#' Generates the c(.05,.5,.95) quantiles for the conditional distribution
#' of term 1 given term 2.
#' @param model A model to be passed to [ergm::ermg.allstats]
#' @param theta Model parameters
conditional_dist <- function(model, theta, theta_edges. = 0) {

  # Checking if it includes edgecounts
  loc <- which(attr(terms(model), "term.labels") == "edges")
  if (!length(loc)) {
    model_new <- update.formula(model, . ~ . + edges)
    loc <- length(attr(terms(model_new), "term.labels"))
    removeit <- TRUE
  } else {
    model_new <- model
    removeit <- FALSE
  }

  S. <- ergm::ergm.allstats(model_new)

  edges_terms <- with(S., cbind(
    edges = statmat[, loc], weights = weights
    ))

  # Removing the edges from statmat
  if (removeit)
    S.$statmat <- S.$statmat[, -loc, drop = FALSE]
  
  seq1 <- sort(unique(S.$statmat[,1]))
  seq2 <- sort(unique(S.$statmat[,2]))
  
  # Computing the conditional probability
  ans <- vector("list", length(seq2))
  for (i in seq_along(ans)) {

    ans[[i]] <- d_ergm(
      a. = seq1, b. = seq2[i], S. = S.,
      theta.       = theta,
      theta_edges. = theta_edges.,
      edges.       = edges_terms[,1]
      )

    # Checking if the probabilities sum to 1
    if (abs(sum(ans[[i]]) - 1) > .00001) {
      stop("The probabilities do not sum to 1")
    }

  }

  # Calculating the quantiles
  quantiles <- sapply(ans, function(a) {
    
    a <- cumsum(a)
    
    # Identifying which are zero
    are0 <- which(a < 0.000000001)
    
    if (length(are0)) {
      a <- a[-are0]
      e <- seq1[-are0]
    } else
      e <- seq1
    
    if (!length(e))
      return(c(NA, NA, NA))
    
    e[c(
      # 5% 
      which.min(abs(a - .05)),
      # 50%
      which.min(abs(a - .5)),
      # %95
      which.min(abs(a - .95))
    )]
    
  })

  quantiles <- t(quantiles)
  
  # Which quantile of the conditional is closer to the 50%
  ergm_probs <- ergmito::exact_loglik(
    if (removeit) 
      cbind(S.$statmat, edges_terms[,1]) else S.$statmat,
    params        = if (removeit)
      c(theta, theta_edges.) else theta,
    stats_weights = S.$weights,
    stats_statmat = if (removeit)
      cbind(S.$statmat, edges_terms[,1]) else S.$statmat,
    as_prob       = TRUE
    ) * S.$weights

  # The sum of ergm_probs should be close to 1
  if (abs(sum(ergm_probs) - 1) > .00001) {
    stop("The probabilities do not sum to 1")
  }
  
  ord <- order(S.$statmat[,2], decreasing = FALSE)
  ergms_probs <- ergm_probs[ord]
  
  b_50pcent <- which.min(abs(cumsum(ergm_probs) - .5))
  b_50pcent <- S.$statmat[ord,][b_50pcent,2]
  b_max     <- S.$statmat[ord,][which.max(ergms_probs),2]

  # Computing weighted OLS
  ols <- if (removeit) {
    with(S., lm(statmat[,1] ~ statmat[,2] + edges_terms[,1], weights = ergm_probs))
  } else {
    with(S., lm(statmat[,1] ~ statmat[,2], weights = ergm_probs))
  }

  # Computing weighted correlation coefficient
  cor. <- with(S., cov.wt(statmat, wt = ergm_probs, cor = TRUE))
  
  # Calculating expected edgecount and density
  n <- network.size(eval(model[[2]]))
  expected_edgecount <- sum(edges_terms[,1] * ergm_probs)
  expected_density   <- expected_edgecount/(n * (n-1))

  structure(
    list(
      a = seq1,
      b = seq2,
      pr = ans,
      quantiles = cbind(
        `5%`= quantiles[,1],
        `50%` = quantiles[,2],
        `95%` = quantiles[,3]
        ),
      par_names = colnames(S.$statmat),
      probs2    = ergm_probs,
      b_50pcent = b_50pcent,
      b_max     = b_max,
      support   = S.,
      ols       = ols,
      cor       = cor.$cor,
      edges     = edges_terms,
      expected_edgecount = expected_edgecount,
      expected_density   = expected_density
      ),
    class = "ergm_conditional"
  )
}

namer <- function(x) {
  switch(
    x, 
    edges = "Edge count",
    ttriple = "Transitive Triads",
    mutual = "Mutual Ties",
    nodeicov.gender = "Gender-Receiver",
    nodematch.gender = "Gender-Homophily"
  )
}

#' Plot method for conditional probabilities on ERGMs
#' @noRd
plot.ergm_conditional <- function(
  x,
  xlab = namer(x$par_names[2]),
  ylab = namer(x$par_names[1]),
  ...
  ) {
  
  scaler <- function(a) a
  
  y_range <- range(x$a) 
  y_range <- y_range + c(-.2, .1) * diff(y_range)/10
  x_range <- scaler(range(x$b)) 
  
  graphics::plot(
    NA,
    xlim = x_range,
    ylim = y_range,
    xlab = xlab,
    ylab = ylab, xaxt = "n", yaxt = "n"
    )

  op <- graphics::par(xpd = FALSE)
  graphics::grid()
  graphics::par(op)

  # Checking the probs2 add up to one
  if (abs(sum(x$probs2) - 1) > .00001) {
    stop("The probabilities do not sum to 1")
  }

  # Sampling
  set.seed(3131)
  sampled <- with(
    x, sample(
      1:length(probs2),
      prob = probs2,
      size = 1000, replace = TRUE
      )
  )

  graphics::smoothScatter(
    x = x$support$statmat[,2][sampled],
    y = x$support$statmat[,1][sampled],
    colramp = colorRampPalette(c("transparent", "black")),
    add = TRUE,
    nrpoints = Inf
    )
  
  # Drawing a polygon
  graphics::polygon(
    x = scaler(c(x$b, rev(x$b))),
    y = c(x$quantile[,1], rev(x$quantile[,3])),
    col    = adjustcolor("lightgray", alpha.f = .5),
    border = "darkgray"
  )
  
  graphics::lines(
    y    = x$quantiles[,2],
    x    = scaler(x$b),
    type = "l", col = "red", lwd = 2, lty = 2
  )

  # Adding the VIF
  vif <- 1/(1 - summary(x$ols)$adj.r.squared)
  text(
    x = x_range[1] + diff(x_range) * .2,
    y = y_range[2] - diff(y_range) * .1,
    # Printing R2 using mathplot
    labels = substitute(
      VIF == .r2,
      list(.r2 = sprintf("%.2f", vif))
      ),
    pos = 4
  )

  rho_val <- sprintf("%.2f", x$cor[1, 2])
  rho_val <- if (rho_val == "1.00") {
    expression(
      rho %->% 1
      )
  } else {
    substitute(
      rho == .rho,
      list(.rho = rho_val)
      )
  }

  text(
    x = x_range[1] + diff(x_range) * .2,
    y = y_range[2] - diff(y_range) * .2,
    # Printing R2 using mathplot
    labels = rho_val,
    pos = 4
  )

  # # Adding the expected density
  # text(
  #   x = x_range[1] + diff(x_range) * .2,
  #   y = y_range[2] - diff(y_range) * .3,
  #   # Printing R2 using mathplot
  #   labels = substitute(
  #     Avg.~Density == .rho,
  #     list(.rho = sprintf("%.2f", x$expected_density))
  #     ),
  #   pos = 4
  # )
  
  # abline(v = x$b_50pcent, lty = 2, lwd = 2)
  opp <- par(tck = .07, las = 1,
    yaxp = c(
      round(c(y_range  + diff(y_range)/5 * c(1, 0))/ 10) * 10 , 5)
    )
  graphics::axis(side = 1, line = -2.5, lwd = 0)
  graphics::axis(side = 1, line = 0, lwd = 0, lwd.ticks = 1, labels = FALSE)
  graphics::axis(side = 2, line = -3, lwd = 0)
  graphics::axis(side = 2, line = 0, lwd = 0, lwd.ticks = 1, labels = FALSE)
  par(opp)
  
  return()
  
}

# Plot the first with a title
plot_first <- function(x, main="", add_xlab = FALSE, xlab = "", ylab = "",...) {
      op <- par(xpd=NA)
  for (i in seq_along(x)) {
    if (i == 1L) {
      plot(x[[i]], xlab = xlab, ylab = ylab, ...)
      # title(main = main, line = 1.25, font.main = 1, cex.main = 1.5, adj = 0)
      graphics::mtext(text = main, side = 2, cex = .75, line = 1)
      par(op)
    } else {
      plot(x[[i]], xlab = xlab, ylab = ylab,  ...)
    }

    if (add_xlab) {

      graphics::mtext(text = namer(x[[i]]$par_names[2]), side = 1, line = 1, cex = .75)

    }
  }

}


plot_vif_cor <- function(coeffs., model., name. = "figures/%s.pdf", ...) { 

  stats_calc <- lapply(
    coeffs., \(coeff.) {
      conditional_dist(model = model., theta = coeff.[-1], theta_edges. = coeff.[1])
    })

  # Extracting the cor elements and putting them in a df
  cor_df <- do.call(rbind, lapply(stats_calc, function(x) x$cor[1,2]))

  ivf_df <- do.call(rbind, lapply(stats_calc, function(x) {
    1 / (1 - summary(x$ols)$r.squared)
    }))

  den_df <- do.call(rbind, lapply(stats_calc, function(x) {
    x$expected_density
  }))

  # Removing the edgecount
  ecoefs <- sapply(coeffs., "[", 1)
  coeffs. <- lapply(coeffs., function(x) x[-1])

  ans_same <- cbind(do.call(rbind, coeffs.), cor_df, ivf_df, den_df, ecoefs)

  colnames(ans_same) <- c(
    sprintf("theta_%s", attr(terms(model.), which = "term.labels")),
    c("cor", "vif", "density", "theta_edges")
  )
  ans_same <- as.data.table(ans_same)

  # Turning density into brackets (5)
  ans_same[, density_num := density]
  ans_same[, density := cut(
    density, breaks = {
      seq(from = floor(min(density) * 10)/10, to = 1, length.out = 7)
  }, include.lowest = TRUE)]

  ans_same[, vif_lab := sprintf("%.2f", vif)]
  ans_same[, cor_lab := sprintf("%.2f", cor)]

  print(ans_same)

  ans_same |>
    ggplot(aes(x = theta_mutual, y = vif, group = theta_ttriad, colour = as.factor(theta_ttriad))) +
      geom_line() +
      geom_point(aes(shape = density)) +
      scale_y_log10() +
      theme_bw() +
      # Color blind colour scale
      scale_colour_brewer(palette = "Dark2") +
      labs(
        x = expression(theta[mutual]),
        y = "VIF (log10 scale)",
        colour = expression(theta[ttriad]),
        shape = "Density" #,
        # tag = "a"
        ) +
      geom_hline(yintercept = 150, linetype = "dashed") +
      annotate(
        "text",
        x = min(ans_same$theta_mutual) + .1, y = 160,
        label = "VIF = 150"
        ) +
      geom_text_repel(
        aes(x = theta_mutual, y = vif, label = vif_lab),
        size = 3, nudge_y = .1, nudge_x = -.05)

  ggsave(sprintf(name., "vif"), width = 5, height = 3)
  ggsave(sprintf(gsub("\\.pdf", ".png", name.), "vif"), width = 5, height = 2.5)

  # Same figure, but instead of VIF use correlation
  ans_same |>
    ggplot(aes(x = theta_mutual, y = cor, group = theta_ttriad, colour = as.factor(theta_ttriad))) +
      geom_line() +
      geom_point(aes(shape = density)) +
      theme_bw() +
      scale_colour_brewer(palette = "Dark2") +
      labs(
        x = expression(theta[mutual]),
        y = "Correlation",
        colour = expression(theta[ttriad]),
        shape = "Density", tag = "(b)"
        ) +
      # geom_hline(yintercept = 0.5, linetype = "dashed") +
      # annotate(
      #   "text",
      #   x = min(ans_same$theta_mutual), y = .55,
      #   label = "Cor = 0.5" 
      #   ) +
      geom_text_repel(
        aes(x = theta_mutual, y = cor, label = cor_lab),
        size = 3) 

  ggsave(sprintf(name., "cor"), width = 3.5, height = 3.5)
  ggsave(sprintf(gsub("\\.pdf", ".png", name.), "cor"), width = 4, height = 3)

  invisible(ans_same)

}