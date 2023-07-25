#' Generates the conditional probability of observing e given (t, theta.)
#' @noRd
d_ergm <- function(a., b., S., theta.) {
  
  if (length(a.) > 1)
    return(sapply(a., d_ergm, b. = b., theta. = theta., S. = S.))
  
  # Which matches
  idx <- which(colSums(t(S.$statmat) == c(a., b.)) == 2)
  if (!length(idx))
    return(0)
  
  idx0 <- which(S.$statmat[,2] == b.)
  
  # Computing the probability
  ans <- tryCatch({
    sum(S.$weights[idx]) * exp(theta. %*% cbind(c(a., b.)))/
      (rbind(S.$weights[idx0]) %*% exp(S.$statmat[idx0, , drop=FALSE] %*% theta.))[1]
      }, warning = function(w) "NO"
      )

  if (inherits(ans, "character")) {
    stop("There was a warning")
  }

  ans
  
}

#' Generates the c(.05,.5,.95) quantiles for the conditional distribution
#' of term 1 given term 2.
#' @param model A model to be passed to [ergm::ermg.allstats]
#' @param theta Model parameters
conditional_dist <- function(model, theta) {

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

  edges_terms <- with(S., cbind(edges = statmat[, loc], weights = weights))

  # Removing the edges from statmat
  if (removeit)
    S.$statmat <- S.$statmat[, -loc, drop = FALSE]
  
  seq1 <- sort(unique(S.$statmat[,1]))
  seq2 <- sort(unique(S.$statmat[,2]))
  
  # Computing the conditional probability
  ans <- vector("list", length(seq2))
  for (i in seq_along(ans))
    ans[[i]] <- d_ergm(a. = seq1, b. = seq2[i], S. = S., theta. = theta)
  
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
    S.$statmat,
    params        = theta,
    stats_weights = S.$weights,
    stats_statmat = S.$statmat,
    as_prob       = TRUE
    ) * S.$weights
  
  ord <- order(S.$statmat[,2], decreasing = FALSE)
  ergms_probs <- ergm_probs[ord]
  
  b_50pcent <- which.min(abs(cumsum(ergm_probs) - .5))
  b_50pcent <- S.$statmat[ord,][b_50pcent,2]
  b_max     <- S.$statmat[ord,][which.max(ergms_probs),2]

  # Computing weighted OLS
  ols <- with(S., lm(statmat[,1] ~ statmat[,2], weights = ergm_probs))

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
  
  # Drawing a polygon
  graphics::polygon(
    x = scaler(c(x$b, rev(x$b))),
    y = c(x$quantile[,1], rev(x$quantile[,3])),
    col = "lightgray",
    border = "darkgray"
      
  )
  
  lines(
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


plot_vif_cor <- function(coeffs., model., name. = "figures/%s.pdf") { 

  stats_calc <- lapply(
    coeffs., conditional_dist, model = model.
    )

  # Extracting the cor elements and putting them in a df
  cor_df <- do.call(rbind, lapply(stats_calc, function(x) x$cor[1,2]))
  ivf_df <- do.call(rbind, lapply(stats_calc, function(x) {
    1 / (1 - summary(x$ols)$r.squared)
    }))
  den_df <- do.call(rbind, lapply(stats_calc, function(x) {
    x$expected_density
  }))
  ans_same <- cbind(do.call(rbind, coeffs.), cor_df, ivf_df, den_df)

  colnames(ans_same) <- c("theta_mutual", "theta_ttriad", "cor", "vif", "density")
  ans_same <- as.data.table(ans_same)

  ans_same[, density := sprintf("%.2f", density)]
  ans_same[, vif_lab := sprintf("%.2f", vif)]
  ans_same[, cor_lab := sprintf("%.2f", cor)]

  ans_same |>
    ggplot(aes(x = theta_mutual, y = vif, group = theta_ttriad, colour = as.factor(theta_ttriad))) +
      geom_line() +
      geom_point(aes(shape = as.numeric(density))) +
      scale_shape_binned() +
      scale_y_log10() +
      theme_bw() +
      # Color blind colour scale
      scale_colour_brewer(palette = "Dark2") +
      labs(
        x = expression(theta[mutual]),
        y = "VIF (log10 scale)",
        colour = expression(theta[ttriad]),
        shape = "Density"
        ) +
      geom_hline(yintercept = 10, linetype = "dashed") +
      annotate(
        "text",
        x = -.75, y = 15,
        label = "VIF = 10"
        ) +
      geom_text_repel(
        aes(x = theta_mutual, y = vif, label = vif_lab),
        size = 3, nudge_y = .1, nudge_x = -.05) +
      # Label the figure as (a)
      annotate(
        "text",
        x = -.75, y = max(ans_same$vif) * 1.1,
        label = "(a)",
        size = 5
        )

  ggsave(sprintf(name., "vif"), width = 3.5, height = 3.5)
  ggsave(sprintf(gsub("\\.pdf", ".png", name.), "vif"), width = 4, height = 3)

  # Same figure, but instead of VIF use correlation
  ans_same |>
    ggplot(aes(x = theta_mutual, y = cor, group = theta_ttriad, colour = as.factor(theta_ttriad))) +
      geom_line() +
      geom_point(aes(shape = as.numeric(density))) +
      scale_shape_binned() +
      theme_bw() +
      scale_colour_brewer(palette = "Dark2") +
      labs(
        x = expression(theta[mutual]),
        y = "Correlation",
        colour = expression(theta[ttriad]),
        shape = "Density"
        ) +
      geom_hline(yintercept = 0.5, linetype = "dashed") +
      annotate(
        "text",
        x = -.75, y = .55,
        label = "Cor = 0.5" 
        ) +
      geom_text_repel(
        aes(x = theta_mutual, y = cor, label = cor_lab),
        size = 3) +
      # Label the figure as (a)
      annotate(
        "text",
        x = -.75, y = max(ans_same$cor) * 1.1,
        label = "(b)",
        size = 5
        )

  ggsave(sprintf(name., "cor"), width = 3.5, height = 3.5)
  ggsave(sprintf(gsub("\\.pdf", ".png", name.), "cor"), width = 4, height = 3)

  invisible(ans_same)

}