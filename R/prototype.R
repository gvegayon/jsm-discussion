library(ergmito)
library(ergm)
library(data.table)
library(ggplot2)
library(ggrepel)

source("R/conditional-prob.R")

# Testing the function
# Conditional probabilities
theta <- c(.5, -.5)

# Stats counts
x <- network(matrix(0, 5, 5))
set.vertex.attribute(x, "gender", c(0,1,0,1,0))


# General plotting parameters
parpar <- list(
  mfrow = c(2, 4), mar = rep(.5, 4), oma = c(3,6,0,0),
  col.axis = "gray20", col.lab = "gray20", col.main = "gray20" 
)
width. <- 7 * .9
height. <- 3.5 * .9


# Mutual --------------------------------------
models_mutual <- list(
  x ~ mutual + edges,
  x ~ mutual + ttriad,
  x ~ mutual + nodematch("gender"),
  x ~ mutual + nodeicov("gender")
)

ans_same <- lapply(models_mutual, conditional_dist, theta = c(0, 0))
ans_diff <- lapply(models_mutual, conditional_dist, theta = c(2, 0))

graphics.off()
pdf("figures/conditional-prob-mutuals.pdf", width = width., height = height.)
op <- do.call(par, parpar)
plot_first(ans_same, main = expression( bquote( (a)~theta[mutual] == 0 )) )
plot_first(ans_diff, main = expression( bquote( (b)~theta[mutual] == 2 )) )
par(op)
title(ylab = c("Number of", "Mutual ties"))
dev.off()


# Triads --------------------------------------
models_ttriad <- list(
  x ~ ttriad + edges,
  x ~ ttriad + mutual,
  x ~ ttriad + nodematch("gender"),
  x ~ ttriad + nodeicov("gender")
)

ans_same <- lapply(models_ttriad, conditional_dist, theta = c(0, 0))
ans_diff <- lapply(models_ttriad, conditional_dist, theta = c(1, 0))

graphics.off()
pdf("figures/conditional-prob-ttriad.pdf", width = width., height = height.)
op <- do.call(par, parpar)
plot_first(ans_same, main = expression(  (a)~theta[ttriad] == 0 ))
plot_first(ans_diff, main = expression(  (b)~theta[ttriad] == 1 ), add_xlab = TRUE)
par(op)
title(ylab = c("Number of", "Transitive Triads"))
dev.off()

graphics.off()
png("figures/conditional-prob-ttriad.png", width = width., height = height., unit = "in", res = 150)
op <- do.call(par, parpar)
plot_first(ans_same, main = expression((a)~theta[ttriad] == 0 ))
plot_first(ans_diff, main = expression((b)~theta[ttriad] == 1 ) , add_xlab = TRUE)
par(op)
title(ylab = c("Number of", "Transitive Triads"))
dev.off()

# Gender homophily --------------------------------------
models_homophily <- list(
  x ~ nodematch("gender") + edges,
  x ~ nodematch("gender") + mutual,
  x ~ nodematch("gender") + ttriad,
  x ~ nodematch("gender") + nodeicov("gender")
)

ans_same <- lapply(models_homophily, conditional_dist, theta = c(0, 0))
ans_diff <- lapply(models_homophily, conditional_dist, theta = c(2, 0))

graphics.off()
pdf("figures/conditional-prob-homophily.pdf", width = width., height = height.)
op <- do.call(par, parpar)
plot_first(ans_same, main = expression((a)~theta[homophily] == 0 ) )
plot_first(ans_diff, main = expression((b)~theta[homophily] == 2 ) )
par(op)
title(ylab = c("Number of", "Gender-Homophilic Ties"))
dev.off()

# Gender homophily --------------------------------------
models_icov <- list(
  x ~ nodeicov("gender") + edges,
  x ~ nodeicov("gender") + mutual,
  x ~ nodeicov("gender") + ttriad,
  x ~ nodeicov("gender") + nodematch("gender")
)

ans_same <- lapply(models_icov, conditional_dist, theta = c(0, 0))
ans_diff <- lapply(models_icov, conditional_dist, theta = c(2, 0))

graphics.off()
pdf("figures/conditional-prob-receiver-effect.pdf", width = width., height = height.)
op <- do.call(par, parpar)
plot_first(ans_same, main = expression( (a)~theta[receiver] == 0 ) )
plot_first(ans_diff, main = expression( (b)~theta[receiver] == 2 ) )
par(op)
title(ylab = "Gender-Receiver Effect")
dev.off()

# Focusing on mutuality with transitivity ----------------
models_mutual2 <- list(
  c(-1, 0), c(0, 0), c(1, 0)
)

ans_same <- lapply(models_mutual2, conditional_dist, model = x~ttriad + mutual)

graphics.off()
op <- par(mfrow = c(1, 3))
plot_first(ans_same, main = expression( (a)~theta[receiver] == 0 ) )

par(op)

# Looking into the different values of VIF under different models 
# (and different values of theta)
# VIFs

# coeffs <- combn(seq(from = -1, to = 1, by = 0.25), m = 2, simplify = FALSE)

coeffs <- expand.grid(
  seq(from = -1, to = 1, by = 0.5),
  c(-2, 0, 2) # seq(from = -2, to = 2, by = 1)
) |> apply(1, c, simplify = FALSE)

# coeffs <- list(
#   c(-2, -2), c(-1, -2), c(0, -2), c(1, -2), c(2, -2),
#   c(-2, -1), c(-1, -1), c(0, -1), c(1, -1), c(2, -1),
#   c(-2, 0), c(-1,0), c(0, 0), c(1, 0), c(2, 0),
#   c(-2, 1), c(-1,1), c(0, 1), c(1, 1), c(2, 1),
#   c(-2, 2), c(-1,2), c(0, 2), c(1, 2), c(2, 2)
# )

# Comparing with networks of size 5 -------------------------------

vif5 <- plot_vif_cor(coeffs, model = x~ttriad + mutual, name. = "figures/%s-n=5.pdf")

# Comparing with networks of size 4 -------------------------------
x4 <- network(matrix(0, 4, 4))
vif4 <- plot_vif_cor(coeffs, model = x4~ttriad + mutual, name. = "figures/%s-n=4.pdf")

x3 <- network(matrix(0, 3, 3))
vif3 <- plot_vif_cor(coeffs, model = x3~ttriad + mutual, name. = "figures/%s-n=3.pdf")

# Combining the datasets
vif5[, n := 5]
vif4[, n := 4]
vif3[, n := 3]
vifs <- rbind(vif5, vif4, vif3)

# All with transitive triad
ggplot(vifs, aes(x = theta_mutual, y = vif)) +
  geom_line(aes(color = factor(n))) +
  geom_point(aes(color = factor(n), shape = as.numeric(density), size = 2)) +
  scale_shape_binned() +
  scale_y_log10() +
  facet_wrap(
    ~ theta_ttriad, ncol = 5, scales = "free",
    labeller = label_bquote(theta[ttriad] == .(theta_ttriad))
    ) +
  theme_bw() +
  scale_colour_brewer(palette = "Dark2") +
  annotate(
    "text",
    x = -.4, y = 11,
    label = "VIF = 10"
    ) +
  geom_text_repel(
    aes(x = theta_mutual, y = vif, label = vif_lab),
    size = 3, nudge_y = .1, nudge_x = -.05) +
  labs(
    x = expression(theta[mutual]),
    y = "Variance Inflation Factor",
    color = "Network Size",
    shape = "Density"
  ) +
  # Suppressing size legend
  guides(size = "none")

ggsave("figures/vif.pdf", width = 10, height = 3)

ggplot(vifs, aes(x = theta_mutual, y = cor)) +
  geom_line(aes(color = factor(n))) +
  geom_point(aes(color = factor(n), shape = as.numeric(density), size = 2)) +
  scale_shape_binned() +
  facet_wrap(
    ~ theta_ttriad, ncol = 5,
    labeller = label_bquote(theta[ttriad] == .(theta_ttriad))
    ) +
  theme_bw() +
  scale_colour_brewer(palette = "Dark2") +
  geom_text_repel(
    aes(x = theta_mutual, y = cor, label = cor_lab),
    size = 3, nudge_y = .1, nudge_x = -.05) +
  labs(
    x = expression(theta[mutual]),
    y = "Correlation",
    color = "Network Size",
    shape = "Density",
    size = NULL
  ) +
  # Suppressing size legend
  guides(size = "none")

ggsave("figures/cor.pdf", width = 10, height = 3)
