library(ergm)
library(network)

# Stats counts
x <- network(matrix(0, 5, 5))
set.vertex.attribute(x, "gender", c(0,1,0,1,0))

set.seed(1231)

nsims <- 1e4
thetas <- cbind(
  runif(nsims, -3, 0), 
  runif(nsims, .5, 1),
  runif(nsims, .5, 1)
)

# Sampling ergms
fn <- "R/edges_calibration_sims.rds"
if (!file.exists(fn)) {

  samples <- parallel::mclapply(
    seq_len(nsims), function(i) {

    # Progress bar every 20 simulations
    if (i %% 20 == 0)
      cat(i, "\n")

    simulate_formula(
      x ~ edges + mutual + ttriad, nsim = 1, 
      coef = thetas[i,], output = "stats"
      )
    
  }, mc.cores = 8)

  samples <- do.call(rbind, samples)

  saveRDS(samples, file = fn)

} else {

  samples <- readRDS(fn)

}

# Subsetting for cases where samples[,1] < 20
idx <- samples[,1] < 20 & samples[,1] > 0
thetas <- thetas[idx,]
samples <- samples[idx,]

hist(samples[,1])

modeldata <- data.frame(
  theta_edges = thetas[,1],
  theta_mutual = thetas[,2],
  theta_ttriad = thetas[,3],
  edges        = samples[,1]
)

# Fitting a model to predict the edges coefficientas a function of thethas and edgecount
fit2 <- glm(theta_edges ~ 
  # theta_mutual + theta_ttriad +
  # edges +
  I(edges^2) + 
  I(1/edges) + 
  I(1/edges^2)+
  I(1/theta_ttriad) +
  # I(1/theta_mutual) + 
  I(theta_ttriad/edges)+ 
  theta_mutual * edges + theta_ttriad * edges,
  # theta_mutual * theta_ttriad * edges,
  data = modeldata,
  family = "gaussian");summary(fit2);mean(abs(predict(fit2) - thetas[,1]))

plot(
  predict(fit2),
  modeldata$theta_edges
)

saveRDS(fit2, file = "R/edges_calibration_ols.rds")

# Calibration using NNEts was too wild
if (FALSE) {
  # Now we fit a neural network to predict the edges coefficient as a function of thetas and edgecount
  library(nnet)

  # initial <- readRDS("R/edges_calibration_weights.rds")

  modeldata <- data.frame(
    theta_edges = thetas[,1],
    theta_mutual = thetas[,2],
    theta_ttriad = thetas[,3],
    edges        = samples[,1]
  )

  fit <- nnet(theta_edges ~
    theta_mutual + theta_ttriad +
    edges + edges ^ 2 + log(edges) + sqrt(edges) +
    # samples[, 3] + samples[,3] ^ 2 + sqrt(samples[,3]) +
    theta_mutual * edges + theta_ttriad * edges +
    theta_mutual * theta_ttriad + theta_mutual * theta_ttriad * edges,
    data = modeldata,
    size = 150, linout = TRUE, MaxNWts = 2e4, maxit = 4000) #, Wts = fit$wts)

  # Computing the mean absolute error
  mean(abs(predict(fit) - thetas[,1]))
  plot(
    predict(fit),
    thetas[,1]
  )

  saveRDS(fit, file = "R/edges_calibration.rds")


  if (FALSE) {
    # Fit a neural network using tensorflow
    library(tensorflow)
    library(keras)

    X <- array_reshape(
      c(thetas[,2], thetas[,3], samples[,1],
        thetas[,2] * samples[,1], thetas[,3] * samples[,1]),
      dim = c(nsims, 5)
    ) 

    y <- thetas[,1]

    # Define the model
    model <- keras_model_sequential() |>
      layer_dense(units = 5, activation = "linear", input_shape = c(5)) |>
      layer_dense(units = 4, activation = "relu") |>
      layer_dense(units = 1, activation = "linear")

    # Compile the model
    model |> compile(
      loss = "mean_squared_error",
      optimizer = optimizer_adam(lr = 0.01)
    )

    # Fit the model
    model |> fit(
      X, y, epochs = 100, batch_size = 32, verbose = 0
    )

    # Plot the results
    plot(
      predict(model, X),
      y
    )

    # Compute the mean absolute error
    mean(abs(predict(model, X) - y))
  }

}