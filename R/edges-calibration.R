library(ergm)
library(network)

# Stats counts
x <- network(matrix(0, 5, 5))
set.vertex.attribute(x, "gender", c(0,1,0,1,0))

set.seed(1231)

nsims <- 1e4
thetas <- cbind(runif(nsims, -4, 0), matrix(runif(2 * nsims, -2, 2), ncol = 2))

# Sampling ergms
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

# Fitting a model to predict the edges coefficientas a function of thethas and edgecount
if (FALSE) {
  fit <- glm(thetas[,1] ~ thetas[,2] + thetas[,3] + 
    samples[,1], family = "gaussian")

  plot(
    predict(fit),
    thetas[,1]
  )
  mean(abs(predict(fit) - thetas[,1]))
}

# Now we fit a neural network to predict the edges coefficient as a function of thetas and edgecount
library(nnet)

fit <- nnet(thetas[,1] ~
  thetas[,2] + thetas[,3] + samples[,1] +
  thetas[,2] * samples[,1] + thetas[,3] * samples[,1] +
  thetas[,2] * thetas[,3] + thetas[,2] * thetas[,3] * samples[,1],
  size = 2, linout = TRUE)



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