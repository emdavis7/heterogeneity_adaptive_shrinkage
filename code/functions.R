#######
# Functions for Simulation
#######
library(Matrix)
library(dplyr)
library(numDeriv)
library(purrr)
library(FEShR)
################################### Auxiliary Functions ##########################################3

### Create auxiliary elements that will help with computation. #########

### Trace: ####
tr <- function(x){
  sum(diag(x))
}


### Combined MSE ####
mse_mle_combined <- function(scaled_XTX, beta, p, D, kron){
  A <-  solve(t(kron) %*% scaled_XTX %*% kron) %*%
    t(kron) %*% scaled_XTX
  combined_bias <- (kron %*% A - diag(1, nrow = p*D)) %*% beta
  combined_var <- solve(t(kron) %*% scaled_XTX %*% kron)
  t(combined_bias) %*% combined_bias + D*tr(combined_var)
}



########################################################
## Functions to create e_i and e_i e_i'
## e_i e_i' : selBlock, a D times D matrix with 1 in the (i,i) position and 0 elsewhere
## e_i : e_i, D times 1 with 1 in ith pos and 0 elsewhere.
########################################################

## e_i e_i' : selBlock, a function returning D times D matrix with 1 in the (i,i) position and 0 elsewhere
selBlock <- function(i, D){
  replace(numeric(D), i, 1) %*% t(replace(numeric(D), i, 1))
}

## e_i : e_i, D times 1 with 1 in ith pos and 0 elsewhere.
e_i <- function(i, D){
  replace(numeric(D), i, 1)
}


### Estimate sigma_sq for each dataset. #####
# If sigma_sq is unknown, estimate for each study using sum of square errors.
est_sigma_sq <- function(data, p){
  x_matrix <- data$X
  y <- data$y
  residual <- (diag(nrow(x_matrix)) - x_matrix %*%
                 solve(t(x_matrix) %*% x_matrix) %*%
                 t(x_matrix)) %*% y
  residual_sq <- residual^2
  colnames(residual_sq) <- "residual_sq"
  tmp <- data.frame(residual_sq, study = data$study)
  tmp2 <- tmp %>% group_by(study) %>% summarize(sse = sum(residual_sq),
                                                adjust = n() - p)
  return(tmp2$sse/tmp2$adjust)
}

####### Lambda Block ####
# Takes a D by 1 lambda vector and turns it into a block diagonal matrix pD by pD
blockLambda_fcn <- function(lambda, p=4, D){
    # Diagonalize - uses hadamard product to keep only the diagonal.
    Lambda <- lambda %*% matrix(data = rep(1, D), nrow = 1, ncol = D) * diag(D)
    kronecker(Lambda, diag(p))
}


###### True MSE #######
mse <- function(par, scaled_XTX, beta, D, p, kron, par_type = "pi"){
  if (par_type == "lambda") {
    pi <- par/(par+1)
  } else {
    pi <- par
  }

  # Block diagonal pi.
  bPi <- blockLambda_fcn(pi, p = p, D = D)


  if (all.equal(pi, rep(0, D))==T){tr(solve(scaled_XTX))
  } else{

    # Theta maker:
    A <- solve( t(kron) %*% bPi %*% scaled_XTX %*% kron) %*%
      t(kron) %*% bPi %*% scaled_XTX
    # bias_comp
    b_comp <- bPi %*% (kron %*% A - diag(1, p*D))
    bias_sq <- t(beta) %*% t(b_comp) %*% b_comp %*% beta

    var_mle <- tr(solve(scaled_XTX))
    var_betahat <- tr( t(b_comp) %*% b_comp %*% solve(scaled_XTX))
    cov_piece <- 2*tr(
      solve(scaled_XTX) %*% t(b_comp)
    )

    # Return the sum.
    var_mle +
      cov_piece +
      var_betahat +
      bias_sq
  }
}


###### SURE MSE #######
# Parameterizing with pi for now.
sure_mse <- function(par, scaled_XTX, indv_mle, D, p, kron, par_type = "pi"){
  if (par_type == "lambda") {
    pi <- par/(par+1)
  } else {
    pi <- par
  }

  # Catch all zero:
  if (all.equal(pi, rep(0, D))==T){tr(solve(scaled_XTX))
  } else{

  # Block diagonal pi.
  bPi <- blockLambda_fcn(pi, p = p, D = D)

  # Theta maker:
  A <- solve( t(kron) %*% bPi %*% scaled_XTX %*% kron) %*%
    t(kron) %*% bPi %*% scaled_XTX
  # bias_comp
  b_comp <- bPi %*% (kron %*% A - diag(1, p*D))

  var_mle <- tr(solve(scaled_XTX))

  diff_sq <- t(indv_mle) %*% t(b_comp) %*% b_comp %*% indv_mle

  cov_piece <- 2*tr(
    solve(scaled_XTX) %*% t(b_comp)
  )

  # Return the sum.
  var_mle +
    diff_sq +
    cov_piece
  }
}


### MSE Hat without C:
mse_hat <- function(par, scaled_XTX, indv_mle, D, p, kron, par_type = "pi", c_comp = NULL){
  if (par_type == "lambda") {
    pi <- par/(par+1)
  } else {
    pi <- par
  }

  # Catch all zero:
  if (all.equal(pi, rep(0, D))==T){tr(solve(scaled_XTX))
  } else{
    # Reparameterize to include scaling:
    c <- max(pi)
    pi2 <- pi/c

    # Block diagonal pi.
    bPi <- blockLambda_fcn(pi2, p = p, D = D)

    # Theta maker:
    A <- solve( t(kron) %*% bPi %*% scaled_XTX %*% kron) %*%
      t(kron) %*% bPi %*% scaled_XTX
    # bias_comp
    b_comp <- bPi %*% (kron %*% A - diag(1, p*D))

    var_mle <- tr(solve(scaled_XTX))

    diff_sq <- t(indv_mle) %*% t(b_comp) %*% b_comp %*% indv_mle

    var_betahat <- tr( t(b_comp) %*% b_comp %*% solve(scaled_XTX))

    cov_piece <- tr(
      solve(scaled_XTX) %*% t(b_comp)
    )

    naive_cstar <- -cov_piece /(diff_sq + var_betahat)
    hc_star <- -cov_piece/diff_sq

    if (is.null(c_comp)){c_comp <-  naive_cstar}# If c_comp is blank, use c_naive.

       c^2*diff_sq -
       2*c*diff_sq*(naive_cstar)
  }
}




# Naive MSE (BMSE):
naive_mse <- function(par, scaled_XTX, indv_mle, D, p, kron, par_type = "pi"){
  if (par_type == "lambda") {
    pi <- par/(par+1)
  } else {
    pi <- par
  }

  # Block diagonal pi.
  bPi <- blockLambda_fcn(pi, p = p, D = D)

  # Theta maker:
  A <- solve( t(kron) %*% bPi %*% scaled_XTX %*% kron) %*%
    t(kron) %*% bPi %*% scaled_XTX

  var_mle <- tr(solve(scaled_XTX))

  diff_sq <- t(indv_mle - kron %*% A %*% indv_mle) %*%
    bPi %*% bPi %*%
    (indv_mle - kron %*% A %*% indv_mle)

  cov_piece <- 2*tr(
    solve(scaled_XTX) %*% (- bPi + bPi %*% kron %*% A
    )
  )

  extra_piece1 <- bPi %*%  (kron %*%A - diag(1, p*D))
  extra_piece <- tr(extra_piece1 %*% solve(scaled_XTX) %*% t(extra_piece1))

  # Return the sum.
  var_mle + diff_sq + cov_piece + extra_piece

}






########### Dataset-wise Functions ##########

studywise_mse <- function(pi, scaled_XTX, scaled_XTX_list, beta, beta_mat, p, D, kron){
  # Block diagonal pi.
  bPi <- blockLambda_fcn(pi, p = p, D = D)

  if (all.equal(pi, rep(0, D))==T){
    sapply(1:D, function(l){
      tr(solve(scaled_XTX_list[[l]]))
    })
  } else{
    # Theta maker:
    A <- solve( t(kron) %*% bPi %*% scaled_XTX %*% kron) %*%
      t(kron) %*% bPi %*% scaled_XTX

    # Returns a vector of MSEs.
    sapply(1:D, function(j){
      S_j <- kronecker(selBlock(j, D), diag(1, p))
      pi_j <- pi[j]
      scaled_XTX_j <- scaled_XTX_list[[j]]
      beta_j <- beta_mat[, j]
      var_comp <- (1-pi_j)^2 * solve(scaled_XTX_j) +
        pi_j^2 * A %*% solve(scaled_XTX) %*% t(A) +
        2*(1-pi_j)*pi_j* t(kron) %*% S_j %*% solve(scaled_XTX) %*% t(A)
      bias_comp <- pi_j*(A %*% beta - beta_j)

      # Return the sum of trace of var and sq bias.
      tr(var_comp) + t(bias_comp) %*% bias_comp
    })
  }
}




###### cstar ###############
c_star <- function(const_bpi, scaled_XTX, p, D, beta, kron){

  # Turn constrained pi vec into block diag matrix;
  bpi <- blockLambda_fcn(const_bpi, p = p, D = D)
  K <- kron
  inv_scaled_XTX <- solve(scaled_XTX)

  # big I, pk times pk;
  big_I  <- diag(1, nrow= p*D, ncol = p*D)

  # Theta Maker;
  A <- solve( t(K) %*% bpi %*% scaled_XTX %*% K ) %*% t(K) %*% bpi %*% scaled_XTX
  B <- K %*% A - big_I

  # trCov;
  trCov <- tr( inv_scaled_XTX %*% t(B) %*% bpi)

  # trVar;
  trVar <- tr(bpi %*% B %*% inv_scaled_XTX %*% t(B) %*% bpi)

  # sqBias;
  sqBias <- t(beta) %*% t(B) %*% t(bpi) %*% bpi %*% B %*% beta

  min( -trCov / (trVar + sqBias), 1 )
}


########## Functions for Ridgelike #######################

## Create the contrast matrix ####
contr_matrix <- function(D, p){
  choices <- c(1, -1)

  # Then, every combination of 1 - d
  main <- combn(c(1:D), 2)

  # Create a d by ncol(main) matrix of zeros.
  mat <- matrix( rep(0, D*ncol(main)), nrow = D, ncol = ncol(main))

  # Then, the first row of main gives us the position of 1, and
  # the second gives position of negative 1.
  ones <- main[1,]
  neg <- main[2,]
  for (j in 1:ncol(main)){
    mat[ones[j], j] <- 1
    mat[neg[j], j] <- -1
  }
  kronecker(mat %*% t(mat), diag(1, p))
}


### Ridge-like Estimator #####
ridgelike <- function(lambda, scaled_XTX, scaled_XTY,
                        p, D){
  cont <- contr_matrix(D, p)
  solve(scaled_XTX + 2*lambda*cont) %*% scaled_XTY
}


rl_mse <- function(lambda, scaled_XTX, p, D, beta){
  cont <- contr_matrix(D, p)
  interior <- solve(scaled_XTX + 2*lambda*cont)
  var_comp <- tr(interior %*% scaled_XTX %*% interior)

  bias <- (interior %*% scaled_XTX - diag(1, p*D)) %*% beta
  bias_comp <- t(bias) %*% bias
  var_comp + bias_comp
}

rl_umse <- function(lambda, scaled_XTX, p, D, indv_mle){
  cont <- contr_matrix(D, p)
  interior <- solve(scaled_XTX + 2*lambda*cont)
  R <- interior %*% scaled_XTX

  bias <- (R - diag(1, nrow = p*D)) %*% indv_mle
  bias_sq <- t(bias) %*% bias

  var_comp <- tr((2*R - diag(1, nrow= p*D))*solve(scaled_XTX))
  bias_sq + var_comp
}


rl_aCI <- function(lambda, scaled_XTX, p, D, indv_mle, alpha = .05){
  cont <- contr_matrix(D, p)
  interior <- solve(scaled_XTX + 2*lambda*cont)
  R <- interior %*% scaled_XTX

  beta_hat <- R %*% indv_mle
  var <- R %*% solve(scaled_XTX) %*% t(R)

  z_crit <- qnorm(1-alpha/2, mean = 0, sd = 1)
  lb <- beta_hat - z_crit * sqrt(diag(var))
  ub <- beta_hat + z_crit * sqrt(diag(var))
  list(lb = lb, ub = ub)
}


################ KLD linear #############################################################


## Alternative expression for linear KLD:
kld_linear_alt <- function(pi, scaled_XTY = NULL, scaled_XTX, indv_mle, p, D, kron){

  # Make block diagonal pi;
  bPi <- blockLambda_fcn(pi, p, D)



  theta_hat <- solve(t(kron) %*% bPi %*%
          scaled_XTX %*% kron) %*%
    t(kron) %*% bPi %*%
    scaled_XTX %*% indv_mle

  beta_hat <- (diag(1, nrow = p*D)-bPi) %*% indv_mle +
    bPi %*% kron %*%
    theta_hat


  list(beta_hat = beta_hat, theta_hat = theta_hat)
}


kld_theta <- function(pi, scaled_XTX, beta, p, D, kron){

  # Make block diagonal pi;
  bPi <- blockLambda_fcn(pi, p, D)


  theta <- solve(t(kron) %*% bPi %*%
                       scaled_XTX %*% kron) %*%
    t(kron) %*% bPi %*%
    scaled_XTX %*% beta

  theta
}



## Asymptotic confidence intervals:
kld_linear_aCI <- function(pi, scaled_XTX, indv_mle, p, D, kron, alpha){
  # Make block diagonal pi;
  if (all.equal(pi, rep(0, D))==T){var <- solve(scaled_XTX)
  } else{
  bPi <- blockLambda_fcn(pi, p, D)
  I <- diag(1, nrow = p*D)
  A <- solve(t(kron) %*% bPi %*% scaled_XTX %*% kron) %*% t(kron) %*% bPi %*% scaled_XTX



  theta_hat <- A %*% indv_mle

  beta_hat <- (I-bPi) %*% indv_mle +
    bPi %*% kron %*%
    theta_hat

  beta_hat_prime <- I - bPi +  bPi %*% kron %*% A
  interior <- solve(scaled_XTX)
  var <- beta_hat_prime %*% interior %*% t(beta_hat_prime)
  }

  z_crit <- qnorm(1-alpha/2, mean = 0, sd = 1)
  lb <- beta_hat - z_crit * sqrt(diag(var))
  ub <- beta_hat + z_crit * sqrt(diag(var))

  list(lb = lb, ub = ub)
}



#### bbeta_m ci ######
ham_bbetam_aci <- function(pi, scaled_XTX, tau, indv_mle, p, D, kron, alpha){
  # Make block diagonal pi;
  bPi <- blockLambda_fcn(pi, p, D)



  theta_hat <- solve(t(kron) %*% bPi %*%
                       scaled_XTX %*% kron) %*%
    t(kron) %*% bPi %*% scaled_XTX %*% indv_mle

  beta_hat <- (diag(1, nrow = p*D)-bPi) %*% indv_mle +
    bPi %*% kron %*%
    theta_hat

  hbbeta_m <- 1/D*t(kron(beeta_hat))

  beta_hat_prime <- diag(1, nrow = p*D) - bPi +
    bPi %*% kron %*% solve(t(kron) %*% bPi %*%
                             scaled_XTX %*% kron) %*%
    t(kron) %*% bPi %*% scaled_XTX
  interior <- solve(scaled_XTX)
  var <- beta_hat_prime %*% interior %*% t(beta_hat_prime)

  z_crit <- qnorm(1-alpha/2, mean = 0, sd = 1)
  lb <- beta_hat - z_crit * sqrt(diag(var))
  ub <- beta_hat + z_crit * sqrt(diag(var))

  list(lb = lb, ub = ub)

}

########### Revised code from FEShR Package ############
# The three functions here are edited versions of functions from the FEShR R
# package, which implements Kwon (2026).
# Revise fe_shrink to use different bounds on the optimal mu.
# We change only code within "gen" option, as this is what we're using.
# Revise opt_mu_Lambda_URE, which minimizes URE wrt mu, given lambda,
# so that mu is bounded between min and max observed mles.
# Revise gen_obj, the objective function, to use opt_mu_Lambda_URE_revise.
fe_shrink_revise <- function (y, M, centering = c("0", "gen", "cov"), W = NULL,
          Z = list(), type = c("URE", "EBMLE"), diag_lam = FALSE,
          tau = 0.95, B = 1000, n_init_vals = 1, all_init_vals = FALSE,
          optim_control = list())
{
  if (!is.matrix(y)) {
    warning("y is not a matrix; converting to a matrix assuming T=1.")
    y <- matrix(y, nrow = 1)
  }
  if (!is.list(M)) {
    stop("M is not a list.")
  }
  if (ncol(y) != length(M)) {
    stop("ncol(y) and length(M) must match.")
  }
  centering <- match.arg(centering)
  type <- match.arg(type)
  print(paste0("method=", type, ", centering=", centering))
  T <- nrow(y)
  J <- ncol(y)
  missing_cells <- any(is.na(y))
  if (missing_cells) {
    O <- vector("list", length = J)
    for (j in 1:J) {
      O[[j]] <- which(!is.na(y[, j]))
      y[-O[[j]], j] <- 0
    }
  }
  else {
    O <- NULL
  }
  if (is.null(W))
    W <- diag(T)
  if (centering == "0") {
    grad <- gen_deriv(y, M, type = type, W, missing_cells,
                      O, diag_lam)
    bounds <- NULL
    Z <- list()
  }
  else if (centering == "gen") {
    bound1 <- sapply(1:T, function(t) min(y[t, ]))
    bound2 <- sapply(1:T, function(t) -max(y[t, ]))
    bounds <- c(bound1, bound2)
    grad <- NULL
    Z <- list()
  }
  else if (centering == "cov") {
    bounds <- B * abs(lm(as.vector(y) ~ do.call("rbind",
                                                Z))$coef)
  }
  obj <- gen_obj_revise(y, M, type, W, missing_cells, O, diag_lam,
                 centering, bounds, Z)
  if (!diag_lam) {
    # Edit here to no borrowing start.
    init_val_mat <- matrix(rnorm(n_init_vals * T * (T + 1)/2),
                           nrow = n_init_vals)
  }
  else {
    # No edits made in diagonal lambda case, as we do not use this.
    init_val_mat <- matrix(rnorm(n_init_vals * T), nrow = n_init_vals)
  }
  opt <- NULL
  all_vals <- NULL
  all_pars <- NULL
  for (j in 1:n_init_vals) {
    init_val <- init_val_mat[j, ]
    opt_res <- optim(par = init_val, obj, method = "BFGS",
                     control = optim_control)
    if (all_init_vals) {
      all_vals <- c(all_vals, opt_res$val)
      all_pars <- rbind(all_pars, opt_res$par)
      if (j == n_init_vals) {
        print(paste0(rep("=", 80), collapse = ""))
        print(all_vals)
        print(all_pars)
      }
    }
    opt <- min(opt_res$val, opt)
    if (opt == opt_res$val) {
      if (!diag_lam) {
        L <- opt_res$par
        Lambda_opt <- make_from_lowertri(L, T)
      }
      else {
        D <- opt_res$par
        Lambda_opt <- diag(D)^2
      }
    }
  }
  if (centering == "0") {
    mu_opt <- matrix(0, T, J)
  }
  else {
    if (type == "URE") {
      mu_opt <- opt_mu_Lambda_URE_revise(Lambda_opt, y, M, bounds,
                                  W, missing_cells, O, Z)
    }
    else if (type == "EBMLE") {
      mu_opt <- opt_mu_Lambda_nll(Lambda_opt, y, M, missing_cells,
                                  O, Z)
    }
  }
  thetahat <- get_thetahat(mu_opt, Lambda_opt, y, M, missing_cells,
                           O)
  return(list(thetahat = thetahat, mu_opt = mu_opt, Lambda_opt = Lambda_opt,
              obj_val = opt, all_vals = all_vals, all_pars = all_pars))
}




# Revised objective function for selecting Lambda
gen_obj_revise<- function (y, M, type, W, missing_cells, O, diag_lam, centering,
                           bounds = NULL, Z = list())
{
  T <- nrow(y)
  J <- ncol(y)
  if (!diag_lam) {
    if (centering == "0") {
      if (type == "URE") {
        obj <- function(L) {
          Lambda <- make_from_lowertri(L, T)
          return(URE(mu = matrix(0, T, J), Lambda, y,
                     M, W, missing_cells, O))
        }
      }
      else if (type == "EBMLE") {
        obj <- function(L) {
          Lambda <- make_from_lowertri(L, T)
          return(nll(mu = matrix(0, T, J), Lambda, y,
                     M, missing_cells, O))
        }
      }
    }
    else if (centering == "gen" | centering == "cov") {
      if (type == "URE") {
        obj <- function(L) {
          # Turn initial value into a matrix
          Lambda <- make_from_lowertri(L, T)
          # With that matrix solve quadratic programming problem that returns the mu that minimizes URE.
          opt_mu <- opt_mu_Lambda_URE_revise(Lambda, y, M,
                                             bounds, W, missing_cells, O, Z)
          # Take that optimal mu and the lambda we tried, and then find the URE using that.
          return(URE(opt_mu, Lambda, y, M, W, missing_cells,
                     O))
        }
      }
      else if (type == "EBMLE") {
        obj <- function(L) {
          Lambda <- make_from_lowertri(L, T)
          opt_mu <- opt_mu_Lambda_nll(Lambda, y, M,
                                      missing_cells, O, Z)
          if (centering == "cov") {
            opt_mu <- gam_to_mu(opt_mu, Z, T)
          }
          return(nll(opt_mu, Lambda, y, M, missing_cells,
                     O))
        }
      }
    }
  }
  else {
    if (centering == "0") {
      if (type == "URE") {
        obj <- function(D) {
          Lambda <- diag(D)^2
          return(URE(mu = matrix(0, T, J), Lambda, y,
                     M, W, missing_cells, O))
        }
      }
      else if (type == "EBMLE") {
        obj <- function(D) {
          Lambda <- diag(D)^2
          return(nll(mu = matrix(0, T, J), Lambda, y,
                     M, missing_cells, O))
        }
      }
    }
    else if (centering == "gen" | centering == "cov") {
      if (type == "URE") {
        obj <- function(D) {
          Lambda <- diag(D)^2
          # No edit made here, as we do not use diagonal Lambda option in our
          # simulation or data analysis.
          opt_mu <- opt_mu_Lambda_URE(Lambda, y, M,
                                      bounds, W, missing_cells, O, Z)
          return(URE(opt_mu, Lambda, y, M, W, missing_cells,
                     O))
        }
      }
      else if (type == "EBMLE") {
        obj <- function(D) {
          Lambda <- diag(D)^2
          opt_mu <- opt_mu_Lambda_nll(Lambda, y, M,
                                      missing_cells, O, Z)
          return(nll(opt_mu, Lambda, y, M, missing_cells,
                     O))
        }
      }
    }
  }
  return(obj)
}

opt_mu_Lambda_URE_revise <- function (Lambda, y, M, bounds, W, missing_cells, O, Z = list())
{
  T <- nrow(y)
  J <- ncol(y)
  Dmat_dvec <- get_Dmat_dvec(Lambda, y, M, W, missing_cells,
                             O, Z)
  if (length(Z) == 0) {
    # Solve without constraints on mu.
    mu_opt <- quadprog::solve.QP(Dmat = Dmat_dvec$Dmat,
                                 dvec = as.numeric(Dmat_dvec$dvec),
                                 Amat = cbind(diag(T),
                                 -diag(T)), bvec = bounds)$solution
  }
  else {
    mu_opt <- solve(Dmat_dvec$Dmat, Dmat_dvec$dvec)
    mu_opt <- gam_to_mu(mu_opt, Z, T)
  }
  return(matrix(mu_opt, T, J))
}
