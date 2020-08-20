####################################################################################
### File containing functions for the Bayesian GP modelling for prediction: with a fixed part
####################################################################################

##### Functions for kernel and posteriors (marginal and integrated out)

## Isomorphic kernel
K_iso_gaussian = function(phi) {
  ## single phi
  return (as.matrix(exp(-(D_all^2)/phi)))
}

K_all = K_iso_gaussian   # by default - Gaussian kernel

K = function(phi) {
  return (K_all(phi)[training_idx, training_idx])
}

S = function(phi) {
  return (K(phi) + diag(tau2, nrow = n_training))
}

##########################################################
## Posterior hyper-parameters (Yiolanda, pp 19-21, (2.20))
##########################################################

## F_ contains model terms (intercept and linear term). n x k matrix
F_ = function(n_training) {
  return (matrix(c(rep(1, n_training), linear_training), byrow = FALSE, ncol = 2))
}

# sigma2 ~ IG(alpha_, gamma_)
alpha_ = function (n_training) { return (alpha + 0.5*n_training) }

gamma_ = function(phi) {
  y_ = as.numeric(dt_energy[training_idx, energy]) - F_(n_training)%*%beta_0
  return (gamma + 
            0.5*(matrix(y_, nrow = 1)%*%solve(S(phi)+F_(n_training)%*%R%*%t(F_(n_training)))
                 %*%matrix(y_, ncol = 1)))
}

S_inv = function(phi) {
  return (solve(S(phi)))
}

Sigma_ = function(phi) {
  return (solve(t(F_(n_training))%*%S_inv(phi)%*%F_(n_training) + solve(R)))
}

mu_ = function(phi) {
  y = as.numeric(dt_energy[training_idx, energy])
  return (Sigma_(phi)%*%(t(F_(n_training))%*%S_inv(phi)%*%y + solve(R)%*%beta_0))
}

####################################
## Posterior distribution for phi
####################################

phi_p_Uniform = function(phi, svd.eps = 10^(-5)) {       ## using SVD, optimising on log scale
  A = S(phi) + F_(n_training)%*%R%*%t(F_(n_training))
  svd_A = svd(A)
  
  A_inv = svd_A$v%*%(diag(1./svd_A$d))%*%t(svd_A$u)
  y_ = as.numeric(dt_energy[training_idx, energy]) - F_(n_training)%*%beta_0
  gg = gamma + 0.5*(matrix(y_, nrow = 1)%*%A_inv%*%matrix(y_, ncol = 1))    
  
  return (-0.5*sum(log(svd_A$d)) - alpha_(n_training)*log(gg))
}

phi_p = phi_p_Uniform  # by default - Uniform prior

## Univariate phi optimisation

phi_iso_optim = function(phi_min = 10^(-5), phi_max = 10^10) {
  
  phi_opt = optimise(phi_p, lower = phi_min, upper = phi_max, maximum = TRUE)
  return (list(phi_post = phi_opt$maximum, phi_value = phi_opt$objective))
}

## Optimisation of phi.
phi_optim = function(phi_low = 10^(-10), n_start = 300) {
  m_phi_start = matrix(runif(n_start, min = 10^(-5), max = 10^3), ncol = 1) 
  phi_dim = ncol(m_phi_start)
  m_phi_post = m_phi_value = err_start = NULL
  for (i in 1:n_start) { 
    phi_start = m_phi_start[i,]
    ## Multivariate optimisation
    tryCatch(
      {
        phi_opt = optim(par = phi_start, fn = phi_p, control=list(fnscale=-1), lower = rep(phi_low, phi_dim),
                        method = "L-BFGS-B")  # SANN BFGS
        m_phi_post = rbind(m_phi_post, phi_opt$par)
        m_phi_value = rbind(m_phi_value, phi_opt$value)
      },
      silent = TRUE,
      error = function(e) {err_start = rbind(err_start, phi_start)}
    )
  }
  return (list(m_phi_post = m_phi_post, m_phi_value = m_phi_value, err_start = err_start))
}
## Posterior mode
phi_post_value = function(phi_opt) { return (phi_opt$m_phi_post[which.max(phi_opt$m_phi_value),]) }

########################################
## Predictive posterior distribution
########################################

## GP kernel between the test and training sets
k = function(phi) {
  return (K_all(phi)[test_idx, training_idx])
}

## y_pred ~ t_(2alpha) (1,a_ + b_y*mu_, gamma_/alpha_ * Sigma_y)
a_ = function(phi) {
  return (k(phi)%*%S_inv(phi))
}

b_t = function(phi) {
  return (matrix(c(rep(1, n_test), linear_test), byrow = FALSE, ncol=2) - a_(phi)%*%F_(n_training))
}

## Posterior predictive mean
mu_y = function(phi) {
  y = as.numeric(dt_energy[training_idx, energy])
  return (a_(phi)%*%y + b_t(phi)%*%mu_(phi))
}

Sigma_y = function(phi) {
  return (diag(1 + tau2, n_test) - a_(phi)%*%t(k(phi)) +
            b_t(phi)%*%Sigma_(phi)%*%t(b_t(phi))
  )
}

y_pred_mean = function(phi_post) {
  return (mu_y(phi_post))
}

# Variance-covariance matrix of posterior predictive, for the lower level
y_pred_var_low = function(phi_post) {
  return ((Sigma_y(phi_post))*as.numeric(gamma_(phi_post))/alpha_(n_training))
}

y_pred_sd = function(phi_post) {
  return (sqrt(diag(Sigma_y(phi_post))*as.numeric(gamma_(phi_post))/alpha_(n_training)))
}

# Variance-covariance matrix of posterior predictive, for the higher level
y_pred_var_high = function(phi_post, low_variance) {
  return (((mu_(phi_post)[2])^2 + as.matrix(Sigma_(phi_post))[2,2]*as.numeric(gamma_(phi_post))/alpha_(n_training))*low_variance 
          + y_pred_var_low(phi_post))
}

# SDs for high level predictive variances 
y_pred_sd_high = function(y_pred_var) {
  return (as.numeric(sqrt(diag(y_pred_var))))
}

y_true = function(test_idx) {
  return (dt_energy[test_idx, energy])
}

## Sampling predicted responses 
n_sample = 10^3
y_pred_sample = function(n_sample = 10^3, phi_post) {
  mm = matrix(rt.scaled(n = n_sample*n_test, df = 2*alpha_(n_training), 
                        mean = as.vector(y_pred_mean(phi_post)), 
                        sd = y_pred_sd(phi_post)),
              nrow = n_test, byrow = FALSE)
  return (mm)
}

#############################
## Choosing the training set
#############################

training_set_random = function(n_training, C_scaled) {
  n_s = nrow(C_scaled)
  training_idx = sort(sample.int(n_s, size = n_training, replace = FALSE))
  return (training_idx)
}

## Maximin design: scaling the sets of coordinates
get_coord_scaling = function(C_all) {
  for (k in 1:length(C_all)) {
    max_norm = max(as.numeric(apply(C_all[[k]], 1, norm, type = "2")))
    C_all[[k]] = C_all[[k]]/max_norm
  }
  return (C_all)
}

## Maximin design: stacking all the coordinates 
full_coord_matrix = function(C_all_scaled) {
  return (as.matrix(do.call(cbind, C_all_scaled)))
}

## C_matrix has n_s rows. Returns sorted training_idx
training_set_maximin = function(n_training, C_scaled) {
  mm_design = maximin.cand.upd(n = n_training, Xcand = C_scaled, Tmax = nrow(C_scaled)) 
  #cat(mm_design$Treached)
  return (sort(as.numeric(mm_design$inds)))
}

test_set = function(training_idx) {
  return (sort(setdiff(1:n_structures, training_idx)))
}

y_train = function(dt_energy, training_idx) {
  return (as.numeric(dt_energy[training_idx, energy]))
}

##########################################
### Obtain distance matrices
##########################################

## Not taking into account the symmetry.
get_distance_matrices = function(C_all) {
  n_dim = length(C_all)
  D_all = vector("list", n_dim)
  for (k in 1:n_dim) {
    D_all[[k]] = as.matrix(dist(C_all[[k]]))
  }
  return (D_all)
}

## Taking into account the symmetry of the molecule. C_all_scaled(!)
sym_oxalic_distance_matrices = function(C_all) {
  n_dim = length(C_all)
  D_all = vector("list", n_dim)
  for (k in 1:n_dim) {
    ## Atom reordering in case of symmetry of the molecule
    # ne = ncol(C_all[[k]])/n_atoms
    # atom_reorder = c((2*ne+1):(4*ne), 1:(2*ne), (5*ne+1):(6*ne),
    #                   (4*ne+1):(5*ne), (7*ne+1):(8*ne), (6*ne+1):(7*ne))
    D = matrix(0, ncol = n_structures, nrow = n_structures)
    for (s1 in 1:n_structures) {
      for (s2 in 1:n_structures) {
        d1 = sqrt(sum((C_all[[k]][s1,] - C_all[[k]][s2,])^2))
        #d2 = sqrt(sum((C_all[[k]][s1,] - C_all[[k]][s2, atom_reorder])^2))
        D[s1,s2] = D[s2,s1] = d1 #min(d1,d2)
      }
    }
    D_all[[k]] = D
  }
  return (D_all)
}

##########################################
### Output functions: plotting, MAE, RMSE
##########################################

## Some plotting

### generate data frame with predicted, true values, SD of the predictions and min distance to the training set
gen_df_pred = function(phi_post) {
  df_pred = data.frame("y_true" = as.numeric(dt_energy[test_idx, energy]),
                       "y_pred" = y_pred_mean(phi_post),
                       "y_sd" = y_pred_sd(phi_post),
                       "min_dist" = apply(Dist_scaled[test_idx, training_idx], 1, min))
  return (df_pred)
}

# as above, but also with the indicator of whether it is a minimum energy structure
gen_df_pred_perturbed = function(phi_post) {
  df_pred = data.frame("y_true" = as.numeric(dt_energy[test_idx, energy]),
                       "y_pred" = y_pred_mean(phi_post),
                       "y_sd" = y_pred_sd(phi_post),
                       "min_dist" = apply(D_all[test_idx, training_idx], 1, min),
                       "minima" = c(rep(0, n_test - n_minima), rep(1, n_minima)))
}

## prediction means and true values
gen_pred_true_plot = function(df_pred) {
  gg = ggplot(data = df_pred, aes(x = y_true)) +
    geom_point(aes(y = y_pred, colour = "y_pred"), pch = 4) +
    geom_point(aes(y = y_true, colour = "y_true"), pch = 4) +
    ggtitle("Predicted responses") +
    ylab("Predicted and true responses") + xlab("True response") +
    theme_minimal()
  return (gg)
}

## prediction means, true values and training points

gen_pred_true_train_plot = function(df_pred, y_training) {
  gg = ggplot() +
    geom_point(aes(x = df_pred$y_true, y = df_pred$y_true, colour = "y_true"), pch = 4) +
    geom_point(aes(x = df_pred$y_true, y = df_pred$y_pred, colour = "y_pred"), pch = 4) +
    geom_point(aes(x = y_training, y = y_training), colour = "black", pch = 17) + 
    ggtitle("Predicted responses") +
    ylab("Predicted and true responses") +
    theme_minimal()
  return (gg)
}

## prediction means and prediciton variance (SD)
gen_pred_sd_plot = function(df_pred) {
  gg = ggplot(data = df_pred, aes(x = y_true)) +
    geom_line(aes(y = y_pred, colour = "y_pred")) +
    geom_line(aes(y = y_sd, colour = "y_sd")) +
    scale_y_continuous(sec.axis = sec_axis(~./1, name = "y_sd")) +
    ggtitle("Predicted responses and prediction SD") +
    theme_minimal()
  return (gg)
}

## prediciton variance (SD) and min distance to the training set
gen_sd_dist_plot = function(df_pred) {
  gg = ggplot(data = df_pred, aes(x = min_dist)) +
    geom_line(aes(y = y_sd, colour = "y_sd")) +
    ggtitle("Prediction SD and min distane to the training set") +
    theme_minimal()
  return (gg)
}

# Mean absolute error (MAE)
get_pred_mae = function(y_true, y_pred) {
  return (mean(abs(y_true - y_pred)))
}

# Root mean squared error (RMSE)
get_pred_rmse = function(y_true, y_pred) {
  return (sqrt(mean((y_true - y_pred)^2)))
}