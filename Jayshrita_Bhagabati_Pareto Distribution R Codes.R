library(VGAM)
library(optimx) 
library(plotly)
library(MASS)

# Question-3

real_data <- read.csv("D:/OneDrive - Indian Institute of Management/statistics-1-term-paper/city-sizes-us.csv")
data_vect <- real_data$Value[!is.na(real_data$Value)]
summary(data_vect)
hist(data_vect, xlab = "population size", main = "Histogram of US city population sizes")
plot(density(data_vect),  xlab = "population size", main = "Histogram of US city population sizes")

# Question-6 Pareto MLE computation function

mle_pareto <- function(x){
  sigma_hat <- min(x)
  lambda_hat <- length(x)/(sum(log(x/sigma_hat)))
  return(
    c(
      sigma_hat,
      lambda_hat
    )
  )
}

# Check if MLE function works

samp_mle_check <- rparetoI(10000, scale = 5, shape = 10)
mle_pareto(samp_mle_check)

# Question-3 Pareto fit for US city population sizes data

real_data_fit <- mle_pareto(data_vect)

# Checking goodness of fit and qq-plot

plot(density(data_vect),  xlab = "population size", main = "Histogram of US city population sizes", xlim = c(0, 9*10^6), ylim = c(0, 7*10^(-6)))
par(new = TRUE)
curve(dpareto(x, scale = real_data_fit[1], shape = real_data_fit[2]), xlim = c(0, 9*10^6), ylim = c(0, 7*10^(-6)), col = "blue", ylab = "", xlab = "")

qqplot(x = qpareto(ppoints(length(data_vect)), scale = real_data_fit[1], shape = real_data_fit[2]), y = data_vect)
abline(0, 1)

ks.test(x = data_vect, y = "ppareto", scale = real_data_fit[1], shape = real_data_fit[2]) 

# Question-4 Pareto contamination with exponential distribtution

r_pareto_exp <- function(size, par_scale, par_shape, exp_rate, doc){
  nr_vect <- c(
    rep(1, floor(doc*size)),
    rep(0, size - floor(doc*size))
  )
  r_vect <- sample(
    nr_vect,
    size = size,
    replace = TRUE
  )
  pareto_vect <- rpareto(n = size, scale = par_scale, shape = par_shape)
  exp_vect <- rexp(n = size, rate = exp_rate)
  mix_vect <- r_vect*exp_vect + (!r_vect)*pareto_vect
  return(
    mix_vect
  )
}

plot(density(r_pareto_exp(size = 1000, par_scale = 1, par_shape = 2, exp_rate = 1, doc = 0.5)), xlab = "x", main = "Contaminated distribution")

# d_pareto_exp <- function(x, par_scale, par_shape, exp_rate, doc){
#   dens <- doc*dexp(x, rate = exp_rate) + (1 - doc)*dpareto(x, scale = par_scale, shape = par_shape)
#   return(dens)
# }

samp_contaminated <- r_pareto_exp(
  size = 10000,
  par_scale = 1,
  par_shape = 2,
  exp_rate = 1,
  doc = 0.5
)

hist(samp_contaminated)

plot(density(samp_contaminated), main = "Contaminated Distribution", xlab = "x")

# curve(d_pareto_exp(
#   x,
#   par_scale = 1,
#   par_shape = 2,
#   exp_rate = 1,
#   doc = 0.5
# ), xlim = c(0, 10), ylab = "density")

# Pareto population functions for characteristics

pareto_mean <- function(scale, shape){
  (scale*shape)/(shape - 1)
}

pareto_variance <- function(scale, shape){
  (shape*scale^2)/(((shape - 1)^2)*(shape - 2))
}

pareto_inv_cdf <- function(scale, shape, prob){
  scale/((1-prob)^(1/shape))
}

pareto_iqr <- function(scale, shape){
  pareto_inv_cdf(scale, shape, prob = 0.75) - pareto_inv_cdf(scale, shape, prob = 0.25)
}

prop_gt_mean_sd <- function(scale, shape){
  1 - ppareto(pareto_mean(scale, shape) + sqrt(pareto_variance(scale, shape)), scale = scale, shape = shape)
}

pareto_skewness <- function(shape){
  ((2*(shape + 1))/(shape - 3))*sqrt((shape - 2)/shape)
}

pareto_kurtosis <- function(shape){
  (6*(shape^3 + shape^2 - 6*shape -2))/(shape*(shape - 3)*(shape - 4))
}

# Variation of characteristics with parameters

# Mean

plot_scale <- seq(0.1, 5, 0.1)
plot_shape <- seq(1.1, 5, 0.1)
plot_mean_data <- outer(plot_scale, plot_shape, pareto_mean)
plot_mean <- plot_ly(x = plot_scale, y = plot_shape, z = plot_mean_data) %>% 
  layout(
    title = "Mean",
    scene = list(
      xaxis = list(title = "scale"),
      yaxis = list(title = "shape"),
      zaxis = list(title = "mean")
    )
  ) %>% 
  add_surface()
plot_mean


# Variance 

plot_scale <- seq(0.1, 5, 0.1)
plot_shape <- seq(2.1, 5, 0.1)
plot_var_data <- outer(plot_scale, plot_shape, pareto_variance)
plot_var <- plot_ly(x = plot_scale, y = plot_shape, z = plot_var_data) %>% 
  layout(
    title = "Variance",
    scene = list(
      xaxis = list(title = "scale"),
      yaxis = list(title = "shape"),
      zaxis = list(title = "variance")
    )
  ) %>% 
  add_surface()
plot_var

# Median

plot_scale <- seq(0.1, 5, 0.1)
plot_shape <- seq(2.1, 5, 0.1)
plot_median_data <- outer(plot_scale, plot_shape, function(x, y){
  pareto_inv_cdf(x, y, 0.5)
})
plot_median <- plot_ly(x = plot_scale, y = plot_shape, z = plot_median_data) %>% 
  layout(
    title = "Median",
    scene = list(
      xaxis = list(title = "scale"),
      yaxis = list(title = "shape"),
      zaxis = list(title = "median")
    )
  ) %>% 
  add_surface()
plot_median

# IQR

plot_scale <- seq(0.1, 5, 0.1)
plot_shape <- seq(2.1, 5, 0.1)
plot_iqr_data <- outer(plot_scale, plot_shape, pareto_iqr)
plot_iqr <- plot_ly(x = plot_scale, y = plot_shape, z = plot_iqr_data) %>% 
  layout(
    title = "IQR",
    scene = list(
      xaxis = list(title = "scale"),
      yaxis = list(title = "shape"),
      zaxis = list(title = "IQR")
    )
  ) %>% 
  add_surface()
plot_iqr

# Prop gt mean + sd

plot_scale <- seq(0.1, 5, 0.1)
plot_shape <- seq(2.1, 5, 0.1)
plot_mean_sd_data <- outer(plot_scale, plot_shape, prop_gt_mean_sd)
plot_mean_sd <- plot_ly(x = plot_scale, y = plot_shape, z = plot_mean_sd_data) %>% 
  layout(
    title = "Proportion greater than mean + sd",
    scene = list(
      xaxis = list(title = "scale"),
      yaxis = list(title = "shape"),
      zaxis = list(title = "proportion")
    )
  ) %>% 
  add_surface()
plot_mean_sd

# Question - 5

sig <- 5
lam <- 10

# Sample mean

set.seed(123)

samp_sizes <- c(
  5, 
  10,
  15, 
  20,
  30,
  40,
  50,
  100,
  200,
  400,
  500,
  1000, 
  2000
)

samp_list <- list()

for(i in 1:length(samp_sizes)){
  samp_list[[i]] <- replicate(
    n = 5000,
    expr = rparetoI(
      n = samp_sizes[i],
      scale = sig,
      shape = lam
    )
  )
}


samp_means <- lapply(
  X = samp_list,
  FUN = function(x){
    mat <- x
    apply(
      X = mat,
      MARGIN = 2,
      FUN = mean
    )
  }
)

mean_ks_df <- data.frame(n = NA, ks_stat = NA, ks_p_val = NA, sw_stat = NA, sw_p_val = NA)
mean_gamma <- data.frame(n = NA, shape = NA, rate = NA)

iter_flag <- 0

lapply(
  X = samp_means,
  FUN = function(x){
    iter_flag <<- iter_flag + 1
    dat <- x
    hist(dat, freq = FALSE, xlab = "xbar", main = paste("n =", samp_sizes[iter_flag]))
    curve(
      dnorm(x, mean = mean(dat), sd = sd(dat)),
      add = TRUE, col = "red"
    )
    qqplot(
      x = qnorm(ppoints(length(dat)), mean = mean(dat), sd = sd(dat)),
      y  = dat,
      main = paste("n =", samp_sizes[iter_flag]),
      xlab = "theoretical",
      ylab = "sample"
    )
    qqline(y = dat, distribution =  function(p){qnorm(p, mean = mean(dat), sd = sd(dat))}, qtype = 5) 
    test_val_ks <- ks.test(x, "pnorm", mean = mean(dat), sd = sd(dat))
    test_val_sw <- shapiro.test(dat)
    mean_ks_df[iter_flag, ] <<- c(
      samp_sizes[iter_flag],
      test_val_ks$statistic,
      test_val_ks$p.value,
      test_val_sw$statistic,
      test_val_sw$p.value
    )
    if(iter_flag <= 7){
      temp <- fitdistr(dat, "gamma")
      mean_gamma[iter_flag, ] <<- c(
        samp_sizes[iter_flag],
        temp$estimate[1],
        temp$estimate[2]
      )
    }
  }
)

iter_flag <- 0

# write.csv(mean_ks_df, "samp_mean.csv")
# write.csv(mean_gamma, "mean_gamma.csv")

# Sample variance

set.seed(123)

samp_sizes <- c(
  5,
  10,
  15,
  20,
  30,
  40,
  50,
  100,
  200,
  500,
  1000,
  2000,
  15000
)

samp_list <- list()

for(i in 1:length(samp_sizes)){
  samp_list[[i]] <- replicate(
    n = 5000,
    expr = rparetoI(
      n = samp_sizes[i],
      scale = sig,
      shape = lam
    )
  )
}


samp_vars <- lapply(
  X = samp_list,
  FUN = function(x){
    mat <- x
    apply(
      X = mat,
      MARGIN = 2,
      FUN = var
    )
  }
)

var_ks_df <- data.frame(n = NA, ks_stat = NA, ks_p_val = NA, sw_stat = NA, sw_p_val = NA)
var_gamma <- data.frame(n = NA, shape = NA, rate = NA)

iter_flag <- 0

lapply(
  X = samp_vars,
  FUN = function(x){
    iter_flag <<- iter_flag + 1
    dat <- x
    hist(dat, freq = FALSE, xlab = "var", main = paste("n =", samp_sizes[iter_flag]))
    curve(
      dnorm(x, mean = mean(dat), sd = sd(dat)),
      add = TRUE, col = "red"
    )
    qqplot(
      x = qnorm(ppoints(length(dat)), mean = mean(dat), sd = sd(dat)),
      y  = dat,
      main = paste("n =", samp_sizes[iter_flag]),
      xlab = "theoretical",
      ylab = "sample"
    )
    qqline(y = dat, distribution =  function(p){qnorm(p, mean = mean(dat), sd = sd(dat))}, qtype = 5) 
    test_val_ks <- ks.test(x, "pnorm", mean = mean(dat), sd = sd(dat))
    test_val_sw <- shapiro.test(dat)
    var_ks_df[iter_flag, ] <<- c(
      samp_sizes[iter_flag],
      test_val_ks$statistic,
      test_val_ks$p.value,
      test_val_sw$statistic,
      test_val_sw$p.value
    )
    if(iter_flag <= length(samp_sizes)){
      temp <- fitdistr(dat, "gamma")
      var_gamma[iter_flag, ] <<- c(
        samp_sizes[iter_flag],
        temp$estimate[1],
        temp$estimate[2]
      )
    }
  }
)

iter_flag <- 0

# write.csv(var_ks_df, "samp_var.csv")
# write.csv(var_gamma, "var_gamma.csv")

# 95th percentile

# Pareto 95th percentile density

dperc19 <- function(x, sigma, lambda){
  380 * (x^(-1 - lambda)) * (lambda) * ((x/sigma)^(-lambda)) * (sigma^lambda) * ((1 - (sigma/x)^lambda)^18)
}

dperc95 <- function(x, sigma, lambda){
  7152314400 * x^(-1 - lambda) * lambda *  (x/sigma)^(-5 * lambda) * sigma^lambda * (1 - (sigma/x)^lambda)^94
}

set.seed(123)

samp_sizes <- c(
  5,
  10,
  15,
  20,
  30,
  40,
  50,
  100,
  200,
  500,
  1000,
  2000
)

samp_list <- list()

for(i in 1:length(samp_sizes)){
  samp_list[[i]] <- replicate(
    n = 5000,
    expr = rparetoI(
      n = samp_sizes[i],
      scale = sig,
      shape = lam
    )
  )
}


samp_perc <- lapply(
  X = samp_list,
  FUN = function(x){
    mat <- x
    apply(
      X = mat,
      MARGIN = 2,
      FUN = function(x){
        quantile(x, 0.95)
      }
    )
  }
)

perc_ks_df <- data.frame(n = NA, ks_stat = NA, ks_p_val = NA, sw_stat = NA, sw_p_val = NA)

iter_flag <- 0

lapply(
  X = samp_perc,
  FUN = function(x){
    iter_flag <<- iter_flag + 1
    dat <- x
    hist(dat, freq = FALSE, xlab = "95th percentile", main = paste("n =", samp_sizes[iter_flag]))
    curve(
      dnorm(x, mean = mean(dat), sd = sd(dat)),
      add = TRUE, col = "red"
    )
    if(iter_flag == 4){
      curve(
        dperc19(x, sig, lam),
        add = TRUE, col = "blue"
      )
    }
    
    if(iter_flag == 8){
      curve(
        dperc95(x, sig, lam),
        add = TRUE, col = "blue"
      )
    }
    
    qqplot(
      x = qnorm(ppoints(length(dat)), mean = mean(dat), sd = sd(dat)),
      y  = dat,
      main = paste("n =", samp_sizes[iter_flag]),
      xlab = "theoretical",
      ylab = "sample"
    )
    qqline(y = dat, distribution =  function(p){qnorm(p, mean = mean(dat), sd = sd(dat))}, qtype = 5) 
    test_val_ks <- ks.test(x, "pnorm", mean = mean(dat), sd = sd(dat))
    test_val_sw <- shapiro.test(dat)
    perc_ks_df[iter_flag, ] <<- c(
      samp_sizes[iter_flag],
      test_val_ks$statistic,
      test_val_ks$p.value,
      test_val_sw$statistic,
      test_val_sw$p.value
    )
  }
)

iter_flag <- 0

# write.csv(perc_ks_df, "samp_95perc.csv")

# IQR

set.seed(123)

samp_sizes <- c(
  5,
  10,
  15,
  20,
  30,
  40,
  50,
  100,
  200,
  300,
  500,
  1000,
  2000
)

samp_list <- list()

for(i in 1:length(samp_sizes)){
  samp_list[[i]] <- replicate(
    n = 5000,
    expr = rparetoI(
      n = samp_sizes[i],
      scale = sig,
      shape = lam
    )
  )
}


samp_iqr <- lapply(
  X = samp_list,
  FUN = function(x){
    mat <- x
    apply(
      X = mat,
      MARGIN = 2,
      FUN = function(x){
        quantile(x, 0.75) - quantile(x, 0.25)
      }
    )
  }
)

iqr_ks_df <- data.frame(n = NA, ks_stat = NA, ks_p_val = NA, sw_stat = NA, sw_p_val = NA)

iter_flag <- 0

lapply(
  X = samp_iqr,
  FUN = function(x){
    iter_flag <<- iter_flag + 1
    dat <- x
    hist(dat, freq = FALSE, xlab = "IQR", main = paste("n =", samp_sizes[iter_flag]))
    curve(
      dnorm(x, mean = mean(dat), sd = sd(dat)),
      add = TRUE, col = "red"
    )
    qqplot(
      x = qnorm(ppoints(length(dat)), mean = mean(dat), sd = sd(dat)),
      y  = dat,
      main = paste("n =", samp_sizes[iter_flag]),
      xlab = "theoretical",
      ylab = "sample"
    )
    qqline(y = dat, distribution =  function(p){qnorm(p, mean = mean(dat), sd = sd(dat))}, qtype = 5) 
    test_val_ks <- ks.test(x, "pnorm", mean = mean(dat), sd = sd(dat))
    test_val_sw <- shapiro.test(dat)
    iqr_ks_df[iter_flag, ] <<- c(
      samp_sizes[iter_flag],
      test_val_ks$statistic,
      test_val_ks$p.value,
      test_val_sw$statistic,
      test_val_sw$p.value
    )
  }
)

iter_flag <- 0

# write.csv(iqr_ks_df, "samp_iqr.csv")


# Maximum Likelihood Estimation

# sigma

set.seed(123)

samp_sizes <- c(
  5,
  10,
  15,
  20,
  30,
  40,
  50,
  100,
  200,
  300,
  500,
  1000,
  2000
)

samp_list <- list()

for(i in 1:length(samp_sizes)){
  samp_list[[i]] <- replicate(
    n = 5000,
    expr = rparetoI(
      n = samp_sizes[i],
      scale = sig,
      shape = lam
    )
  )
}


samp_sigma <- lapply(
  X = samp_list,
  FUN = function(x){
    mat <- x
    apply(
      X = mat,
      MARGIN = 2,
      FUN = function(x){
       temp <- mle_pareto(x)
       temp[1]
      }
    )
  }
)

sigma_ks_df <- data.frame(n = NA, ks_stat = NA, ks_p_val = NA, sw_stat = NA, sw_p_val = NA)
sigma_asym_var <- data.frame(n = NA, mse = NA, var = NA, bias = NA)

iter_flag <- 0

lapply(
  X = samp_sigma,
  FUN = function(x){
    iter_flag <<- iter_flag + 1
    dat <- x
    hist(dat, freq = FALSE, xlab = "sigma", main = paste("n =", samp_sizes[iter_flag]))
    curve(
      dnorm(x, mean = mean(dat), sd = sd(dat)),
      add = TRUE, col = "red"
    )
    qqplot(
      x = qnorm(ppoints(length(dat)), mean = mean(dat), sd = sd(dat)),
      y  = dat,
      main = paste("n =", samp_sizes[iter_flag]),
      xlab = "theoretical",
      ylab = "sample"
    )
    qqline(y = dat, distribution =  function(p){qnorm(p, mean = mean(dat), sd = sd(dat))}, qtype = 5) 
    test_val_ks <- ks.test(x, "pnorm", mean = mean(dat), sd = sd(dat))
    test_val_sw <- shapiro.test(dat)
    sigma_ks_df[iter_flag, ] <<- c(
      samp_sizes[iter_flag],
      test_val_ks$statistic,
      test_val_ks$p.value,
      test_val_sw$statistic,
      test_val_sw$p.value
    )
    sigma_asym_var[iter_flag, ] <<- c(
      samp_sizes[iter_flag],
      mean((dat - sig)^2),
      var(dat),
      mean(dat) - sig
    )
  }
)

iter_flag <- 0

# write.csv(sigma_ks_df, "samp_sigma_mle.csv")
# write.csv(sigma_asym_var, "samp_sigma_params.csv")


# lambda

set.seed(123)

samp_sizes <- c(
  5,
  10,
  15,
  20,
  30,
  40,
  50,
  100,
  200,
  300,
  500,
  1000,
  2000
)

samp_list <- list()

for(i in 1:length(samp_sizes)){
  samp_list[[i]] <- replicate(
    n = 5000,
    expr = rparetoI(
      n = samp_sizes[i],
      scale = sig,
      shape = lam
    )
  )
}


samp_lambda <- lapply(
  X = samp_list,
  FUN = function(x){
    mat <- x
    apply(
      X = mat,
      MARGIN = 2,
      FUN = function(x){
        temp <- mle_pareto(x)
        temp[2]
      }
    )
  }
)

lambda_ks_df <- data.frame(n = NA, ks_stat = NA, ks_p_val = NA, sw_stat = NA, sw_p_val = NA)
lambda_asym_var <- data.frame(n = NA, mse = NA, var = NA, bias = NA, asym_var = NA)

iter_flag <- 0

lapply(
  X = samp_lambda,
  FUN = function(x){
    iter_flag <<- iter_flag + 1
    dat <- x
    hist(dat, freq = FALSE, xlab = "lambda", main = paste("n =", samp_sizes[iter_flag]))
    curve(
      dnorm(x, mean = mean(dat), sd = sd(dat)),
      add = TRUE, col = "red"
    )
    qqplot(
      x = qnorm(ppoints(length(dat)), mean = mean(dat), sd = sd(dat)),
      y  = dat,
      main = paste("n =", samp_sizes[iter_flag]),
      xlab = "theoretical",
      ylab = "sample"
    )
    qqline(y = dat, distribution =  function(p){qnorm(p, mean = mean(dat), sd = sd(dat))}, qtype = 5) 
    test_val_ks <- ks.test(x, "pnorm", mean = mean(dat), sd = sd(dat))
    test_val_sw <- shapiro.test(dat)
    lambda_ks_df[iter_flag, ] <<- c(
      samp_sizes[iter_flag],
      test_val_ks$statistic,
      test_val_ks$p.value,
      test_val_sw$statistic,
      test_val_sw$p.value
    )
    lambda_asym_var[iter_flag, ] <<- c(
      samp_sizes[iter_flag],
      mean((dat - lam)^2),
      var(dat),
      mean(dat) - lam,
      (lam^2)/samp_sizes[iter_flag]
    )
  }
)

iter_flag <- 0

# write.csv(lambda_ks_df, "samp_lambda_mle.csv")
# write.csv(lambda_asym_var, "samp_lambda_params.csv")


# MLE of characteristics general

set.seed(123)

samp_sizes <- c(
  5,
  10,
  15,
  20,
  30,
  40,
  50,
  100,
  200,
  300,
  500,
  1000,
  2000
)

samp_list <- list()

for(i in 1:length(samp_sizes)){
  samp_list[[i]] <- replicate(
    n = 5000,
    expr = rparetoI(
      n = samp_sizes[i],
      scale = sig,
      shape = lam
    )
  )
}


samp_gen <- lapply(
  X = samp_list,
  FUN = function(x){
    mat <- x
    apply(
      X = mat,
      MARGIN = 2,
      FUN = function(x){
        temp <- mle_pareto(x)
        # pareto_mean(temp[1], temp[2])
        # pareto_inv_cdf(temp[1], temp[2], 0.5)
        # pareto_iqr(temp[1], temp[2])
        sqrt(pareto_variance(temp[1], temp[2]))
      }
    )
  }
)

gen_ks_df <- data.frame(n = NA, ks_stat = NA, ks_p_val = NA, sw_stat = NA, sw_p_val = NA)
gen_asym_var <- data.frame(n = NA, mse = NA, var = NA, bias = NA, asym_var = NA)
conf_int <- data.frame(n = NA, ll = NA, ul = NA)

iter_flag <- 0

lapply(
  X = samp_gen,
  FUN = function(x){
    iter_flag <<- iter_flag + 1
    dat <- x
    hist(dat, freq = FALSE, xlab = "lambda", main = paste("n =", samp_sizes[iter_flag]))
    curve(
      dnorm(x, mean = mean(dat), sd = sd(dat)),
      add = TRUE, col = "red"
    )
    qqplot(
      x = qnorm(ppoints(length(dat)), mean = mean(dat), sd = sd(dat)),
      y  = dat,
      main = paste("n =", samp_sizes[iter_flag]),
      xlab = "theoretical",
      ylab = "sample"
    )
    qqline(y = dat, distribution =  function(p){qnorm(p, mean = mean(dat), sd = sd(dat))}, qtype = 5) 
    test_val_ks <- ks.test(x, "pnorm", mean = mean(dat), sd = sd(dat))
    test_val_sw <- shapiro.test(dat)
    gen_ks_df[iter_flag, ] <<- c(
      samp_sizes[iter_flag],
      test_val_ks$statistic,
      test_val_ks$p.value,
      test_val_sw$statistic,
      test_val_sw$p.value
    )
    # gen_asym_var[iter_flag, ] <<- c(
    #   samp_sizes[iter_flag],
    #   mean((dat - paret)^2),
    #   var(dat),
    #   mean(dat) - sqrt(pareto_variance(sig, lam)),
    #   # sig^2 / (-1 + lam)^4 * (lam^2) / samp_sizes[iter_flag]
    #   # (400^(1/lam) * sig^2 * log(20)^2)/lam^4 * (lam^2) / samp_sizes[iter_flag]
    #   # (4^(1/lam) * sig^2 * log(2)^2)/lam^4 * (lam^2) / samp_sizes[iter_flag]
    #   # (sig^2 * ((4/3)^(1/lam) * log(4/3) - 4^(1/lam) * log(4))^2) / lam^4 * (lam^2) / samp_sizes[iter_flag]
    #   # ((1 + lam - lam^2)^2 * sig^2) / ((-2 + lam)^3 * (-1 + lam)^4 * lam) * (lam^2) / samp_sizes[iter_flag]
    # )
    conf_int[iter_flag, ] <<- c(
      samp_sizes[iter_flag], 
      mean(dat - 1.96*sd(dat)),
      mean(dat + 1.96*sd(dat))
    )
  }
)

iter_flag <- 0

# write.csv(lambda_ks_df, "samp_lambda_mle.csv")
# write.csv(lambda_asym_var, "samp_lambda_params.csv")



# MOM of characteristics general

mom_est <- function(x){
  a <- mean(x)
  b <- mean(x^2)
  c(
    -((-b + sqrt(b * (-a^2 + b)))/a),
    -((-a^2 + b + sqrt(b * (-a^2 + b)))/(a^2 - b))
  )
}


set.seed(123)

samp_sizes <- c(
  5,
  10,
  15,
  20,
  30,
  40,
  50,
  100,
  200,
  300,
  500,
  1000,
  2000
)

samp_list <- list()

for(i in 1:length(samp_sizes)){
  # samp_list[[i]] <- replicate(
  #   n = 5000,
  #   expr = rparetoI(
  #     n = samp_sizes[i],
  #     scale = sig,
  #     shape = lam
  #   )
  # )
  samp_list[[i]] <- replicate(
    n = 5000,
    expr = r_pareto_exp(
      size = samp_sizes[i],
      par_scale = sig,
      par_shape = lam,
      exp_rate = 5,
      doc = 0.75
    )
  )
}


samp_gen <- lapply(
  X = samp_list,
  FUN = function(x){
    mat <- x
    apply(
      X = mat,
      MARGIN = 2,
      FUN = function(x){
        # temp <- mom_est(x)
        # temp <- mle_pareto(x)
        # pareto_mean(temp[1], temp[2])
        # pareto_inv_cdf(temp[1], temp[2], 0.5)
        # pareto_iqr(temp[1], temp[2])
        # IQR(x)
        # sqrt(pareto_variance(temp[1], temp[2]))
        median(x)
      }
    )
  }
)

gen_ks_df <- data.frame(n = NA, ks_stat = NA, ks_p_val = NA, sw_stat = NA, sw_p_val = NA)
# gen_asym_var <- data.frame(n = NA, mse = NA, var = NA, bias = NA, asym_var = NA)
gen_asym_var <- data.frame(n = NA, mse = NA, var = NA, bias = NA)

iter_flag <- 0

lapply(
  X = samp_gen,
  FUN = function(x){
    iter_flag <<- iter_flag + 1
    dat <- x
    hist(dat, freq = FALSE, xlab = "lambda", main = paste("n =", samp_sizes[iter_flag]))
    curve(
      dnorm(x, mean = mean(dat), sd = sd(dat)),
      add = TRUE, col = "red"
    )
    qqplot(
      x = qnorm(ppoints(length(dat)), mean = mean(dat), sd = sd(dat)),
      y  = dat,
      main = paste("n =", samp_sizes[iter_flag]),
      xlab = "theoretical",
      ylab = "sample"
    )
    qqline(y = dat, distribution =  function(p){qnorm(p, mean = mean(dat), sd = sd(dat))}, qtype = 5) 
    test_val_ks <- ks.test(x, "pnorm", mean = mean(dat), sd = sd(dat))
    test_val_sw <- shapiro.test(dat)
    gen_ks_df[iter_flag, ] <<- c(
      samp_sizes[iter_flag],
      test_val_ks$statistic,
      test_val_ks$p.value,
      test_val_sw$statistic,
      test_val_sw$p.value
    )
    gen_asym_var[iter_flag, ] <<- c(
      samp_sizes[iter_flag],
      mean((dat - pareto_inv_cdf(sig, lam, 0.5))^2),
      var(dat),
      mean(dat) - pareto_inv_cdf(sig, lam, 0.5)
      # 1/(4*samp_sizes)
      # sig^2 / (-1 + lam)^4 * (lam^2) / samp_sizes[iter_flag]
      # (400^(1/lam) * sig^2 * log(20)^2)/lam^4 * (lam^2) / samp_sizes[iter_flag]
      # (4^(1/lam) * sig^2 * log(2)^2)/lam^4 * (lam^2) / samp_sizes[iter_flag]
      # (sig^2 * ((4/3)^(1/lam) * log(4/3) - 4^(1/lam) * log(4))^2) / lam^4 * (lam^2) / samp_sizes[iter_flag]
      # ((1 + lam - lam^2)^2 * sig^2) / ((-2 + lam)^3 * (-1 + lam)^4 * lam) * (lam^2) / samp_sizes[iter_flag]
    )
  }
)

iter_flag <- 0

# write.csv(lambda_ks_df, "samp_lambda_mle.csv")
# write.csv(lambda_asym_var, "samp_lambda_params.csv")

1/(4*c(200, 300, 500, 1000, 2000)*dparetoI(pareto_inv_cdf(sig,lam, 0.5), scale = sig, shape = lam)^2)


# How often does box plot detect outlier. 
set.seed(123)

samp_box <- replicate(
  n = 5000, 
  expr = rparetoI(2000, scale = sig, shape = lam)
)

prob_detect <- apply(
  X = samp_box,
  MARGIN = 2,
  FUN = function(x){
    sum(
      x >= quantile(x, 0.75) + 1.5*IQR(x) | x <= quantile(x, 0.25) - 1.5*IQR(x)
    )
  }
)

mean(prob_detect/2000)

set.seed(123)

samp_box <- replicate(
  n = 5000, 
  expr = rparetoI(1000, scale = sig, shape = 5)
)

prob_detect <- apply(
  X = samp_box,
  MARGIN = 2,
  FUN = function(x){
    sum(
      x >= quantile(x, 0.75) + 1.5*IQR(x) | x <= quantile(x, 0.25) - 1.5*IQR(x)
    )
  }
)

mean(prob_detect/1000)

