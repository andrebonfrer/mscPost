# run_SC_wls.R

library(mscPost)
library(data.table)
library(Matrix)

rm(list=ls())
gc()


# read in SC object here

# source data which contains weights
a <- readRDS("~/Dropbox/Raiz/Results/ssc_AmountDeposit_l28_lambda_0.3.rds")
#a <- readRDS("~/Dropbox/Raiz/Results/ssc_signin_week_l28_nu0_lambda_0.3.rds")


#######


# add selection
flags <- list()
# for X
flags$cov.vars <- c("tvg.dummy","GR")
# for GR/selection
flags$cov.var.first.stage <- c(			"age",
                                  "gender_male",
                                  "gender_female",
                                  "income",
                                  "netwealth",
                                  "portfolio_2",
                                  "portfolio_3",
                                  "portfolio_4",
                                  "portfolio_5",
                                  "portfolio_6",
                                  "employment_Self.Employed",
                                  "employment_Full.Time",
                                  "employment_Student",
                                  "employment_Unemployed",
                                  "employment_Retired",
                                  "employment_NA",
                                  "reasonforinvest_Long_term",
                                  "reasonforinvest_Children",
                                  "reasonforinvest_Short_term",
                                  "reasonforinvest_Major_Purchase",
                                  "reasonforinvest_NA",
                                  "Users.Source.ID_referral",
                                  "Users.Source.ID_nonsocial_media",
                                  "Users.Source.ID_social_media",
                                  "timeframe_5",
                                  "timeframe_10",
                                  "timeframe_25",
                                  "timeframe_15",
                                  "timeframe_100"
)

# GR calculation instruments
flags$columns_instruments <- c("bank.recoded",
                               #  "timetoupdate.log",
                               "Devices.Platform",
                               "goal.set.count.by.zip",
                               #"deposit.count.by.zip",
                               "signin.count.by.zip")

dta <- a$data
#### any data work here
dta[,engagement := referrals_week+signin_week-inactive]

####### create any new conctrol variables for Z or X here
dta[,portfolioriskpreference := ifelse(portfolio_1,1,
                                       ifelse(portfolio_2,2,
                                              ifelse(portfolio_3,3,
                                                     ifelse(portfolio_4,4,
                                                            ifelse(portfolio_5,5,
                                                                   ifelse(portfolio_6,6,NA
                                                                   ))))))]
# prepare goal characteristics variables here
dta[, goal_difficulty := 100*winsorize(goal_difficulty)]
dta[, goal_commitment := 100*winsorize(goalcommitment2)]
dta[, initialgoalstance := 100*winsorize(initialgoalstance)]

# add back to data object
a$data <- dta

# Step 1: generate data object

gdata <- prepare_data(dta=a$data,
                      res=a$res,
                      f.X = "AmountDeposit ~ 1 + tvg.dummy",
                      f.Z = "~ 1 + goal_difficulty + goal_commitment + initialgoalstance +
      I(goal_difficulty^2) + I(goal_commitment^2) + I(initialgoalstance^2) ",
                      flags = flags
)

# Step 2: sample from posteriors of each set of parameters
# Perform Gibbs sampling
samples <- gibbs_sampling(gdata,
                          n_iter = 5,
                          burn_in = 2)

# Step 3: Extract samples and summarise
beta_samples <- samples$beta_samples
gamma_samples <- samples$gamma_samples
sigma2_samples <- samples$sigma2_samples
tau_samples <- samples$tau_samples

# Summary statistics
beta_mean_est <- apply(beta_samples, 2, mean)
beta_sd_est <- apply(beta_samples, 2, sd)

gamma_mean_est <- apply(gamma_samples, 2, mean)
gamma_sd_est <- apply(gamma_samples, 2, sd)

gamma_res <- cbind(gamma_mean_est,gamma_sd_est)
colnames(gamma_res) <- c("Mean", "Std Dev")

sigma2_mean_est <- mean(sigma2_samples)
tau_mean_est <- apply(tau_samples, 2, mean)

# Print summary statistics
print("Estimated beta:")
print(beta_mean_est)

print("Estimated gamma:")
print(gamma_mean_est)
print(gamma_sd_est)

print("Estimated sigma2:")
print(sigma2_mean_est)

print("Estimated tau:")
print(tau_mean_est)

#### try doing parallel
library(parallel)
library(coda)

# Function to run a single Gibbs sampling chain
run_single_chain <- function(seed) {
  set.seed(seed)
  samples <- gibbs_sampling(gdata = gdata,
                            n_iter = 2000, burn_in = 1000)
  return(samples)
}

# Number of parallel chains
num_chains <- 4

# Seeds for reproducibility
seeds <- 123 + 0:(num_chains - 1)

# Run the chains in parallel
chains <- mclapply(seeds, run_single_chain, mc.cores = num_chains)

# Extract the beta samples from each chain
beta_chains <- lapply(chains, function(chain) mcmc(chain$gamma_samples))

# Combine into an mcmc.list object
combined_mcmc <- mcmc.list(gamma_chains)

# Convergence diagnostics
gelman_diag <- gelman.diag(combined_mcmc)
print(gelman_diag)
