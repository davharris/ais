# This example uses an "easy" t distribution with many degrees of freedom as
# a baseline for estimating the area under a "harder" t distribution with fewer
# degrees of freedom, in N dimensions.

# I compare the AIS approach with a standard importance sampler and usually find
# that AIS works better, depending on various details and how the random number
# generator behaves.

set.seed(1)       # Seed the random number generator for reproducibility

df = 2            # Degrees of freedom in the "hard" t distribution
dfa = 20          # Degrees of freedom in the "easy" t distribution
n = 10            # Number of dimensions of distribution
K = 10000         # Number of annealed transitions per run
replicates = 100  # Number of AIS runs


# Log-probability of "easy" distribution
fa = function(x){
  sum(dt(x, df = dfa, log = TRUE))
}

# Log-probability of "hard" distribution
fb = function(x){
  sum(dt(x, df = df, log = TRUE))
}

# Inverse temperatures
betas = cooling(K, exponent = -8)

# List of samples from the "easy" distribution
samples = replicate(replicates, rt(n, df = dfa), simplify = FALSE)

# Collect importance weights using annealed importance sampling
ais_weights = AIS(
  samples = samples, 
  betas = betas, 
  fa = fa, 
  fb = fb, 
  transition = metropolis, 
  jump = function(){rt(n, (df + dfa) / 2)}
)

# Collect importance weights using non-annealed importance sampling.
# Note that this sampler collects K times more samples than the AIS sampler,
# which should make them roughly similar in terms of total computational complexity.
# The samples from AIS tend to be higher quality, even though they're each more 
# difficult to compute.
s = matrix(rt(K * n * replicates, df = dfa), ncol = n)
is_weights = exp(rowSums(dt(s, df = df, log = TRUE) - dt(s, df = dfa, log = TRUE)))


# Results -----------------------------------------------------------------

mean(ais_weights)  # Should be close to 1
mean(is_weights)   # Should also be close to 1, but will usually be too low
########
sd(ais_weights) / sqrt(length(ais_weights)) # Estimated standard error (hopefully reliable)
sd(is_weights) / sqrt(length(is_weights))   # Estimated standard error (unreliable)
