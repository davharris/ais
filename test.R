df = 2
dfa = 20
n = 5
K = 1000      # Number of transitions
replicates = 100


fa = function(x){
  sum(dt(x, df = dfa, log = TRUE))
}

fb = function(x){
  sum(dt(x, df = df, log = TRUE))
}

betas = cooling(K, exponent = -8)

samples = replicate(replicates, rt(n, df = dfa), simplify = FALSE)


ais_weights = AIS(
  samples = samples, 
  betas = betas, 
  fa = fa, 
  fb = fb, 
  transition = metropolis, 
  jump = function(){rnorm(n)},
  parallel = FALSE
)

s = matrix(rt(K * n * replicates, df = dfa), ncol = n)
is_weights = exp(rowSums(dt(s, df = df, log = TRUE) - dt(s, df = dfa, log = TRUE)))

mean(ais_weights)  # Should be close to 1
mean(is_weights)   # Should be close to 1, but will usually be too low
########
sd(ais_weights) / sqrt(length(ais_weights))
sd(is_weights) / sqrt(length(is_weights))
