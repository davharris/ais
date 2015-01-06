AIS = function(samples, betas, fa, fb, transition, parallel = TRUE, ...){
  if(parallel){
    f = mclapply
  }else{
    f = lapply
  }
  
  simplify2array(
    f(
      samples,
      function(x){
        run(x, betas = betas, fa = fa, fb = fb, transition = transition, ...)
      }
    )
  )
}

run = function(x, betas, fa, fb, transition, ...){
  K = length(betas)
  
  assert_that(all(betas <= 1))
  assert_that(all(betas >= 0))
  assert_that(all(betas == sort(betas)))
  assert_that(betas[K] == 1)
  
  # Empty vectors for storing each sample's negative energy under 
  # both parent distributions. Throughout, "a" refers to the prior/simple distribution
  # and "b" refers to the intractable one.
  f_as = numeric(K)
  f_bs = numeric(K)
  
  for(k in 1:K){
    # Sample at new temperature
    x = transition(x, fa, fb, betas[k], ...)
    
    # save negative energies under both distributions
    f_as[k] = fa(x)
    f_bs[k] = fb(x)
  }
  
  # Betas in numerator goes from 1:K
  # Betas in denominator go from 0:(K-1)
  w = exp(
    sum(
      log_pstar(f_as, f_bs, betas) - log_pstar(f_as, f_bs, c(0, betas[-K]))
    )
  )
  w
}
