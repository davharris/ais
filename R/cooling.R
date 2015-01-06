cooling = function(K, exponent = -8){
  assert_that(exponent < 0)
  assert_that(K == as.integer(K))
  assert_that(K > 0)
  
  out = 1 - .95^(1:K / K * (exponent * log(10) / log(.95)))
  out[K] = 1
  
  out
}