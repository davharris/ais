metropolis = function(x, fa, fb, beta, jump, ...){
  proposal = x + jump(...)
  
  old_p = exp(fa(x))^(1-beta) * exp(fb(x))^beta
  new_p = exp(fa(proposal))^(1-beta) * exp(fb(proposal))^beta
  
  if(new_p / old_p > runif(1)){
    return(proposal)
  }else{
    return(x)
  }
}
