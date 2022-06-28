##### This script is to store the main functions that is needed by GLM #######

# calculate e(-k.yji)
calculate_expNdot <- function(k, Yji){
  power <- Yji %*% k 
  #head(power)
  expNdot <- exp(-1 * power)
  #head(expNdot)
  return(expNdot)
}

# calculate UBj
calculate_UBj <- function(expNdot, gene_order){
  UBj <- expNdot %>% 
    Matrix.utils::aggregate.Matrix(., groupings = gene_order, fun = 'sum')
  #head(UBj)
  
  return(UBj)
}

# calculate alphaj
calculate_alphaj <- function(lambda, SBj, UBj){
  alphaj <- SBj$score / (lambda * UBj)
  #head(alphaj)
  
  return(alphaj)
}

# calculate VBj 
calculate_VBj <- function(expNdot, Yji, gene_order){
  VBj <- (Yji * as.vector(expNdot)) %>% 
    Matrix.utils::aggregate.Matrix(., groupings = gene_order, fun = 'sum') 
  #head(VBj)
  
  return(VBj)
}

# calculate simplified likelihood
# calculate_likelihood <- function(SBj, k, TBj, UBj){
#   item1 <- (-1)*SBj$score*log(UBj)
#   #head(item1)
#   
#   item2 <- TBj %*% k 
#   #head(TBj)
#   #head(item2)
#   
#   likelihood <- sum(item1-item2)
#   
#   return(likelihood)
# }

# calculate penalized likelihood, USE FULL LIKELIHOOD EQUATION
calculate_likelihood <- function(SBj, k, TBj, UBj, lambda1, lambda2, n){
  item1 <- SBj$score*(log(SBj$score) - log(UBj))
  item2 <- TBj %*% k
  
  penalty <- n*(lambda1 * (lambda2 * sum(abs(k)) + (1-lambda2)/2 * sum(k^2))) # takes number of bins, n
  
  likelihood <- sum(item1-item2-SBj$score)
  p_likelihood <- likelihood - penalty
  
  #print(likelihood/(penalty + 0.01)) # prevent when penalty could be 0
  
  return(p_likelihood)
}


# calculate gradient 
calculate_gradient <- function(lambda, alphaj, VBj, TBj, lambda1, lambda2, k, n){
  #head(VBj)
  #head(alphaj)
  item1 <- as.vector(lambda * alphaj) * VBj
  #head(item1)
  #head(TBj) 
  
  ## change gradient due to penalized likelihood
  penalty_g <- function(nth_k, lambda1, lambda2, n){
    if(nth_k == 0){
      p_g = 0 # treat derivative(|k1|, k1=0) = 0
    }else{
      p_g = nth_k/abs(nth_k) * lambda1 * lambda2 + lambda1 * (1 - lambda2) * nth_k
      p_g = p_g * n # take in the number of bins
    }
    return(p_g)
  }
  
  p_gradient <- sapply(k, penalty_g, lambda1, lambda2, n) 
  head(p_gradient) 
  
  gradient <- colSums(item1 - TBj) - p_gradient
  #head(gradient)
  
  return(gradient)
}
