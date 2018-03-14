load('final_dat.Robj')

loglikelihood <- function(data,theta)
{
  alpha1 <- theta[1]
  alpha2 <- theta[2]
  beta1 <- theta[3]
  beta2 <- theta[4]
  p.hat <- theta[5]
  return(sum(p.hat*(-1000*alpha1*log(beta1) - 1000*log(gamma(alpha1)) + (alpha1-1)*log(data) - data/beta1) + (1-p.hat)*(-1000*alpha2*log(beta2) - 1000*log(gamma(alpha2)) + (alpha2-1)*log(data) - data/beta2)))
}
em <- function(data,niter)
{
  # initialization
  kmean <- kmeans(data,2) # Initialization
  kmean$cluster[kmean$cluster==2] <- 0
  z1.hat <- kmean$cluster
  z2.hat <- 1-z1.hat
  p.hat <- sum(1-kmean$cluster)/length(data)
  beta1 <- var(data[kmean$cluster==1])/mean(data[kmean$cluster==1])
  beta2 <- var(data[kmean$cluster==0])/mean(data[kmean$cluster==0])
  alpha1 <- mean(data[kmean$cluster==1])/beta1
  alpha2 <- mean(data[kmean$cluster==0])/beta2
  theta0 <- c(alpha1=alpha1,alpha2=alpha2,beta1=beta1,beta2=beta2,p.hat=p.hat)

  for (i in 1:niter)
  {
    # E step - Compute z.hat, beta1, beta2, alpha1, alpha2

    alpha1 <- abs(-log(beta1) + mean(z1.hat*log(data)))
    alpha2 <- abs(-log(beta2) + mean(z2.hat*log(data)))
    beta1 <- abs(mean(data)/alpha1)
    beta2 <- abs(mean(data)/alpha2)
    
    
    # M step - Update z1
    z1.hat <- (p.hat*(1/(gamma(alpha1))*(beta1^alpha1))*(data^(alpha1-1))*exp(-data/beta1)) / ((p.hat*(1/(gamma(alpha1))*(beta1^alpha1))*(data^(alpha1-1))*exp(-data/beta1)) + ((1-p.hat)*(1/(gamma(alpha2))*(beta2^alpha2))*(data^(alpha2-1))*exp(-data/beta2)))
    z2.hat <- 1-z1.hat
    p.hat <- sum(1-z1.hat)/length(data)
    
    # Termination
    theta1 <- c(alpha1=alpha1,alpha2=alpha2,beta1=beta1,beta2=beta2,p.hat=p.hat)
    if ((sum(theta1-theta0)^2 < 1e-8) || (sum(loglikelihood(data,theta1) - loglikelihood(data,theta0))^2 < 1e-8))
    {
      print('YES!!')
      return(theta1)
    } 
    else
    {
      print('No...')
      theta0 <- theta1
    }
  }
}

em(dat,1000)

