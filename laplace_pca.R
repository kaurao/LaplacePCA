laplace_pca <- function(data=NULL, e=NULL, d=NULL, n=NULL, alpha=0, beta=0, center=TRUE, scale=FALSE) {
  
  # Arguments
  # data: data matrix
  # e   : eigen values
  # d   : data dimensionality
  # n   : number of data points
  #       if data is not provided then e, d and n must be provided
  
  # Returns a list with elements
  # k   : estimated dimensionality
  # p   : selection criterion values for different dimensionalities
  # ks  : dimensionalities tested
  # i   : selected idex from ks, i.e. ks[i]==k
  # pmax: maximum of selection criterion
  
  # Transcribed from Tom Minka's MATLAB code
  # see: Automatic choice of dimensionality for PCA, Tom Minka
  # Note that center=TRUE is necessary to hold the theoretical results
  #
  # function [k,p,e] = laplace_pca(data, e, d, n, alpha, beta)
  # % LAPLACE_PCA   Estimate latent dimensionality by Laplace approximation.
  # %
  # % k = LAPLACE_PCA([],e,d,n) returns an estimate of the latent dimensionality
  # % of a dataset with eigenvalues e, original dimensionality d, and size n.
  # % LAPLACE_PCA(data) computes (e,d,n) from the matrix data 
  # % (data points are rows)
  # % [k,p] = LAPLACE_PCA(...) also returns the log-probability of each 
  # % dimensionality, starting at 1.  k is the argmax of p.
  # % Written by Tom Minka
  
  if(!is.null(data)) {
    n <- nrow(data)
    d <- ncol(data)
    e <- prcomp(data, center = center, scale. = scale)
    e <- e$sdev^2
  } else {
    stopifnot(!is.null(e))
    stopifnot(!is.null(d))
    stopifnot(!is.null(n))
  }
  
  # break off the eigenvalues which are identically zero
  i <- e < .Machine$double.eps
  e <- e[!i]
  
  logediff <- rep(0,length(e));
  for(i in 1:(length(e)-1)) {
    j <- (i+1):length(e)
    logediff[i] <- sum(log(e[i] - e[j])) + (d-length(e))*log(e[i])
  }
  cumsum_logediff <- cumsum(logediff)
  
  n1 <- n-1+alpha
  ehat <- (e + beta/n)*n/n1
  inve <- 1/ehat
  
  invediff <- replicate(length(e), inve) - t(replicate(length(e), inve))
  invediff[invediff <= 0] <- 1
  invediff <- log(invediff)
  
  cumsum_invediff <- apply(invediff,2,cumsum)
  row_invediff <- apply(cumsum_invediff,1,sum)
  
  loge <- log(ehat)
  cumsum_loge <- cumsum(loge)
  
  cumsum_e <- cumsum(ehat)
  
  dn <- length(e)
  kmax <- length(e)
  
  ks <- 1:kmax
  
  z <- log(2) + (d-ks+1)/2*log(pi) - lgamma((d-ks+1)/2)
  cumsum_z <- cumsum(z)
  p <- rep(NA, length(e))
  for(i in 1:length(ks)) {
    k <- ks[i]
    v <- (tail(cumsum_e,1) - cumsum_e[k])/(d-k)
    p[i] <- -n1/2*cumsum_loge[k] + (-n1*(d-k)/2)*log(v)
    p[i] <- p[i] - cumsum_z[k] - k/2*log(n1)
    
    h <- row_invediff[k] + cumsum_logediff[k]
    
    h <- h + (d-k)*sum(log(1/v - inve[1:k]));
    m <- d*k-k*(k+1)/2;
    h <- h + m*log(n);
    p[i] <- p[i] + (m+k)/2*log(2*pi) - h/2;
    
    p[i] <- p[i] + 1.5*k*log(2);
    p[i] <- p[i] - 0.5*log(d-k);
    if(alpha > 0) {
      ck <- alpha*(d-k)/2*log(beta*(d-k)/2) - lgamma(alpha*(d-k)/2) + 
        k*(alpha/2*log(beta/2) - lgamma(alpha/2));
      p[i] <- p[i] - n1*d/2 - 0.5*log(n1) + ck;
    }
  }
  
  pmax <- max(p, na.rm = TRUE)
  i <- which(p==pmax)
  i <- i[1]
  k <- ks[i]
  
  v0 <- tail(cumsum_e,1)/length(cumsum_e);
  p0 <- -n1*d/2*log(v0) - 0.5*log(d);
  
  if(alpha > 0) { 
    # must inline ck and put at end
    p0 <- p0 - n1*d/2 - 0.5*log(n1) + alpha*d/2*log(beta*d/2) - lgamma(alpha*d/2);
  }
  if(p0 >= pmax) k <- 0
  p <- c(p0, p)
  ks <- c(0, ks)
    
  return(list(k=k, p=p, e=e, ks=ks, i=i, pmax=pmax))
}
