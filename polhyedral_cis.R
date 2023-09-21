# out truncation 
cdf_outer <- function(y, mu, sigma, a, b) {
  W = pnorm( (b - mu) / sigma) - pnorm( (a - mu) / sigma)  
  if (y <= a) {
    return(pnorm( (y - mu) / sigma) / (1 - W))
  }
  if (y > a & y <= b) {
    return(pnorm( (a - mu) / sigma) / (1 - W))
  }
  if (y > b) {
    return((pnorm( (y - mu) / sigma) - W) / (1 - W))
  }
}

# solve for CI 
quantile_finder <- function(y, p, sigma, a, b) {
  fu <- function(mu) { cdf_outer(y, mu, sigma, a, b) - (1 - p / 2)}
  fl <- function(mu) { cdf_outer(y, mu, sigma, a, b) - (p / 2)}
  int <- c(a - abs(y) + 2 * qnorm(p / 2), b + abs(y) - 2 * qnorm(p / 2))
  u   <- uniroot(fu, int)$root
  l   <- uniroot(fl, int)$root
  return(c(u, l))
} 

# examples 
x <- rnorm(1000000, 3, 2)
x <- x[x > 1.2 | x < -1.2]

mean(x < 2)
cdf_outer(y=2, mu=3, sigma=2, a=-1.2, b=1.2)


quantile_finder(y=2, p=0.05, sigma=2, a=-1.2, b=1.2)
quantile(x, c(0.025, 0.975))
