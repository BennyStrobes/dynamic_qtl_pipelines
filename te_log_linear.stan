data {
  int<lower=0> N; # Number of samples
  int<lower=0> P; # Number of covariates
  real library_size[N]; 
  matrix[N,P] x_1; 
  int<lower=0> gene_counts[N]; 
  real<lower=0> concShape; 
  real<lower=0> concRate;  
}
parameters {
  real<lower=0> nb_conc; 
  vector[P] beta;
}
model {
  vector[N] xb_1; 
  xb_1 <- x_1 * beta;
  for (n in 1:N) {
    real te; 
    te = exp(xb_1[n]);
    if (gene_counts[n]>0) {
      gene_counts[n] ~ neg_binomial_2( te * library_size[n], nb_conc );
    }
  }
  nb_conc ~ gamma(concShape, concRate);
}