data {
  int<lower=0> N; # Number of samples
  int<lower=0> P; # Number of covariates
  int<lower=0> K; # Number of heterozygous, exonic snps
  real library_size[N]; 
  matrix[N,P] x_1; 
  matrix[N,P] x_2; 
  int<lower=0> ys[N,K];
  int<lower=0> ns[N,K];
  int<lower=0> gene_counts[N]; 
  real<lower=0> concShape; 
  real<lower=0> concRate; 
  real<lower=0> as_overdispersion_parameter; 
}
parameters {
  vector[P] beta;
}
model {
  vector[N] xb_1; 
  vector[N] xb_2; 
  xb_1 <- x_1 * beta;
  xb_2 <- x_2 * beta;
  for (n in 1:N) {
    real te; 
    real allele1; 
    real allele2; 
    real p; 
    allele1 = exp(xb_1[n]);
    allele2 = exp(xb_2[n]);
    te = allele1 + allele2;
    p = allele1/te;
    for (k in 1:K) {
      if (ns[n,k]>0) {
        ys[n,k] ~ beta_binomial(ns[n,k], as_overdispersion_parameter*p, as_overdispersion_parameter*(1.0-p));
      }
    }
  }
}