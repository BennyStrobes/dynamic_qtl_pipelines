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
}
parameters {
  real<lower=0> nb_conc; 
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
    allele1 = exp(xb_1[n]);
    allele2 = exp(xb_2[n]);
    te = allele1 + allele2;
    if (gene_counts[n]>0) {
      gene_counts[n] ~ neg_binomial_2( te * library_size[n], nb_conc );
    }
  }
  nb_conc ~ gamma(concShape, concRate);
}