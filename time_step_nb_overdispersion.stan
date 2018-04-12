data {
  int<lower=0> N; # num individuals
  int<lower=0> P; # num genes
  int<lower=0> T; # num time steps
  real library_size[N]; 
  int time_step[N];
  int<lower=0> gene_counts[N,P]; 
  real<lower=0> concShape; 
  real<lower=0> concRate;  
}
parameters {
  real<lower=0> conc[P,T];
  real<lower=0> gene_means[P];
}
model {
  for (n in 1:N) {
    for (p in 1:P) {
      gene_counts[n,p] ~ neg_binomial_2( gene_means[p] * library_size[n], conc[p, time_step[n]] ); 
    }
  }
  for (p in 1:P) {
    for (t in 1:T) {
      conc[p,t] ~ gamma(concShape, concRate);
    }
  }
}