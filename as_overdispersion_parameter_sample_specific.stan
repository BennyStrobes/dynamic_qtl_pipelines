data {
  int<lower=0> N; # Number of samples
  int<lower=0> K; # Number of individuals
  int<lower=0> ys[N,1];
  int<lower=0> ns[N,1];
  int sample_names[N];  # sample indices
  real<lower=0> concShape; 
  real<lower=0> concRate;  

}
parameters {
  real<lower=0> conc[K]; 
}
model {
  for (n in 1:N) {
    if (ns[n,1]>0) {
      ys[n,1] ~ beta_binomial(ns[n,1], conc[sample_names[n]]*.5, conc[sample_names[n]]*(.5));
    }
  }
  conc ~ gamma(concShape, concRate);
}