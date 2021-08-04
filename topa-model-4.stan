data{
    int<lower=1> N;
    int<lower=1> N_site;
    int abund[N];
    int jsite[N];
    real xyear[N];
    //Maxn model
    int<lower=1> M; //number of obs in maxn
    int maxn[M];
    int xyear_m[M]; 
    int xyear_h[M];
    real nyears;
    real nyears_h;
}
parameters{
    real<lower=0> scale;
    real a;
    real b_n;
    real b_h;
    vector[N_site] aj;
    real<lower=0> sigma_aj;
    //Maxn model
    // real<lower=0> scale_m;
    real a_m;
    real b_m;
}
transformed parameters{
  real b_d;
  real b_dh;
  b_d = b_m + b_n;
  b_dh = (b_m * (nyears_h/nyears)) + b_h; 
}
model{
  vector[N] mu;
  vector[M] mu_m;

  // UVC model 
    scale ~ cauchy( 0 , 2 );
    sigma_aj ~ exponential( 1 );
    aj ~ normal( 0 , sigma_aj );
    b_n ~ normal(0 , 10);
    a ~ normal( 0 , 100);
    for ( i in 1:N ) {
        mu[i] = a + aj[jsite[i]] + b_n * xyear[i];
        mu[i] = exp(mu[i]);
    }
    abund ~ neg_binomial_2( mu , scale );
    
  //Maxn model 
    b_m ~ normal(0.0, 3);
    b_h ~ normal(0, 10);
    a_m ~ normal(0 , 100);
    for (j in 1:M){
      mu_m[j] = exp(a_m + b_d * xyear_m[j] + b_dh* xyear_h[j]);
      maxn[j] ~ poisson(mu_m[j]);
    }
  
}
generated quantities{
    vector[N] mu;
    vector[M] mu_m;
    real mu_hist; 
    real mu2000;
    real mu2018;
    real dev;
    dev = 0;

    for ( i in 1:N ) {
        mu[i] = a + aj[jsite[i]] + b_n * xyear[i];
        mu[i] = exp(mu[i]);
    }
    
    mu2000 = exp(a);
    mu2018 = exp(a + b_n);
    mu_hist = exp(a + b_h);
    
    for ( j in 1:M ) {
        mu_m[j] = a_m + b_d * xyear_m[j] + b_dh* xyear_h[j];
        mu_m[j] = exp(mu_m[j]);
    }
    dev = dev + (-2)*poisson_lpmf( maxn | mu_m );
    dev = dev + (-2)*neg_binomial_2_lpmf( abund | mu , scale );
}
