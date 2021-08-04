data{
    int<lower=1> N;
    int<lower = 1> N_survey;
    int<lower=1> N_site;
    real lnsize[N];
    int jsite[N];
    int jsurvey[N];
    real xyear[N];
}
parameters{
    real a;
    real b;
    vector[N_site] asite;
    vector[N_survey] asurvey;
    real<lower=0> sigma_asite;
    real<lower=0> sigma_asurvey;
    real<lower=0> sigma;
}

model{
  vector[N] mu;

  // UVC model 
    sigma ~ exponential( 1 );
    sigma_asite ~ exponential( 1 );
    sigma_asurvey ~ exponential( 1 );
    asite ~ normal( 0 , sigma_asite );
    asurvey ~ normal( 0 , sigma_asurvey );
    b ~ normal(0 , 10);
    a ~ normal( 0 , 10);
    for ( i in 1:N ) {
        mu[i] = a + asurvey[jsurvey[i]] + asite[jsite[i]] + b * xyear[i];
    }
    lnsize ~ normal( mu, sigma );
}

generated quantities{
    real mu2000;
    real mu2018;
    real size2000;
    real size2018;
    real mupred[N];
    
     for ( i in 1:N ) {
        mupred[i] = a + asurvey[jsurvey[i]] + asite[jsite[i]] + b * xyear[i];
    }
    
    mu2000 = exp(a)*exp(sigma);
    mu2018 = exp(a + b)*exp(sigma);
    size2000 = exp(normal_rng(a, sigma));
    size2018 = exp(normal_rng(a+b, sigma));
    
}

