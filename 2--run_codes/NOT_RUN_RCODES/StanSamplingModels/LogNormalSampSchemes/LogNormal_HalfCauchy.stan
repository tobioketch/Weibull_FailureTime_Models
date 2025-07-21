functions{
  matrix weib2Par_fisherInfo(vector d, real shape, real scale){
    // returns a fisher information matrix of a 2parameter Weibull distribution
    int n = num_elements(d);
    matrix[2,2] Out;
    Out[1][1] = 1.823680 * (n / shape ^ 2);
    Out[2][2] = n * (shape / scale) ^ 2;
    Out[1][2] =  -0.422784 * n / scale;
    Out[2][1] = Out[1][2];
    return Out;
  }
  matrix weib2par_varCov (vector d, real shape, real scale){
    // returns variance-covariance matrix of a 2parameter Weibull distribution
    matrix[2,2] Out ;
    int n= num_elements(d);
    Out[2,2] = 1.1087*(scale/ shape)^2;
    Out[1,2] = 0.2570*scale;
    Out[2,1] =Out[1,2];
    Out[1,1]= 0.6079*shape^2;
    Out = Out/ n ;
    return Out;
  }
}
data{
    int n; 
    vector[n] y_n;
    real mu_shape;    // empirical/ bootstrapped mean shape param
    real var_shape;   // empirical var shape param
    real mu_scale;    // empirical mean scale param
    real var_scale;   // empirical var scale param
    real fctor;      // sensitivity of hyper - parameters
    real m;            // parameters of hyper-priors
    real<lower = 0> v;
}
transformed data{
    real mu_scale_param = mu_scale* fctor;
    real sigma_scale = sqrt(var_scale)* fctor;
    real mu_lnorm_beta = (log(mu_shape / sqrt( (var_shape / square(mu_shape) ) +1 ) ) )* fctor ;
    real<lower= 0> sigma_lnorm_beta = (sqrt( log( (var_shape / square(mu_shape) ) +1 )  ) )* fctor;
}
parameters{
    real <lower=0> beta;  // unknown parameters
    real <lower=0>  alpha;
}
model{
    mu_lnorm_beta ~ normal(m,v);   // hyper parameter priors   
    sigma_lnorm_beta ~ lognormal(m, v);
    mu_scale_param ~ normal(m,v);           // location cauchy distribution (real)
    sigma_scale ~ lognormal(m, v);        // scale Cauchy distribution (positive real)
    beta ~ lognormal(mu_lnorm_beta, sigma_lnorm_beta);  // parameter priors
    alpha ~ cauchy(mu_scale_param, sigma_scale);
    y_n ~ weibull(beta, alpha);  // likelihood
}
generated quantities{
  matrix[2,2] fisherInfo = weib2Par_fisherInfo(y_n, beta, alpha);
  matrix[2,2] varCov= weib2par_varCov(y_n, beta, alpha);
}


