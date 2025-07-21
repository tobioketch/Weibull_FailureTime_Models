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
    real mu_shape;
    real var_shape;
    real mu_scale;
    real var_scale;
    real fctor;      // sensitivity of hyper - parameters
    real m;            // parameters of hyper-priors
    real<lower = 0> v;
}
transformed data{
    real sigma_shape_param = sqrt(var_shape)* fctor;
    real mu_shape_param = mu_shape * fctor;
    real<lower=0> mu_lnorm_alpha = (log(mu_scale / sqrt( (var_scale / square(mu_scale) ) +1 ) ) )* fctor;
    real<lower=0> sigma_lnorm_alpha = ( sqrt( log( (var_scale / square(mu_scale) ) +1 )  ) )* fctor;
}
parameters{
    real <lower=0> beta;          // shape 
    real <lower=0> alpha;        // scale 
}
model{
    mu_shape_param ~ normal(m, v);  //priors hyper parameters
    sigma_shape_param ~ lognormal(m, v);
    mu_lnorm_alpha ~ normal(m, v);
    sigma_lnorm_alpha ~ lognormal(m, v);
    beta ~ normal(mu_shape_param, sigma_shape_param);   //priors parameters
    alpha ~ lognormal(mu_lnorm_alpha, sigma_lnorm_alpha);
    target += weibull_lpdf(y_n | beta, alpha);  // likelihood
}
generated quantities{
  matrix[2,2] fisherInfo = weib2Par_fisherInfo(y_n, beta, alpha);
  matrix[2,2] varCov= weib2par_varCov(y_n, beta, alpha);
}

