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
    int n; // number of observed failure times
    vector[n] y_n;      // the obeserved failure times
    real mu_shape;    // prior mean shape param
    real var_shape;   // prior var shape param
    real mu_scale;    // prior mean scale param
    real fctor;      // sensitivity of hyper - parameters
    real m;            // parameters of hyper-priors
    real<lower = 0> v;
}
transformed data{
    real <lower=0> beta_gamma_shape = (square(mu_shape) / var_shape)* fctor;
    real <lower=0> beta_gamma_scale = (mu_shape / var_shape)* fctor ;
    real <lower=0> lambda_alpha = (inv(mu_scale))*fctor ;
}
parameters{
    real <lower=0> beta;
    real <lower=0>  alpha;
}
model{
    beta_gamma_shape ~ lognormal(m, v);  // hyper - parameter priors
    beta_gamma_scale ~ lognormal(m, v);
    lambda_alpha ~ lognormal(m, v);
    beta ~ gamma(beta_gamma_shape, beta_gamma_scale);  // parameter priors
    alpha ~ exponential(lambda_alpha);
    y_n ~ weibull(beta, alpha);   // likelihood observed failure time
}
generated quantities{
  matrix[2,2] fisherInfo = weib2Par_fisherInfo(y_n, beta, alpha);
  matrix[2,2] varCov= weib2par_varCov(y_n, beta, alpha);
}

