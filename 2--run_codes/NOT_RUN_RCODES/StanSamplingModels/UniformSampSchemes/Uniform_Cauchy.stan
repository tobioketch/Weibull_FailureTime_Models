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
    int  p; 
    vector[n] y_n;
    real  mu_shape;
    real  var_shape;
    real  mu_scale;
    real  var_scale;
}
transformed data{
    real lower_bound = mu_shape - sqrt(var_shape)*6;
    real upper_bound = mu_shape + sqrt(var_shape)*6;
}
parameters{
    real <lower=0> beta;
    real <lower=0>  alpha;
    real <lower=0> y_p[p];
}
model{
    // priors
    mu_scale ~ cauchy(0,1);
    var_scale ~ cauchy(0,1);
    beta ~ uniform((lower_bound < 0 ? 0 : lower_bound), upper_bound);
    alpha ~ cauchy(mu_scale, var_scale );
    target += weibull_lpdf(y_n | beta, alpha); // likelihood
    target += weibull_lpdf(y_p | beta, alpha);
}
generated quantities{
  matrix[2,2] fisherInfo = weib2Par_fisherInfo(y_n, beta, alpha);
  matrix[2,2] varCov= weib2par_varCov(y_n, beta, alpha);
}


