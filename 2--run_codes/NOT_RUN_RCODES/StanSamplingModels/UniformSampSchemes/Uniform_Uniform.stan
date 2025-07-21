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
  // int  p; 
  vector[n] y_n;
  real mu_shape;
  real var_shape;
  real mu_scale;
  real var_scale;
 }
 transformed data{
    // real lbd_shp = mu_shape - 0.5* sqrt(12*var_shape);
    // real ubd_shp = mu_shape + 0.5* sqrt(12*var_shape);
    // real lbd_scl = mu_scale - 0.5* sqrt(12*var_scale);
    // real ubd_scl = mu_scale + 0.5* sqrt(12*var_scale);
}
 parameters{
   real <lower=0> beta;
   real <lower=0>  alpha;
   // real <lower=0> y_pred[p];
 }
model{
  beta ~ uniform(0, 1);
  alpha ~ uniform(0, 1);
  // target += uniform_lpdf(beta | 0, 1); // priors
  // target += uniform_lpdf(alpha | 0, 1);
  target += weibull_lpdf(y_n | beta, alpha); // likelihood
  // target += weibull_lpdf(y_pred | beta, alpha) ;
}
generated quantities{
  matrix[2,2] fisherInfo = weib2Par_fisherInfo(y_n, beta, alpha);
  matrix[2,2] varCov= weib2par_varCov(y_n, beta, alpha);
}
