function parameters = variance_by_a_factor(parameters, variance_factor)
% VARIANCE_BY_A_FACTOR multiplies the standard deviation (and all related parameters) by a factor. Used to check robustness and for the ZLB experiments

parameters.sigma = parameters.sigma*variance_factor;
parameters.sigeps = parameters.sigeps*variance_factor;
parameters.thmin = parameters.thmin*variance_factor;
parameters.thmax = parameters.thmax*variance_factor;
[parameters.QuadrPoints,parameters.QuadrWeights] = qnwnorm(parameters.nQuadr,0,parameters.sigma^2);

end