function [bcoeff, prediction_error, allocations, forecast_errors, forecast_errors_squared] = mc_allocations(parpolicy, fspace, bench_p, model, ...
    gaintype, shocktype, sim_parameters)

% MC_ALLOCATIONS simulates the model with a Montecarlo and returns the average series over time for the b's (learning coefficients), the allocations x_t and pi_t, and the prediction errors for them. The outputs are the series for the average value of the b's, the allocations (pi and output gap), and the respective forecast errors. 

% Inputs:

% parpolicy: coefficients for the basis functions
% fspace: functional space for the basis function
% bench_p: parameters for the model
% model: the model we are simulating
% gaintype: either 'decr' or 'const'
% shocktype: if series or irf
% sim_parameters: parameters for the simulations

[gap_lag_mc, pi_mc, gap_mc, b_pi_mc, b_gap_mc, ~ ,~ , ~ , ~ ,  ~ , ~ , ~  , ...
    ~ , ~ , ~ , ~ ,  ~  ] = simul_fast(parpolicy, fspace, bench_p, model, ...
        gaintype, shocktype, sim_parameters);

prediction_error_b_pi = b_pi_mc(:,2:end) - b_pi_mc(:,1:end-1); 
prediction_error_b_gap = b_gap_mc(:,2:end) - b_gap_mc(:,1:end-1); 

forecast_errors_gap = gap_mc(:,2:end) - b_gap_mc(:,1:end-2).*gap_lag_mc(:,1:end-2); 
forecast_errors_pi = pi_mc(:,2:end) - b_pi_mc(:,1:end-2).*gap_lag_mc(:,1:end-2); 

forecast_errors_gap_squared = forecast_errors_gap.^2;
forecast_errors_pi_squared = forecast_errors_pi.^2;


% average out
b_pi_average = mean(b_pi_mc, 1);
b_gap_average = mean(b_gap_mc, 1);

prediction_error_b_pi_average = mean(prediction_error_b_pi,1);
prediction_error_b_gap_average = mean(prediction_error_b_gap,1);

gap_mc_average = mean(gap_mc, 1);
pi_mc_average = mean(pi_mc, 1);

forecast_errors_gap_average = mean(forecast_errors_gap, 1) ; 
forecast_errors_pi_average = mean(forecast_errors_pi, 1); 

forecast_errors_gap_squared_average = mean(forecast_errors_gap_squared, 1) ; 
forecast_errors_pi_squared_average = mean(forecast_errors_pi_squared, 1); 


% collect series in individual outputs
bcoeff = [b_pi_average; b_gap_average];
prediction_error = [prediction_error_b_pi_average; prediction_error_b_gap_average]; 
allocations = [pi_mc_average; gap_mc_average ];
forecast_errors = [forecast_errors_pi_average; forecast_errors_gap_average ];
forecast_errors_squared = [forecast_errors_pi_squared_average; forecast_errors_gap_squared_average ];

end