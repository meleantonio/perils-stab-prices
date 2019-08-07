%% Forecast errors


%% HOUSEKEEPING
clear all; close all; clc;

%% Parameters for the simulation algorithm 

%% Gain type: decreasing or constant gain
gaintype = 'const'; % 'decr'; %

%% In some simulations we shut down learning (setting this parameter to 1)
shutdownlearning = 0; % 1; %

%% This parameter selects if we are calculating a IRF or one (or many) series for montecarlo analysis
shocktype = 'series';%'irf'; %

%% This is the model we simulate: optimal policy ('MMS'), or Evans and Honkapoihja either with commitment ('EHCOMM') or with discretion ('EHDISCR'), rational expectations either with commitment ('RECOMM') or with discretion ('REDISCR')
model = 'MMS'; % 'RECOMM'; % 'EHCOMM'; %

%% benchmark parameter for the model
bench_p = benchmarkparameters(gaintype);

%% projection parameters: can use 'cheb', 'lin' or 'spli'; if using splines, the default spline order is cubic; if not using splines, or for default spline order, use [] as second input
projection_parameters = projectionparameters('cheb', []);

%% Montecarlo parameters: draws and length of the simulation
intdraws = 1;
draws = 1;
periods_simulations =100;


%% SOLVE THE MODEL
disp('Solving the model...')
[parpolicy, fspace, Grid, max_test] = main_solver(bench_p,projection_parameters, gaintype);

 
% initial conditions
gap_zero = 0;
b_gap_zero = bench_p.b_x_comm;
b_pi_zero =  bench_p.b_pi_comm; 
gamma_zero = 1;
% seed for the random number generator
seedzero = 13;

% collect Montecarlo parameters into a structure
sim_parameters = mcparameters(intdraws,draws,periods_simulations, seedzero,...
    gap_zero,b_gap_zero, b_pi_zero, gamma_zero, shutdownlearning);
disp(' ');
% OPTIMAL POLICY
disp('Simulating optimal policy...')
[bcoeff(:,:), prediction_errors_bcoeff(:,:), ...
allocations(:,:), forecast_errors_allocations(:,:) forecast_errors_allocations_squared(:,:) ] = mc_allocations(parpolicy, fspace, bench_p, model, ...
    gaintype, shocktype, sim_parameters);


%% Figure 4
figure(1);
plot(1:periods_simulations-1, forecast_errors_allocations(1,:) , 'r--', ...
1:periods_simulations-1, forecast_errors_allocations(2, :), 'b-', 'LineWidth', 2);
legend('FE \pi_t', 'FE x_t')
title('Forecast Errors for $x_t$ and $\pi_t$ from PLT beliefs', 'interpreter', 'latex');

%% Save results
cd('../Results')
print -f1 -r600 -depsc2 figure4;
saveas(1,'figure4','fig');
