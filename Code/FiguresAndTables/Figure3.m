%% Squared forecast errors, inflation and output gap


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
intdraws = 10000;
draws = 1;
periods_simulations =5000;


%% SOLVE THE MODEL
disp('Solving the model...')
[parpolicy, fspace, Grid, max_test] = main_solver(bench_p,projection_parameters, gaintype);

%% Initial beliefs' vector
scaler_b_pi = [1; 0]; 
scaler_b_gap =[1; 0]; 
initial_beliefs_b_gap = scaler_b_gap.*bench_p.b_x_comm;  
initial_beliefs_b_pi  = scaler_b_pi.*bench_p.b_pi_comm; 


%% LOOP MONTECARLO EXERCISE
for i = 1: length(initial_beliefs_b_gap)   
    
    % initial conditions
    gap_zero = 0;
    b_gap_zero = initial_beliefs_b_gap(i);
    b_pi_zero =  initial_beliefs_b_pi(i);
    gamma_zero = 1;
    % seed for the random number generator
    seedzero = 13;
    
    % collect Montecarlo parameters into a structure
    sim_parameters = mcparameters(intdraws,draws,periods_simulations, seedzero,...
        gap_zero,b_gap_zero, b_pi_zero, gamma_zero, shutdownlearning);
    disp(' ');
    disp(sprintf('Initial beliefs : b_{x,0} = %g, b_{pi,0} = %g ', b_gap_zero, b_pi_zero ));
    % OPTIMAL POLICY
    disp('Simulating optimal policy...')
    [bcoeff_MMS(:,:,i), prediction_errors_bcoeff(:,:,i), ...
    allocations(:,:,i), forecast_errors_allocations(:,:,i) forecast_errors_allocations_squared(:,:,i) ] = mc_allocations(parpolicy, fspace, bench_p, model, ...
        gaintype, shocktype, sim_parameters);


    % EVANS-HONKAPOHJA PLT
    [bcoeff_PLT(:,:,i), prediction_errors_bcoeff_PLT(:,:,i), ...
    allocations_PLT(:,:,i), forecast_errors_allocations_PLT(:,:,i) forecast_errors_allocations_squared_PLT(:,:,i) ] = mc_allocations(parpolicy, fspace, bench_p, 'EHCOMM', ...
        gaintype, shocktype, sim_parameters);
    
end


%% Figure 3a
figure(1);
plot(1:periods_simulations-1, forecast_errors_allocations_squared(1,:,1) , 'r', ...
1:periods_simulations-1, forecast_errors_allocations_squared(1,:, end), 'b', ...
1:periods_simulations-1, forecast_errors_allocations_squared_PLT(1,:, 1), 'g', 'LineWidth', 2);
legend('OPT from PLT beliefs', 'OPT from IT beliefs', 'PLT from PLT beliefs')
title('Squared Forecast Error for $\pi_t$', 'interpreter', 'latex');


%% Figure 3b
figure(2);
plot(1:periods_simulations-1, forecast_errors_allocations_squared(2,:,1) , 'r', ...
1:periods_simulations-1, forecast_errors_allocations_squared(2,:, end), 'b', ...
1:periods_simulations-1, forecast_errors_allocations_squared_PLT(2,:, 1), 'g', 'LineWidth',2);
legend('OPT from PLT beliefs', 'OPT from IT beliefs', 'PLT from PLT beliefs')
title('Squared Forecast Error for $x_t$', 'interpreter', 'latex');

rwerq
%% Save results
cd('../Results')
print -f1 -r600 -depsc2 figure3a;
saveas(1,'figure3a','fig');

print -f2 -r600 -depsc2 figure3b;
saveas(2,'figure3b','fig');
