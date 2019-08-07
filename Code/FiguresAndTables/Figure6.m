% File to generate Figures 6a, 6b and 6c
% IRFs to one standard deviation shock

%% HOUSEKEEPING
clear all; close all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters for the simulation algorithm %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Gain type: decreasing or constant gain
gaintype = 'const'; %'decr'; %

% In some simulations we shut down learning (setting this parameter to 1)
shutdownlearning = 0; % 1; %

% This parameter selects if we are calculating a IRF or one (or many)
% series for montecarlo analysis
shocktype = 'irf'; %'series';%
% This is the model we simulate: optimal policy ('MMS'), or Evans and
% Honkapoihja either with commitment ('EHCOMM') or with discretion
% ('EHDISCR'), rational expectations either with commitment ('RECOMM') or
% with discretion ('REDISCR')
model = 'MMS'; % 'RECOMM'; % 'EHCOMM'; %

% benchmark parameter for the model
bench_p = benchmarkparameters(gaintype);

% projection parameters: can use 'cheb', 'lin' or 'spli';
% if using splines, the default spline order is cubic;
% if not using splines, or for default spline order, use []
% as second input
projection_parameters = projectionparameters('cheb', []);

%% Montecarlo parameters
% draws and length of the simulation
intdraws = 1;
draws = 1;
periods_simulations = 8;

%% SOLVE THE MODEL
disp('Solving the model...')
[parpolicy, fspace, Grid, max_test] = main_solver(bench_p,projection_parameters, gaintype);


%% IRFs starting from PLT beliefs

% initial conditions
gap_zero = 0;
b_gap_zero = bench_p.b_x_comm;
b_pi_zero =  bench_p.b_pi_comm;
gamma_zero = 1;
% seed for the random number generator
seedzero = 1;

% collect Montecarlo parameters into a structure
sim_parameters = mcparameters(intdraws,draws,periods_simulations, seedzero,...
    gap_zero,b_gap_zero, b_pi_zero, gamma_zero, shutdownlearning);


% create IRF for optimal policy
[gap_lag_simu_MMS, pi_simu_MMS , gap_simu_MMS,b_pi_simu_MMS,b_gap_simu_MMS,...
    gamma_t_simu_MMS,lambda1_simu_MMS, welfare_cumul_MMS, welfare_x_cumul_MMS, ...
    welfare_pi_cumul_MMS, welfare_cumul_final_MMS, welfare_x_cumul_final_MMS , ...
    welfare_pi_cumul_final_MMS, costpushshock_MMS, welfare_inst_MMS, welfare_x_inst_MMS, ...
    welfare_pi_inst_MMS ] = simul_fast(parpolicy, fspace, bench_p, model, ...
    gaintype, shocktype, sim_parameters);
price_MMS = pricefrominflation(pi_simu_MMS,periods_simulations); % this calculate the price series from inflation

% create IRF for EH IT
[gap_lag_simu_EH_IT, pi_simu_EH_IT , gap_simu_EH_IT,b_pi_simu_EH_IT,b_gap_simu_EH_IT,...
    gamma_t_simu_EH_IT,lambda1_simu_EH_IT, welfare_cumul_EH_IT, welfare_x_cumul_EH_IT, ...
    welfare_pi_cumul_EH_IT, welfare_cumul_final_EH_IT, welfare_x_cumul_final_EH_IT , ...
    welfare_pi_cumul_final_EH_IT, costpushshock_EH_IT, welfare_inst_EH_IT, welfare_x_inst_EH_IT, ...
    welfare_pi_inst_EH_IT ] = simul_fast(parpolicy, fspace, bench_p, 'EHDISCR', ...
    gaintype, shocktype, sim_parameters);
price_EH_IT = pricefrominflation(pi_simu_EH_IT,periods_simulations); % this calculate the price series from inflation

% create IRF for EH PLT
[gap_lag_simu_EH_PLT, pi_simu_EH_PLT , gap_simu_EH_PLT,b_pi_simu_EH_PLT,b_gap_simu_EH_PLT,...
    gamma_t_simu_EH_PLT,lambda1_simu_EH_PLT, welfare_cumul_EH_PLT, welfare_x_cumul_EH_PLT, ...
    welfare_pi_cumul_EH_PLT, welfare_cumul_final_EH_PLT, welfare_x_cumul_final_EH_PLT , ...
    welfare_pi_cumul_final_EH_PLT, costpushshock_EH_PLT, welfare_inst_EH_PLT, welfare_x_inst_EH_PLT, ...
    welfare_pi_inst_EH_PLT ] = simul_fast(parpolicy, fspace, bench_p, 'EHCOMM', ...
    gaintype, shocktype, sim_parameters);
price_EH_PLT = pricefrominflation(pi_simu_EH_PLT,periods_simulations); % this calculate the price series from inflation


% create IRF for EH PLT
[gap_lag_simu_MMS_gamma0, pi_simu_MMS_gamma0 , gap_simu_MMS_gamma0,b_pi_simu_MMS_gamma0,b_gap_simu_MMS_gamma0,...
    gamma_t_simu_MMS_gamma0,lambda1_simu_MMS_gamma0, welfare_cumul_MMS_gamma0, welfare_x_cumul_MMS_gamma0, ...
    welfare_pi_cumul_MMS_gamma0, welfare_cumul_final_MMS_gamma0, welfare_x_cumul_final_MMS_gamma0 , ...
    welfare_pi_cumul_final_MMS_gamma0, costpushshock_MMS_gamma0, welfare_inst_MMS_gamma0, welfare_x_inst_MMS_gamma0, ...
    welfare_pi_inst_MMS_gamma0 ] = simul_fast(parpolicy, fspace, bench_p, 'MMSgamma0', ...
    gaintype, shocktype, sim_parameters);
price_MMS_gamma0 = pricefrominflation(pi_simu_MMS_gamma0,periods_simulations); % this calculate the price series from inflation

%% Figure 6
%% Panel (a)
figure(2);
plot(0:length(pi_simu_MMS),[0 pi_simu_MMS],'b',...
    0:length(pi_simu_MMS_gamma0),[0 pi_simu_MMS_gamma0],'g-*',...
    0:length(pi_simu_EH_PLT),[0 pi_simu_EH_PLT],'r--');
legend('Optimal','Optimal \gamma = 0','EH PLT');
axis([0 length(b_pi_simu_EH_PLT) -.2 .2])

%% Panel (b)
figure(3);
plot(0:length(gap_simu_MMS),[0 gap_simu_MMS],'b',...
    0:length(gap_simu_MMS_gamma0),[0 gap_simu_MMS_gamma0],'g-*',...
    0:length(gap_simu_EH_PLT),[0 gap_simu_EH_PLT],'r--');
legend('Optimal','Optimal \gamma = 0','EH PLT');
axis([0 length(b_pi_simu_EH_PLT) -.2 .2])

%% Panel (c)
figure(4);
plot(0:length(price_MMS),[0 price_MMS],'b',...
    0:length(price_MMS_gamma0),[0 price_MMS_gamma0],'g-*',...
    0:length(price_EH_PLT),[0 price_EH_PLT],'r--');
legend('Optimal','Optimal \gamma = 0','EH PLT');
axis([0 length(b_pi_simu_EH_PLT) -.2 .2])

%% Save graphs
cd('../Results')
print -f2 -r600 -depsc2 figure6a;
saveas(2,'figure6a','fig');

print -f3 -r600 -depsc2 figure6b;
saveas(3,'figure6b','fig');

print -f4 -r600 -depsc2 figure6c;
saveas(4,'figure6c','fig');