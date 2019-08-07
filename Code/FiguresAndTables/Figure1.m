% File to generate Figure 1
% Consumption equivalent on a rolling window

% In Figure 1 we produce a consumption equivalent (CE) measure of welfare losses (percentage of steady-state consumption) over a rolling window, for the optimal policy (OP- red line).
% We simulate 10000 draws of 2000-periods long series, starting from beliefs corresponding to PLT at time 0, and we calculate the CE welfare loss.
% We take the beliefs in period 1 for each one of the 10000 draws, and from those beliefs we simulate 10000 draws of 2000-periods long series, then calculating the CE welfare loss. We repeat this process for 8000 periods.
% Therefore, the red line at time t is the expected consumption equivalent loss from period-t average beliefs, starting from PLT beliefs at time zero. For comparison, we plot the same CE measure for two Taylor-type rules, that have been proven to drive beliefs respectively to PLT (black line) and IT (purple line) equilibria.
% In order to illustrate the long run  welfare implications of keeping expectations in the PLT and IT equilibriums respectively, for the PLT rule we set the initial beliefs at PLT, and the IT policy is simulated starting from IT beliefs.



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
shocktype = 'series';%'irf'; %

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
intdraws = 10000;
draws = 1;
periods_simulations =8000; 
% length of the Montecarlo welfare experiment
periods_simulations_mc =2000; 


%% SOLVE THE MODEL
disp('Solving the model...') 
[parpolicy, fspace, Grid, max_test] = main_solver(bench_p,projection_parameters, gaintype);



%% MONTECARLO EXERCISE
%% Starting from PLT beliefs

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
sim_parameters_IT = mcparameters(intdraws,draws,periods_simulations, seedzero,...
    gap_zero,0, 0, gamma_zero, shutdownlearning);


% create a Montecarlo simulation for each model
[~,~ , ~,b_pi_simu_MMS,b_gap_simu_MMS,~,~,~,~,~,~,~,~,~,~,~,~] = ...
    simul_fast(parpolicy, fspace, bench_p, model, ...
    gaintype, shocktype, sim_parameters);
[~,~ , ~,b_pi_simu_EH_PLT,b_gap_simu_EH_PLT,~,~,~,~,~,~,~,~,~,~,~,~] = ...
    simul_fast(parpolicy, fspace, bench_p, 'EHCOMM', ...
    gaintype, shocktype, sim_parameters);
[~,~ , ~,b_pi_simu_EH_IT,b_gap_simu_EH_IT,~,~,~,~,~,~,~,~,~,~,~,~] = ...
    simul_fast(parpolicy, fspace, bench_p, 'EHDISCR', ...
    gaintype, shocktype, sim_parameters_IT);


% initialize consumption equivalent series
mcce_MMS_from_PLT = zeros(1, periods_simulations);
mcce_EHCOMM_from_PLT = zeros(1, periods_simulations);
mcce_EHDISCR_from_IT = zeros(1, periods_simulations);

for i = 1: periods_simulations


    % collect Montecarlo parameters into a structure
    sim_parameters_MMS = mcparameters(intdraws,draws,periods_simulations_mc, seedzero,...
        gap_zero,b_gap_simu_MMS(:,i), b_pi_simu_MMS(:,i), gamma_zero, shutdownlearning);

    sim_parameters_EH_PLT = mcparameters(intdraws,draws,periods_simulations_mc, seedzero,...
        gap_zero,b_gap_simu_EH_PLT(:,i), b_pi_simu_EH_PLT(:,i), gamma_zero, shutdownlearning);

    sim_parameters_EH_IT = mcparameters(intdraws,draws,periods_simulations_mc, seedzero,...
        gap_zero,b_gap_simu_EH_IT(:,i), b_pi_simu_EH_IT(:,i), gamma_zero, shutdownlearning);


    disp(['PERIOD ' num2str(i)]);
    % calculate consumption equivalent from simulation
    % OPTIMAL POLICY
    disp('Simulating optimal policy from PLT beliefs...')
    mcce_MMS_from_PLT(:,i) = mc_cons_equiv(parpolicy, fspace, bench_p, model, ...
            gaintype, shocktype, sim_parameters_MMS);
    % EVANS-HONKAPOHJA COMMITMENT
    disp('Simulating EH PLT-inducing policy from PLT beliefs...')
    mcce_EHCOMM_from_PLT(:,i) = mc_cons_equiv(parpolicy, fspace, bench_p, 'EHCOMM', ...
            gaintype, shocktype, sim_parameters_EH_PLT);
    % EVANS-HONKAPOHJA DISCRETION       
    disp('Simulating EH IT-inducing policy from IT beliefs...') 
    mcce_EHDISCR_from_IT(:,i) = mc_cons_equiv(parpolicy, fspace, bench_p, 'EHDISCR', ...
            gaintype, shocktype, sim_parameters_EH_IT);

end
    
%% Figure 1
cd('../Results');

figure(1);
plot(1:length(mcce_MMS_from_PLT),mcce_MMS_from_PLT,'r',...
    1:length(mcce_EHCOMM_from_PLT), mcce_EHCOMM_from_PLT,'k',...
    1:length(mcce_EHDISCR_from_IT), mcce_EHDISCR_from_IT, 'm',...
    'LineWidth',2);
legend('Optimal','EH PLT', 'EH IT');

print -f1 -r600 -depsc2 figure1;
saveas(1,'figure1','fig');
