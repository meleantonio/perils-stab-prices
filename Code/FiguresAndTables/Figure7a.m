% File to generate Figure 7a
% Consumption equivalent on a rolling window (as in Figure 1)
% Sensitivity analysis: we vary the learning parameter gam to check for consistency

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
periods_simulations =4000;
% length of the Montecarlo welfare experiment
periods_simulations_mc =1000;



% Set the vector containing the gammas we want to test
gamma_list = [.01; .08];

%% LOOP OVER GAMMAS
for index_gamma = 1:length(gamma_list)
    
    bench_p.gam = gamma_list(index_gamma);
    
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

    % create a Montecarlo simulation for each model
    [~,~ , ~,b_pi_simu_MMS,b_gap_simu_MMS,~,~,~,~,~,~,~,~,~,~,~,~] = ...
        simul_fast(parpolicy, fspace, bench_p, model, ...
        gaintype, shocktype, sim_parameters);
    [~,~ , ~,b_pi_simu_EH_PLT,b_gap_simu_EH_PLT,~,~,~,~,~,~,~,~,~,~,~,~] = ...
        simul_fast(parpolicy, fspace, bench_p, 'EHCOMM', ...
        gaintype, shocktype, sim_parameters);
    
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
        
    end
    
    %% Figure 4a: build it recursively
    if index_gamma == 1
        figure(1);
    end
    
    plot(1:length(mcce_MMS_from_PLT),mcce_MMS_from_PLT,...
        'DisplayName',['Opt, \gamma = ' num2str(gamma_list(index_gamma))], 'LineWidth',2);
    hold on;

    plot(1:length(mcce_EHCOMM_from_PLT), mcce_EHCOMM_from_PLT,...
        'DisplayName',['EH PLT, \gamma = ' num2str(gamma_list(index_gamma))], 'LineWidth',2);
    axis([0 periods_simulations 0 0.001]);   
    hold on;
    
end
hold off;
legend('show', 'Location','southeast')

cd('../Results')

print -f1 -r600 -depsc2 figure7a;
saveas(1,'figure7a','fig');
