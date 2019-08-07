%% Figure 5a
% Impulse effect on output gap as a function of beliefs, for OP and PLT policies


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
periods_simulations = 1;

%% SOLVE THE MODEL
disp('Solving the model...')
[parpolicy, fspace, Grid, max_test] = main_solver(bench_p,projection_parameters, gaintype);


%% SET INITIAL conditions
numberInitialConditions = 20;
b_pi_zero_vector = linspace(0,bench_p.b_pi_comm,numberInitialConditions);
b_gap_zero_vector = bench_p.b_x_comm.*ones(1,numberInitialConditions); %  linspace(0,bench_p.b_x_comm,numberInitialConditions); % 


%% ALLOCATE MEMORY
gap_MMS = zeros(1,numberInitialConditions);
gap_PLT = zeros(1,numberInitialConditions);
gap_IT = zeros(1,numberInitialConditions);


%% IRFs starting from different beliefs

for i = 1: numberInitialConditions

    % initial conditions
    gap_zero = 0;
    b_gap_zero = b_gap_zero_vector(i);
    b_pi_zero =  b_pi_zero_vector(i);
    % b_gap_zero = -bench_p.b_x_comm + 2*bench_p.b_x_comm.*rand(intdraws,1) ;
    % b_pi_zero =  bench_p.b_pi_min + (bench_p.b_pi_max-bench_p.b_pi_min).*rand(intdraws,1) ;
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

    % create IRF for PLT
    [gap_lag_simu_EH_PLT, pi_simu_EH_PLT , gap_simu_EH_PLT,b_pi_simu_EH_PLT,b_gap_simu_EH_PLT,...
        gamma_t_simu_EH_PLT,lambda1_simu_EH_PLT, welfare_cumul_EH_PLT, welfare_x_cumul_EH_PLT, ...
        welfare_pi_cumul_EH_PLT, welfare_cumul_final_EH_PLT, welfare_x_cumul_final_EH_PLT , ...
        welfare_pi_cumul_final_EH_PLT, costpushshock_EH_PLT, welfare_inst_EH_PLT, welfare_x_inst_EH_PLT, ...
        welfare_pi_inst_EH_PLT ] = simul_fast(parpolicy, fspace, bench_p, 'EHCOMM', ...
        gaintype, shocktype, sim_parameters);


    % create IRF for OPTgamma0
    [gap_lag_simu_MMSgamma0, pi_simu_MMSgamma0 , gap_simu_MMSgamma0,b_pi_simu_MMSgamma0,b_gap_simu_MMSgamma0,...
    gamma_t_simu_MMSgamma0,lambda1_simu_MMSgamma0, welfare_cumul_MMSgamma0, welfare_x_cumul_MMSgamma0, ...
    welfare_pi_cumul_MMSgamma0, welfare_cumul_final_MMSgamma0, welfare_x_cumul_final_MMSgamma0 , ...
    welfare_pi_cumul_final_MMSgamma0, costpushshock_MMSgamma0, welfare_inst_MMSgamma0, welfare_x_inst_MMSgamma0, ...
    welfare_pi_inst_MMSgamma0 ] = simul_fast(parpolicy, fspace, bench_p, 'MMSgamma0', ...
    gaintype, shocktype, sim_parameters);

    gap_MMS(i) = gap_simu_MMS;
    gap_PLT(i) = gap_simu_EH_PLT;
    gap_MMSgamma0(i) = gap_simu_MMSgamma0;

end

%% Figure 5
%% Panel (a), optimal policy and PLT
figure(30);
plot(b_pi_zero_vector,gap_MMS,'b-',...
    b_pi_zero_vector,gap_PLT,'r--',...
    b_pi_zero_vector,gap_MMSgamma0,'g-*',...
    'LineWidth',2);

    
xlabel('b_0^\pi', 'FontSize',16,...
       'FontWeight','bold'); ylabel('x_1', 'FontSize',16,...
       'FontWeight','bold');
set(gca,'FontSize',20);
legend('OP','PLT', 'OP w/ \gamma = 0');


cd('../Results')
print -f30 -r600 -depsc2 figure5a;
saveas(30,'figure5a','fig');

