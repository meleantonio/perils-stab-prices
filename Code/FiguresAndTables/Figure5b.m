%% Figure 5b
% Impact of a one-standard deviation shock on output gap, from different initial conditions, and different gammas



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

% gamma vector
gamma_vector =  [.01 .02 0.05 .1 .2 ];

%% SET INITIAL conditions
numberInitialConditions = 10;
b_pi_zero_vector = linspace(0,bench_p.b_pi_comm,numberInitialConditions);
b_gap_zero_vector = bench_p.b_x_comm.*ones(1,numberInitialConditions); %  linspace(0,bench_p.b_x_comm,numberInitialConditions); %

%% ALLOCATE MEMORY
gap_MMS = zeros(numberInitialConditions,length(gamma_vector));
gap_PLT = zeros(numberInitialConditions,length(gamma_vector));
gap_IT =  zeros(numberInitialConditions,length(gamma_vector));

%% LOOP OVER GAMMAS
for l = 1:length(gamma_vector)
    
    bench_p.gam = gamma_vector(l);
    
    %% SOLVE THE MODEL
    disp('Solving the model...')
    [parpolicy, fspace, Grid, max_test] = main_solver(bench_p,projection_parameters, gaintype);
    
    
    
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
        [gap_lag_simu_MMS,...
            pi_simu_MMS  , ...
            gap_simu_MMS ,...
            b_pi_simu_MMS,...
            b_gap_simu_MMS,...
            gamma_t_simu_MMS,...
            lambda1_simu_MMS, ...
            welfare_cumul_MMS, ...
            welfare_x_cumul_MMS, ...
            welfare_pi_cumul_MMS, ...
            welfare_cumul_final_MMS, ...
            welfare_x_cumul_final_MMS , ...
            welfare_pi_cumul_final_MMS, ...
            costpushshock_MMS, ...
            welfare_inst_MMS, ...
            welfare_x_inst_MMS, ...
            welfare_pi_inst_MMS] = simul_fast(parpolicy, fspace, bench_p, model, ...
            gaintype, shocktype, sim_parameters);
        
        % create IRF for PLT
        [gap_lag_simu_EH_PLT, pi_simu_EH_PLT , gap_simu_EH_PLT,b_pi_simu_EH_PLT,b_gap_simu_EH_PLT,...
            gamma_t_simu_EH_PLT,lambda1_simu_EH_PLT, welfare_cumul_EH_PLT, welfare_x_cumul_EH_PLT, ...
            welfare_pi_cumul_EH_PLT, welfare_cumul_final_EH_PLT, welfare_x_cumul_final_EH_PLT , ...
            welfare_pi_cumul_final_EH_PLT, costpushshock_EH_PLT, welfare_inst_EH_PLT, welfare_x_inst_EH_PLT, ...
            welfare_pi_inst_EH_PLT ] = simul_fast(parpolicy, fspace, bench_p, 'EHCOMM', ...
            gaintype, shocktype, sim_parameters);
        
        gap_MMS(i,l) = gap_simu_MMS;
        % gap_IT(i) = gap_simu_EH_IT;
        gap_PLT(i,l) = gap_simu_EH_PLT;
        
    end
end

%% Figure 5
%% Panel (b), gammas and initial conditions
figure(41);

symbols = {'bo-','r*-', 'k^-', 'gd-','ms-'};

for gamma_i = 1:length(gamma_vector)
    plot(b_pi_zero_vector,gap_MMS(:,gamma_i), symbols{gamma_i},'LineWidth',1)
    legendInfo{gamma_i}=['\gamma = ' num2str(gamma_vector(gamma_i))];
    hold on;
end
hold off;
legend(legendInfo);

xlabel('b_0^\pi', 'FontSize',16,...
    'FontWeight','bold'); ylabel('x_1', 'FontSize',16,...
    'FontWeight','bold');
set(gca,'FontSize',20);

cd('../Results')

print -f41 -r600 -depsc2 figure5b;
saveas(41,'figure5b','fig');

