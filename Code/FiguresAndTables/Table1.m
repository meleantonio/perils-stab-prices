%% Creating table 1: welfare losses

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
periods_simulations =2000; 


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


% calculate consumption equivalent from simulation
% OPTIMAL POLICY
disp('Simulating optimal policy from PLT beliefs...')
mcce_MMS_from_PLT = mc_cons_equiv(parpolicy, fspace, bench_p, model, ...
        gaintype, shocktype, sim_parameters);
% EVANS-HONKAPOHJA COMMITMENT
disp('Simulating EH PLT-inducing policy from PLT beliefs...')
mcce_EHCOMM_from_PLT = mc_cons_equiv(parpolicy, fspace, bench_p, 'EHCOMM', ...
        gaintype, shocktype, sim_parameters);

% RATIOS:
ratio_PLT_OPT_from_PLT = mcce_EHCOMM_from_PLT/mcce_MMS_from_PLT;


%% Starting from IT beliefs

% initial conditions
gap_zero = 0;
b_gap_zero = 0;
b_pi_zero =  0;
gamma_zero = 1;
% seed for the random number generator
seedzero = 1;

% collect Montecarlo parameters into a structure
sim_parameters = mcparameters(intdraws,draws,periods_simulations, seedzero,...
gap_zero,b_gap_zero, b_pi_zero, gamma_zero, shutdownlearning);


% calculate consumption equivalent from simulation
% OPTIMAL POLICY
disp('Simulating optimal policy from IT beliefs...')
mcce_MMS_from_IT = mc_cons_equiv(parpolicy, fspace, bench_p, model, ...
        gaintype, shocktype, sim_parameters);
% EVANS-HONKAPOHJA COMMITMENT
disp('Simulating EH PLT-inducing policy from IT beliefs...')
mcce_EHCOMM_from_IT = mc_cons_equiv(parpolicy, fspace, bench_p, 'EHCOMM', ...
        gaintype, shocktype, sim_parameters);

% RATIOS:
ratio_PLT_OPT_from_IT = mcce_EHCOMM_from_IT/mcce_MMS_from_IT;


%% EXPORT TO LATEX

cd('../Results');
disp('Creating LaTeX file table1.tex...') 

% name rows and columns
rowLabels = { 'PLT'; 'IT'};
columnLabels = { 'OPT', 'PLT', 'ratioPLTtoOPT'};

% collect data in the right format
OPT = [mcce_MMS_from_PLT ; mcce_MMS_from_IT ];
PLT = [mcce_EHCOMM_from_PLT; mcce_EHCOMM_from_IT];
ratio = [ratio_PLT_OPT_from_PLT ; ratio_PLT_OPT_from_IT]; 

% create table
table1 = table(OPT, PLT, ratio, 'RowNames', rowLabels, 'VariableNames', columnLabels );

% save it in Results as a tex file (needs some formatting)
table2latex(table1, 'table1.tex')

disp('Done!') 

