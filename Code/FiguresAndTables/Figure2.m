%% Assess the impact of different initial beliefs on the dynamics of the learning coefficients


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
scaler_b_pi = [1;.25;0]; 
scaler_b_gap =[1;.25;0]; 
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
    % calculate consumption equivalent from simulation
    % OPTIMAL POLICY
    disp('Simulating optimal policy...')
    [bcoeff_MMS(:,:,i), ~, ~, ~, ~ ] = mc_allocations(parpolicy, fspace, bench_p, model, ...
        gaintype, shocktype, sim_parameters);
    
end

%% Produce graphs
% first, the lines!
colorsused = {'r', 'b', 'k'};
linesdrawn = {'-', '--', '.'};

%% Figure 2a
figure(1);
for i = 1: length(initial_beliefs_b_gap)
    textleg(i)={['Beliefs: b^{\pi} = ' num2str(initial_beliefs_b_pi(i))  ' , b^x = ' num2str(initial_beliefs_b_gap(i))]}; 
    plot(0:periods_simulations,bcoeff_MMS(1,:,i), [colorsused{i}, linesdrawn{i}], 'LineWidth', 2); 
    hold on;

end
title('$b_t^{\pi}$ dynamics for different initial beliefs', 'interpreter', 'latex');
hold off;

legend(textleg);


%% Figure 2b
figure(2);
for i = 1: length(initial_beliefs_b_gap)
    textleg(i)={['Beliefs b^x: b^{\pi} = ' num2str(initial_beliefs_b_pi(i))  ' , b^x = ' num2str(initial_beliefs_b_gap(i))]}; 
    plot(0:periods_simulations, bcoeff_MMS(2,:,i), [colorsused{i}, linesdrawn{i}], 'LineWidth', 2);
    hold on;

end
title('$b_t^{x}$ dynamics for different initial beliefs', 'interpreter', 'latex');
hold off;

legend(textleg);


%% Save results
cd('../Results')
print -f1 -r600 -depsc2 figure2a;
saveas(1,'figure2a','fig');

print -f2 -r600 -depsc2 figure2b;
saveas(2,'figure2b','fig');
