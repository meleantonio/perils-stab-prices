function [welfare_cumul_MMS, welfare_x_cumul, ...
            welfare_pi_cumul, welfare_inst, ...
            welfare_x_inst, welfare_pi_inst, gap_lag_simu, ...
            pi_simu, gap_simu, b_pi_simu, b_gap_simu, gamma_t_simu, lambda1_simu] = ...
            allocate_memory_for_series(intdraws,periods_simulations)

% ALLOCATE_MEMORY_FOR_SERIES allocates memory for the series to be generated

welfare_cumul_MMS = zeros(intdraws, periods_simulations);
welfare_x_cumul = zeros(intdraws, periods_simulations);
welfare_pi_cumul=zeros(intdraws, periods_simulations);
welfare_inst = zeros(intdraws, periods_simulations);
welfare_x_inst  = zeros(intdraws, periods_simulations);
welfare_pi_inst =zeros(intdraws, periods_simulations);
gap_lag_simu = zeros(intdraws, periods_simulations+1);
pi_simu =zeros(intdraws, periods_simulations);
gap_simu = zeros(intdraws, periods_simulations);
b_pi_simu = zeros(intdraws, periods_simulations+1);
b_gap_simu = zeros(intdraws, periods_simulations+1);
gamma_t_simu = zeros(intdraws, periods_simulations+1);
lambda1_simu =zeros(intdraws, periods_simulations);
  

end