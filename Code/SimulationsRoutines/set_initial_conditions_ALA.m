function [b_pi_simu_zero , b_gap_simu_zero, ...
            gamma_t_simu_zero, gap_lag_simu_zero, R_simu_zero] = ...
            set_initial_conditions_ALA(model, b_pi_comm, b_x_comm, ...
                b_pi_zero, b_gap_zero, gamma_t_zero, gap_zero, OLS_flag, R_zero)

%% SET_INITIAL_CONDITIONS_ALA sets the inital conditions for simulations, depending on the model and alternative learning algorithm
% Inputs:
% model: which model we are simulating
% b_pi_comm: RE commitment value for b^{\pi}
% b_x_comm: RE commitment value for b^{x}
% b_pi_zero: initial conditions used for other models (no RE)
% b_gap_zero: initial conditions used for other models (no RE)
% gamma_t_zero: initial conditions for gamma
% gap_zero: initial conditions for output gap 
% OLS_flag: flag for OLS learning ('yes', 'no')
% R_zero: used for OLS learning

%% Set initial conditions
if strcmp(model,'RECOMM')
    b_pi_simu_zero = b_pi_comm;
    b_gap_simu_zero = b_x_comm;
elseif strcmp(model,'REDISCR')
    b_pi_simu_zero = 0;
    b_gap_simu_zero = 0;
else
    b_pi_simu_zero = b_pi_zero;
    b_gap_simu_zero = b_gap_zero;
end
gamma_t_simu_zero = gamma_t_zero;
gap_lag_simu_zero = gap_zero;
if strcmp(OLS_flag, 'yes')
    R_simu_zero = R_zero ;
else
    R_simu_zero = 1;
end

end
