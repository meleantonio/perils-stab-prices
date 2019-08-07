function [states, gamma_t_simu] = states_in_period_i_ALA(gap_lag_simu,b_pi_simu,gamma_t_simu,gam,costpushshock,gaintype,shutdownlearning, ~)

% STATES_IN_PERIOD_I_ALA collects the states and the gain parameter for each period, depending on the gaintype (alternative learning algorithms). [It is in fact identical to states_in_period_i, so in a next version can be dropped]

        % set state variables for the different cases
        if strcmp(gaintype,'decr')
            states = [gap_lag_simu ...
                b_pi_simu  gamma_t_simu  costpushshock ];
        elseif strcmp(gaintype,'const')
            if shutdownlearning ==1
                gamma_t_simu  = 0;
            else
                gamma_t_simu  = gam;
            end
            states = [gap_lag_simu  ...
                b_pi_simu   costpushshock ];
        end
        
end
