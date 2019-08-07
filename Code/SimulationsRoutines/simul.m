function [gap_lag_simu, pi_simu , gap_simu,b_pi_simu,b_gap_simu,...
    gamma_t_simu,lambda1_simu, welfare_cumul_MMS, welfare_x_cumul, ...
    welfare_pi_cumul, welfare_cumul_MMS_final, welfare_x_cumul_final , ...
    welfare_pi_cumul_final, costpushshock, welfare_inst, welfare_x_inst, ...
    welfare_pi_inst ] = simul(parpolicy, fspace, bench_p, model, ...
        gaintype, shocktype, sim_parameters)

% SIMUL simulates series or IRFs for the MMS paper.
%
% INPUTS:
%     parpolicy: vector of coefficients for the interpolation
%     fspace: functional space created with Miranda Fackler Compecon Toolbox
%     bench_p: benchmark parameters
%     model: this is the model we are running, can be 'MMS', 'EHCOMM', 'EHDISCR', 'RECOMM',
%     'REDISCR'
%     gaintype: can be set to 'decr' for decreasing gain, or 'const'
%         for constant gain
%     shocktype: can be set to 'irf' for impulse response function or
%         'series' for general montecarlo simulation
%     sim_parameters: parameters for the simulation
%
% OUTPUT: SERIES FOR
%
%     gap_lag_simu: output gap lagged
%     pi_simu: inflation
%     gap_simu: output gap
%     b_pi_simu: learning coefficient for inflation
%     b_gap_simu: learning coefficient for output gap
%     gamma_t_simu: learning parameter series (this willbe constant if
%         constant gain is chosen
%     lambda1_simu: Lagrange multiplier on the law of motion for learning
%         coefficient on output gap
%     welfare_cumul_MMS: series of the welfare
%     welfare_x_cumul: series of the welfare from output gap
%     welfare_pi_cumul: series of the welfare from inflation
%     welfare_cumul_MMS_final: average series of the welfare
%     welfare_x_cumul_final: average series of the welfare from output gap
%     welfare_pi_cumul_final: average series of the welfare from inflation
%     costpushshock: the cost push shock
%     welfare_inst: instantanous undiscounted loss
%     welfare_x_inst : instantanous undiscounted loss from output gap
%     welfare_pi_inst : instantanous undiscounted loss form inflation

%% PARAMETERS
% give individual names to benchmark parameters
[alpha, betta, ~, ~, gam, ~, ~, ~, ...
    ~, ~, ~, ~, b_x_comm, ...
    b_pi_comm, ~, ~, ~, ...
    ~, ~, ~, ...
    ~, ~, ~, ~, ...
    ~, ~, ~, ~, ...
    sigeps, ~, ~, ~, ~, imp] = translate_parameters(bench_p);

% give individual names to montecarlo parameters
[intdraws,draws,periods_simulations,seedzero,...
    gap_zero, b_gap_zero, b_pi_zero , gamma_t_zero, ...
    shutdownlearning] = translate_mcparameters(sim_parameters);

%% SHOCKS
% reset random number generator
rng(seedzero);

% create the cost push shock u_t
costpushshock = zeros(intdraws, periods_simulations);
if strcmp(shocktype,'irf')
    % IF IRF (one standard deviation shock)
    costpushshock(:,1)= sigeps;
elseif strcmp(shocktype,'series')
    costpushshock = sigeps.*randn(intdraws, periods_simulations);
end

%% ALLOCATE MEMORY FOR WELFARE 
welfare_cumul_MMS_final=zeros(1, periods_simulations);
welfare_x_cumul_final = zeros(1, periods_simulations);
welfare_pi_cumul_final = zeros(1, periods_simulations);


% when calculating the ergodic distribution for learning coefficients,
% draws is larger than 1 (it's a trick to generate very long series without 
% running into memory problems)
for jjj = 1:draws
    
    %% ALLOCATE MEMORY FOR SIMULATED SERIES
    [welfare_cumul_MMS, welfare_x_cumul, ...
            welfare_pi_cumul, welfare_inst, ...
            welfare_x_inst, welfare_pi_inst, gap_lag_simu, ...
            pi_simu, gap_simu, b_pi_simu, b_gap_simu, gamma_t_simu, lambda1_simu] = ...
            allocate_memory_for_series(intdraws,periods_simulations);
        
    %% Set initial conditions
     [b_pi_simu(:,1), b_gap_simu(:,1), gamma_t_simu(:,1), gap_lag_simu(:,1)] = ...
        set_initial_conditions(model, b_pi_comm, b_x_comm, ...
                b_pi_zero, b_gap_zero, gamma_t_zero, gap_zero);

    % Initiate welfare values at zero
    welfare_x = 0; welfare_pi = 0; welfare = 0;

    %% SIMULATION 
    for i=1: periods_simulations
        % set state variables 
        [states, gamma_t_simu(:,i)] = states_in_period_i(gap_lag_simu(:,i),...
            b_pi_simu(:,i),gamma_t_simu(:,i),gam,costpushshock(:,i),...
            gaintype,shutdownlearning);
        
        % set allocations
        [gap_simu(:,i), pi_simu(:,i), gap_lag_simu(:,i+1), ...
            b_pi_simu(:,i+1), b_gap_simu(:,i+1), ...
            lambda1_simu(:,i), gamma_t_simu(:,i+1)] = ...
            allocations_in_period_i(model, gaintype, states, fspace, parpolicy,...
            bench_p, b_gap_simu(:,i),sim_parameters);
        
        % calculate several welfare measures
        welfare = welfare + (betta.^(i-1)).*(pi_simu(:,i).^2 + ...
            alpha.*(gap_simu(:,i).^2));
        welfare_cumul_MMS(:,i) = welfare;
        
        welfare_x = welfare_x + (betta.^(i-1)).*(alpha.*(gap_simu(:,i).^2));
        welfare_pi = welfare_pi + (betta.^(i-1)).*(pi_simu(:,i).^2);
        welfare_x_cumul(:,i) = welfare_x ;
        welfare_pi_cumul(:,i) =welfare_pi;
        
        welfare_inst(:,i) = (pi_simu(:,i).^2 + alpha.*(gap_simu(:,i).^2));
        welfare_x_inst(:,i)  = alpha.*(gap_simu(:,i).^2);
        welfare_pi_inst(:,i) = pi_simu(:,i).^2 ;
        
        
    end
    
%     if ergodic==1
%         gap_zero = gap_simu(:,end);
%         b_pi_zero= b_pi_simu(:,end);
%         b_gap_zero= b_gap_simu(:,end);
%     end
    
    
    welfare_cumul_MMS_final = ((jjj-1).*welfare_cumul_MMS_final + ...
        mean(welfare_cumul_MMS,1))./jjj;
    welfare_x_cumul_final = ((jjj-1).*welfare_x_cumul_final + ...
        mean(welfare_x_cumul,1))./jjj;
    welfare_pi_cumul_final = ((jjj-1).*welfare_pi_cumul_final + ...
        mean(welfare_pi_cumul,1))./jjj;
    
end