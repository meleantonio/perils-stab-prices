function [mcce_MMS_from_PLT, mcce_EHCOMM_from_PLT, mcce_EHDISCR_from_IT] = ...
    mc_simul_ALA(parpolicy, fspace, bench_p, model, ...
        gaintype, shocktype, sim_parameters, OLS_flag, ALA, sim_parameters_IT, ...
        periods_simulations_mc, gamma_zero, ~)    

% give individual names to montecarlo parameters
[intdraws,draws,periods_simulations,seedzero,...
    gap_zero, ~, ~ , ~, ...
    shutdownlearning] = translate_mcparameters(sim_parameters);


% create a Montecarlo simulation for each model
[~,~ , ~,b_pi_simu_MMS,b_gap_simu_MMS,~,~,~,~,~,~,~,~,~,~,~,~,~] = ...
    simul_fast_ALA(parpolicy, fspace, bench_p, model, ...
    gaintype, shocktype, sim_parameters, OLS_flag, ALA);
[~,~ , ~,b_pi_simu_EH_PLT,b_gap_simu_EH_PLT,~,~,~,~,~,~,~,~,~,~,~,~,~] = ...
    simul_fast_ALA(parpolicy, fspace, bench_p, 'EHCOMM', ...
    gaintype, shocktype, sim_parameters, OLS_flag, ALA);
[~,~ , ~,b_pi_simu_EH_IT,b_gap_simu_EH_IT,~,~,~,~,~,~,~,~,~,~,~,~,~] = ...
    simul_fast_ALA(parpolicy, fspace, bench_p, 'EHDISCR', ...
    gaintype, shocktype, sim_parameters_IT, OLS_flag, ALA);


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
    mcce_MMS_from_PLT(:,i) = mc_cons_equiv_ALA(parpolicy, fspace, bench_p, model, ...
            gaintype, shocktype, sim_parameters_MMS, OLS_flag, ALA);
    % EVANS-HONKAPOHJA COMMITMENT
    disp('Simulating EH PLT-inducing policy from PLT beliefs...')
    mcce_EHCOMM_from_PLT(:,i) = mc_cons_equiv_ALA(parpolicy, fspace, bench_p, 'EHCOMM', ...
            gaintype, shocktype, sim_parameters_EH_PLT, OLS_flag, ALA);
    % EVANS-HONKAPOHJA DISCRETION       
    disp('Simulating EH IT-inducing policy from IT beliefs...') 
    mcce_EHDISCR_from_IT(:,i) = mc_cons_equiv_ALA(parpolicy, fspace, bench_p, 'EHDISCR', ...
            gaintype, shocktype, sim_parameters_EH_IT, OLS_flag, ALA);

end
