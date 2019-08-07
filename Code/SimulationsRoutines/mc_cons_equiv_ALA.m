function mcce = mc_cons_equiv_ALA(parpolicy, fspace, bench_p, model, ...
            gaintype, shocktype, sim_parameters, OLS_flag, ALA)

% MC_CONS_EQUIV_ALA simulates the model at hand and returns the consumption equivalent welfare measure (Alternative learning algorithms)

[~ , ~  , ~ ,~ ,~ , ~ ,~ , welfare_cumul , ~ ,  ~ , ~ , ~  , ...
            ~ , ~ , ~ , ~ ,  ~  ] = simul_fast_ALA(parpolicy, fspace, bench_p, model, ...
                gaintype, shocktype, sim_parameters, OLS_flag, ALA);
            
welfare_losses = mean(welfare_cumul,1);
% Calculate consumption equivalent measure
mcce  = consumption_equivalent(welfare_losses(:,end), bench_p);

end