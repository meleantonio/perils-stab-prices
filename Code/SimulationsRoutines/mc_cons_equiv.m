function mcce = mc_cons_equiv(parpolicy, fspace, bench_p, model, ...
            gaintype, shocktype, sim_parameters)

% MC_CONS_EQUIV simulates the model at hand and returns the consumption equivalent welfare measure

[~ , ~  , ~ ,~ ,~ , ~ ,~ , welfare_cumul , ~ ,  ~ , ~ , ~  , ...
            ~ , ~ , ~ , ~ ,  ~  ] = simul_fast(parpolicy, fspace, bench_p, model, ...
                gaintype, shocktype, sim_parameters);
            
welfare_losses = mean(welfare_cumul,1);
% Calculate consumption equivalent measure
mcce  = consumption_equivalent(welfare_losses(:,end), bench_p);

end