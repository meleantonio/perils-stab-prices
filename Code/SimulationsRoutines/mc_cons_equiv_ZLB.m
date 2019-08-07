function mcce = mc_cons_equiv_ZLB(parpolicy, fspace, bench_p, model, ...
            gaintype, shocktype, sim_parameters)

% MC_CONS_EQUIV_ZLB simulates the model at hand under the ZLB constraint and returns the consumption equivalent welfare measure

[~ , ~  , ~ ,~ ,~ , ~ ,~ , welfare_cumul , ~ ,  ~ , ~ , ~  , ...
            ~ , ~ , ~ , ~ ,  ~  ] = simul_fast_ZLB(parpolicy, fspace, bench_p, model, ...
                gaintype, shocktype, sim_parameters);
            
welfare_losses = mean(welfare_cumul,1);
% Calculate consumption equivalent measure
mcce  = consumption_equivalent(welfare_losses(:,end), bench_p);

end