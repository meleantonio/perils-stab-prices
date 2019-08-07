function [intdraws,draws,periods_simulations,seedzero,...
    gap_zero, b_gap_zero, b_pi_zero , gamma_t_zero, ...
    shutdownlearning] = translate_mcparameters(p)

% TRANSLATE_MCPARAMETERS translates the Montecarlo parameters from a structure to individual variables
%

% montecarlo parameters
intdraws  = p.intdraws;
draws = p.draws;
periods_simulations = p.periods_simulations;
seedzero = p.seedzero;
gap_zero = p.gap_zero;
b_gap_zero = p.b_gap_zero ;
b_pi_zero = p.b_pi_zero ; 
gamma_t_zero = p.gamma_t_zero ; 
shutdownlearning = p.shutdownlearning; 

end