function p = mcparameters(i,d,ps, seed0,x0,bg0, bp0,g0,sdl)

% MCPARAMETERS generates the Montecarlo parameters 
% Input: 
%       i:      number of draws 
%       d:      number of repetitions of the draws
%       ps:     number of periods for the simulation
%       seed:   random number generator seed
%       x0:     output gap initial condition
%       bg0:    b_gap initial condition
%       bp0:    b_pi initial condition
%       g0:     gamma initial condition
%       sdl:    shutdownlearning
%
% Output: 
%       p:              a structure with all parameters


% montecarlo parameters
p.intdraws = i;
p.draws = d;
p.periods_simulations = ps;
p.seedzero = seed0;
p.gap_zero = x0;
p.b_gap_zero = bg0;
p.b_pi_zero = bp0; 
p.gamma_t_zero = g0; 
p.shutdownlearning = sdl; 

end