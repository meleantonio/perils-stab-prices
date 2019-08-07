function [alpha, betta, kappa, sig, gam, nphi, sigmaa, rho, ...
            omega, nQuadr, QuadrPoints, QuadrWeights, b_x_comm, ...
            b_pi_comm, c_x_comm, c_pi_comm, c_x_discr, ...
            c_pi_discr, rounds_approx, Order_vector, ...
            gap_lag_bar, gap_lag_sd, gap_lag_min, gap_lag_max, ...
            b_pi_min, b_pi_max, gamma_t_min, gamma_t_max, ...
            sigeps, thmax, thmin, LowerBound, UpperBound, imp, uncertainty_shock_sd] = translate_parameters(p)

% This function translates a parameters' structure into single parameters for the montecarlo simulations. The parameters' structure has the following values

alpha = p.alpha;
betta = p.betta;
kappa = p.kappa;
sig = p.sig; 
gam = p.gam; 
nphi = p.nphi;
sigmaa = p.sigma;
rho = p.rho;
omega = p.omega;
nQuadr = p.nQuadr;
QuadrPoints = p.QuadrPoints;
QuadrWeights = p.QuadrWeights;
b_x_comm = p.b_x_comm;
b_pi_comm = p.b_pi_comm;
c_x_comm = p.c_x_comm;
c_pi_comm = p.c_pi_comm;
c_x_discr = p.c_x_discr;
c_pi_discr = p.c_pi_discr;
rounds_approx = p.rounds_approx;
Order_vector = p.Order_vector;
gap_lag_bar = p.gap_lag_bar;
gap_lag_sd = p.gap_lag_sd;
gap_lag_min = p.gap_lag_min;
gap_lag_max = p.gap_lag_max;
b_pi_min = p.b_pi_min;
b_pi_max = p.b_pi_max;
gamma_t_min = p.gamma_t_min;
gamma_t_max = p.gamma_t_max;
sigeps = p.sigeps;
thmax = p.thmax;
thmin = p.thmin;
LowerBound = p.LowerBound;
UpperBound = p.UpperBound;
imp = p.imp;
uncertainty_shock_sd = p.uncertainty_shock_sd;

end
