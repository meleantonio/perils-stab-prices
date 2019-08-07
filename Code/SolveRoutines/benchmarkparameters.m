function p = benchmarkparameters(gaintype)

% BENCHMARKPARAMETERS generates benchmark parameters and the grid. 
% Input: 
%       gaintype:       either 'decr' or 'const'
%
% Output: 
%       p:              a structure with all parameters, including extrema for the grid

%% Benchmark parameters

p.alpha =.04;  % how much CB cares for output gap
p.betta =   .99; % discount factor
p.kappa =   .024; % parameter for the Phillips Curve
p.sig = .147; % for IS curve
p.gam =  .05; % gain parameter in case of constant gain (benchmark)
p.nphi = 10000^(1/4); % for testing the accuracy of the solution
p.sigma = 0.07; % standard deviation of the shock
p.rho = 0;% .95; % Persistence of the shock
p.omega = .47; % we take this parameter from Woodford (2003)

% Implementability parameter
p.imp = 0.001; % implementability constant
p.uncertainty_shock_sd = .01; % shock on beliefs for implementation

p.nQuadr = 150; %100;% %number of quadrature points;
% we choose nQuadr high to get smoothness;
[p.QuadrPoints,p.QuadrWeights] = qnwnorm(p.nQuadr,0,p.sigma^2);


% Define commitment and discretion solution under RE
p.b_x_comm = (p.kappa^2 + p.alpha*(1 + p.betta) - sqrt((p.kappa^2 + p.alpha*(1 + p.betta))^2 ...
    - 4*p.betta*(p.alpha^2)  ))/(2*p.betta*p.alpha);
p.b_pi_comm = - p.alpha*(p.b_x_comm-1)/p.kappa;
p.c_x_comm = - p.kappa*p.b_x_comm/p.alpha;
p.c_pi_comm = - p.alpha*p.c_x_comm/p.kappa;
p.c_x_discr = - p.kappa/(p.kappa^2 + p.alpha);
p.c_pi_discr = p.alpha/(p.kappa^2 + p.alpha);

%% Grid parameters for solving the model %

% Number of rounds of approximation
p.rounds_approx = 1; % 2; % 3;%

% Order of the polynomials, in a matrix for doing several rounds of
% approximation
if strcmp(gaintype,'decr')
    p.Order_vector = [[4;4;4;4]  [5;5;5;5]   [6;6;6;6]  ];
else
    p.Order_vector = [  [5;5;5]   [6;6;6]  ];
end

% Extrema for x_{t-1}
p.gap_lag_bar = 0;
p.gap_lag_sd =.04;
p.gap_lag_min = p.gap_lag_bar - p.gap_lag_sd;
p.gap_lag_max = p.gap_lag_bar +p.gap_lag_sd;

% Extrema for b^{pi}_{t-1}
p.b_pi_min =  -.0251; % 0; %-(p.b_pi_comm + .0251); % -0.02;% -b_pi_comm;%0; % -.251; %
p.b_pi_max = p.b_pi_comm + .0251;%

% Extrema for gamma_t
p.gamma_t_min = 0;%0.001;%  .9251; %
p.gamma_t_max = 1;

p.sigeps = p.sigma/sqrt(1-p.rho^2);
% Very large range for shock:
p.thmax = 3*p.sigeps;
p.thmin = -3*p.sigeps;


% Range on which we approximate the solution:
if strcmp(gaintype,'decr')
    p.LowerBound = [p.gap_lag_min; p.b_pi_min; p.gamma_t_min; p.thmin  ];
    p.UpperBound = [p.gap_lag_max; p.b_pi_max; p.gamma_t_max ;p.thmax   ];
    
elseif strcmp(gaintype,'const')
    p.LowerBound = [p.gap_lag_min; p.b_pi_min;  p.thmin  ];
    p.UpperBound = [p.gap_lag_max; p.b_pi_max; p.thmax   ];
end

end