function equ = focs_dg(collocation_coefficients, Grid,fspace, ResWeights,parameters)

% FOCS_DG contains the first order conditions of the planner problem
% for decreasing gain model. It is used in the main file that solves for
% the coefficients collocation_coefficients that make equ as close to zero as possible
%
% INPUT:
%       collocation_coefficients: coefficients of the interpolants
%       Grid: the grid over which we solve the model
%       fspace: the functional space created with the Miranda-Fackler
%               Compecon Toolbox
%       ResWeights: this is non-zero if we use Galerkin method, it is zero
%       for collocation
%
% OUTPUT:
%       equ: a vector that contains the values of the two first order
%       conditions for every point of the grid, vertically stacked (i.e.
%       [equ1; equ2]

%% PARAMETERS
% give individual names to benchmark parameters
[alpha, betta, kappa, ~, ~, ~, ~, rho, ...
    ~, nQuadr, QuadrPoints, QuadrWeights, ~, ...
    ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ...
    ~, ~, ~, LowerBound, UpperBound,~] = translate_parameters(parameters);


% rename projection coefficients (just for convenience)
par0 = reshape(collocation_coefficients,length(collocation_coefficients)/2,2 );
parlambda1 = par0(:,1);
pargap = par0(:,2);

% evaluate current policy functions
lambda1 = funeval(parlambda1 , fspace, Grid)  ;
gap         = funeval(pargap , fspace, Grid)  ;

%rename grid
gap_lag = Grid(:,1);
b_pi    = Grid(:,2);
gamma_t    = Grid(:,3);
theta = Grid(:,4);

% generate inflation, next period's b^pi and decreasing gain parameter
pi = (betta.*b_pi +kappa).*gap + theta;
b_pi_next = b_pi  + gap_lag.*(pi - gap_lag.*b_pi).*gamma_t;
gamma_t_next = gamma_t./(gamma_t  + 1);

n = length(gap_lag);
% generate nQuadr replications of the Grid, one for each realization of shock:
Grid_gap = kron(gap,ones(nQuadr,1));
Grid_b_pi_next  = kron(b_pi_next ,ones(nQuadr,1));
Grid_gamma_t_next  = kron(gamma_t_next ,ones(nQuadr,1));

% Expected value of next theta, corresponding to Grid:
ExpTh = rho*theta;
% all realizations of next theta:
GridThNext = kron(ExpTh,ones(nQuadr,1)) + ...
    kron(ones(n,1),QuadrPoints);
% truncate it to state space:
GridThNext = min(max(GridThNext,LowerBound(4)),UpperBound(4));
GridNext = [Grid_gap Grid_b_pi_next Grid_gamma_t_next GridThNext];


% calculate variables at t+1
gap_next = funeval(pargap , fspace, GridNext) ;
lambda1_next =funeval(parlambda1 , fspace,GridNext) ;
pi_next = (betta.*Grid_b_pi_next + kappa).*gap_next+ GridThNext ;

% calculate expectations with quadrature
exp_pi_next_gap_next = (QuadrWeights'*reshape(pi_next.*gap_next,nQuadr,n))';
exp_lambda1_next = (QuadrWeights'*reshape(lambda1_next,nQuadr,n) )';
exp_lambda1_next_pi_next = (QuadrWeights'*reshape(lambda1_next.*pi_next ,nQuadr,n) )';
exp_lambda1_next_gap_next = (QuadrWeights'*reshape(lambda1_next.*gap_next ,nQuadr,n) )';


%%  EQUATIONS
equ1 = -alpha*gap - (betta.*b_pi + kappa).*pi - ...
    lambda1.*gamma_t.*gap_lag.*(betta.*b_pi + kappa) -...
    betta.*exp_lambda1_next_pi_next.*gamma_t_next +...
    2.*betta.*exp_lambda1_next.*b_pi_next.*gap.*gamma_t_next;

equ2 = lambda1 - betta.*exp_lambda1_next.*...
    (1-(gap.^2).*gamma_t_next) - ...
    (betta.^2).*(exp_pi_next_gap_next  + ...
    exp_lambda1_next_gap_next.*gamma_t_next.*gap);


% collect all equations in a vector
equ = [equ1; equ2] ;

% return weighted residuals if we use Galerkin method:
if(ResWeights~=0)
    equ = [ResWeights ResWeights ]*equ;
end

end
