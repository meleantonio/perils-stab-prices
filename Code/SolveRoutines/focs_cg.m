function equ = focs_cg(collocation_coefficients, Grid,fspace, ResWeights, parameters)

% FOCS_CG contains the first order conditions of the planner problem
% for the constant gain model. It is used in the main file that solves for
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
[alpha, betta, kappa, ~, gam, ~, ~, rho, ...
    ~, nQuadr, QuadrPoints, QuadrWeights, ~, ...
    ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ...
    ~, ~, ~, LowerBound, UpperBound, ~] = translate_parameters(parameters);


% rename projection coefficients (just for convenience)
par0 = reshape(collocation_coefficients,length(collocation_coefficients)/2,2 );
parlambda1 = par0(:,1);
pargap = par0(:,2);

% evaluate policy functions
lambda1 = funeval(parlambda1 , fspace, Grid)  ;
gap         = funeval(pargap , fspace, Grid)  ;

%rename grid
gap_lag = Grid(:,1);
b_pi    = Grid(:,2);
theta = Grid(:,3);

% generate inflation and next period's b^pi
pi = (betta.*b_pi +kappa).*gap + theta;
b_pi_next = b_pi  + gam.*gap_lag.*(pi - gap_lag.*b_pi);


n = length(gap_lag);
% generate nQuadr replications of the Grid, one for each realization of shock:
Grid_gap = kron(gap,ones(nQuadr,1));
Grid_b_pi_next  = kron(b_pi_next ,ones(nQuadr,1));

% Expected value of next theta, corresponding to Grid:
ExpTh = rho*theta;
% all realizations of next theta:
GridThNext = kron(ExpTh,ones(nQuadr,1)) + ...
    kron(ones(n,1),QuadrPoints);
% truncate it to state space:
GridThNext = min(max(GridThNext,LowerBound(3)),UpperBound(3));
GridNext = [Grid_gap Grid_b_pi_next GridThNext];


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
    lambda1.*gam.*gap_lag.*(betta.*b_pi + kappa) -...
    betta.*exp_lambda1_next_pi_next.*gam +...
    2.*betta.*exp_lambda1_next.*b_pi_next.*gap.*gam;

equ2 = lambda1 - betta.*exp_lambda1_next.*...
    (1-(gap.^2).*gam) - ...
    (betta.^2).*(exp_pi_next_gap_next  + ...
    exp_lambda1_next_gap_next.*gam.*gap);


% collect all equations in a vector
equ = [equ1; equ2] ;

% return weighted residuals if we use Galerkin method:
if(ResWeights~=0)
    equ = [ResWeights ResWeights ]*equ;
end

end
