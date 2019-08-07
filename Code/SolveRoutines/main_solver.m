function [parpolicy, fspace, Grid,  max_test] = main_solver(parameters,projection_parameters, gaintype)

% MAIN_SOLVER solves the planner problem for a specific version of the learning (decreasing or constant)
%
% Inputs: 
% parameters: the model's parameters
% projection_parameters: parameters for the collocation method
% gaintype: string indicating which type of learning ('decr' or 'const')
%
% Outputs:
% parpolicy: coefficients of the optimal interpolants saved in a vector
% fspace: functional space for the basis functions
% Grid: the grid over which we solve the problem
% max_test: accuracy statistics (the max of the residuals' errors in a larger grid)


%% PARAMETERS
% give individual names to benchmark parameters
[~, ~, ~, ~, ~, nphi, ~, ~, ...
    ~, ~, ~, ~, ~, ...
    ~, ~, ~, ~, ...
    ~, rounds_approx, Order_vector, ...
    ~, ~, ~, ~, ...
    ~, ~, ~, ~, ...
    ~, ~, ~, LowerBound, UpperBound,~] = translate_parameters(parameters);

% give indivisual names to projection parameters
[approxtype, splineorder] = translate_projectionparameters(projection_parameters);


tic;
% we solve iteratively over finer grids and higher orders of approximations
for oo = 1: rounds_approx
    RoundAppr = oo;     % round of approximation is stored for convenience
    disp(sprintf('  RoundAppr = %d ',RoundAppr));
    
    % take the order and number of gridpoints from Order_vector
    Order = Order_vector(:,oo);
    
    % save the functional space from previous iteration
    % (this is used to generate a new guess in the new functional space)
    if oo >= 2
        fspace_old = fspace;
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % generate the space of basis functions
    % by using CompEcon toolbox for function
    % approximation (see Miranda-Fackler book)
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % the following lines generate basis function space
    % we can choose among chebychev polynomials, splines of different
    % orders and piecewise linear functions. Baseline code uses Chebychev
    % polynomials
    
    
    %     approxtype = 'lin';  % piecewise linear
    % approxtype = 'cheb'; % chebychev polynomials
    %         approxtype = 'spli'; % splines
    % splineorder = []; % splines' order, default are cubic splines
    if(strcmp(approxtype,'spli'))
        fspace = fundefn(approxtype,Order,LowerBound,UpperBound,splineorder);
        nodes = funnode(fspace);
    else
        fspace = fundefn(approxtype,Order,LowerBound,UpperBound,[]);
        nodes = funnode(fspace);
    end
    
    % create gridpoints
    Grid = gridmake(nodes);
    
    % Initial conditions for the nonlinear equation solver
    if (RoundAppr == 1) 
        if strcmp(gaintype, 'decr')
            lambda1 = .1+ zeros(length(Grid),1);
            gap        =.1+ zeros(length(Grid),1);
        else
            %% for first round, initial conditions from an external file
            load init_cond parlambda2 pargap fspace_init;
            lambda1 =funeval(parlambda2 , fspace_init, Grid)  ;
            gap         =funeval(pargap , fspace_init, Grid)  ;
        end        
    else    % if we are at second approx round, we use the solution of the first round
        % as initial conditions on the new larger grid
        lambda1 =funeval(parlambda1 , fspace_old, Grid)  ;
        gap         =funeval(pargap , fspace_old, Grid)  ;
    end;
    
    
    % generate basis functions at Grid and residual weights (non-zero only
    % if we use Galerking method):
    Basis = funbas(fspace,Grid);
    ResWeights = 0;
    
    % set initial value for coefficients of the interpolation
    parlambda1 =  Basis\lambda1;
    pargap         = Basis\gap;
    
    % save all coefficients in the vector parpolicy
    parpolicy = [parlambda1; pargap];
    %%              M A I N   L O O P       %%
    if strcmp(gaintype,'decr')
        % solve Lagrangean FOCs by Broyden method for nonlinear equations
        [opt_vec,info] =  broydn('focs_dg',parpolicy,1e-8,0,1,Grid,...
            fspace, ResWeights, parameters);
    else 
        % solve Lagrangean FOCs by Broyden method for nonlinear equations
        [opt_vec,info] =  broydn('focs_cg',parpolicy,1e-8,0,1,Grid,...
            fspace, ResWeights, parameters);
    end

    disp(sprintf(' info = %d',info)); % if info=0, everything went fine,
    % o/w the Broyden algorithm didn't converge
    disp(sprintf('    '));
    
    %rename optimized parameters
    par00 = opt_vec;
    par0 = reshape(par00,length(Grid),2 );
    parlambda1 =                par0(:,1);
    pargap =      par0(:,2);
    
    % update vector of coefficients
    parpolicy = opt_vec;
    
end

toc;

% Calculate time needed for solution
time_hours = fix(toc/3600);
time_minutes = fix((toc - time_hours*3600)/60);
time_seconds = fix(toc -  time_hours*3600 - time_minutes*60);

disp(sprintf('    '));disp(sprintf('Time for solving the model was    '));
disp(sprintf('%d hours     %d  minutes %d  seconds',time_hours,time_minutes,time_seconds));

%%     T E S T I N G    %%
% Testing the approximation: compute residuals at points off the grid:

% Create a large grid with nphi points
grid_nodes_test1 = linspace(LowerBound(1),UpperBound(1),nphi)';%
grid_nodes_test2 = linspace(LowerBound(2),UpperBound(2),nphi)';%
grid_nodes_test3 = linspace(LowerBound(3),UpperBound(3),nphi)';%
if strcmp(gaintype,'decr')
    grid_nodes_test4 = linspace(LowerBound(4),UpperBound(4),nphi)';%
    GridTest = gridmake(grid_nodes_test1,grid_nodes_test2,grid_nodes_test3,grid_nodes_test4);
    % calculate residuals on the test grid with optimized parameters
    test_residuals =  focs_dg(parpolicy,GridTest,fspace,ResWeights, parameters);
else
    GridTest = gridmake(grid_nodes_test1,grid_nodes_test2,grid_nodes_test3);
    %calculate residuals on the new grid with optimized parameters
    test_residuals =  focs_cg(parpolicy,GridTest,fspace,ResWeights,parameters);
end

% generate statistics to measure the accuracy
norm_test  = norm(test_residuals);
max_test   = max(abs(test_residuals));

% disp(sprintf(' max_test = %d',testing_max));
disp(sprintf(' max_test = %d',max_test));
disp(sprintf('    '));
