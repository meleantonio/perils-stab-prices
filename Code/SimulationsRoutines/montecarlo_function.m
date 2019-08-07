
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% NOTE: switching the montecarlo codes from having separated scripts to having a montecarlo function 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [outargs] = montecarlo_function(gaintype, ...
                        shutdownlearning, shocktype, model, ...
                        intdraws, draws,periods_simulations, ergodic, ...
                        gap_zero, b_gap_zero, b_pi_zero, gamma_t_zero, parameters )


% MONTECARLO_FUNCTION produces a Montecarlo simulation for the model
% and produces a graph of the ergodic distribution of the coefficient
% b^{\pi}
%
% inputs:
%     gaintype: can be set to 'decr' for decreasing gain, or 'const'
%         for constant gain
%     shutdownlearning:           , 
%     shocktype: can be set to 'irf' for impulse response function or
%         'series' for general montecarlo simulation

%       shocktype:                  can be 'irf', model, intdraws, draws,periods_simulations, ergodic, gap_zero, b_gap_zero, b_pi_zero, gamma_t_zero,parameters 
%
% outputs: 



[alpha, betta, kappa, sig, gam, nphi, sigma, rho, ...
            omega, nQuadr, QuadrPoints, QuadrWeights, b_x_comm, ...
            b_pi_comm, c_x_comm, c_pi_comm, c_x_discr, ...
            c_pi_discr, rounds_approx, Order_vector, ...
            gap_lag_bar, gap_lag_sd, gap_lag_min, gap_lag_max, ...
            b_pi_min, b_pi_max, gamma_t_min, gamma_t_max, ...
            sigeps, thmax, thmin, LowerBound, UpperBound, imp] = translate_parameters(parameters)


distribution_bi_pi = zeros(intdraws,length(gam_vector));
distribution_bi_gap = zeros(intdraws,length(gam_vector));

distribution_bi_pi_zero = zeros(intdraws,length(gam_vector));
distribution_bi_gap_zero = zeros(intdraws,length(gam_vector));

% this is for the simulation file simul_MMS_2eqs.m: 1 is for the ergodic distribution, 0
% otherwise
ergodic =1;


% solve the model first
if strcmp(gaintype,'const')
    main_cg;
elseif strcmp(gaintype,'decr')
    main_dg;
end
    
    
    
    % Initial conditions for the simulation
    gap_zero = 0;
    b_gap_zero = -b_x_comm + 2*b_x_comm.*rand(intdraws,1) ;%-.02 + (b_x_comm+.02).*rand(intdraws,1) ;%0;%.25*0;%0;% .5*.5*
    b_pi_zero =  b_pi_min + (b_pi_max-b_pi_min).*rand(intdraws,1) ;%-.02 + (b_pi_comm+.02).*rand(intdraws,1) ;%0;%.25*  0;%0;%.5*
    gamma_t_zero = 1;
    
    
    % store the initial distributions
    distribution_bi_pi_zero(:,k) = b_pi_zero;
    distribution_bi_gap_zero(:,k) = b_gap_zero;
    
    
    % simulation

[gap_lag_simu, pi_simu , gap_simu,b_pi_simu,b_gap_simu, gamma_t_simu,lambda1_simu,...
            welfare_cumul_MMS, welfare_x_cumul, welfare_pi_cumul,...
            welfare_cumul_MMS_final , welfare_x_cumul_final , welfare_pi_cumul_final, ...
            costpushshock, welfare_inst , welfare_x_inst, welfare_pi_inst ] = ...
            simul(parpolicy, fspace,betta, alpha, kappa, model, ...
                gaintype, shocktype, gap_zero, b_gap_zero, b_pi_zero, ...
                gamma_t_zero, shutdownlearning)


%     [gap_lag_simu, pi_simu , gap_simu,b_pi_simu,b_gap_simu,...
%         gamma_t_simu,lambda1_simu,...
%         welfare_cumul_MMS, welfare_x_cumul,welfare_pi_cumul,...
%         welfare_cumul_MMS_final ,welfare_x_cumul_final , ...
%         welfare_pi_cumul_final,costpushshock ] = simul_MMS_2eqs(parpolicy, fspace,betta, ...
%         alpha, kappa, model, ...
%         gaintype, shocktype, ...
%         gap_zero,b_gap_zero,b_pi_zero,gamma_t_zero,shutdownlearning);
    
    % store the distributions
%    distribution_bi_pi(:,k) = b_pi_simu(:,end);
%    distribution_bi_gap(:,k) = b_gap_simu(:,end);
    
end

% 
% figure(1);
% subplot(2,3,1);
% xax = -.2:0.001:.2;
% mm = length(xax);
% [n,xout] = hist(distribution_bi_pi(:,1),xax);
% bar(xout,n./intdraws);
% title('\gamma = .01');
% subplot(2,3,2);
% [n,xout] = hist(distribution_bi_pi(:,2),xax);
% bar(xout,n./intdraws);
% title('\gamma = .02');
% subplot(2,3,3);
% [n,xout] = hist(distribution_bi_pi(:,3),xax);
% bar(xout,n./intdraws);
% title('\gamma = .05');
% subplot(2,3,4);
% [n,xout] = hist(distribution_bi_pi(:,4),xax);
% bar(xout,n./intdraws);
% title('\gamma = .1');
% subplot(2,3,5);
% [n,xout] = hist(distribution_bi_pi(:,5),xax);
% bar(xout,n./intdraws);
% title('\gamma = .2');
% subplot(2,3,6);
% [n,xout] = hist(distribution_bi_pi(:,6),xax);
% bar(xout,n./intdraws);
% title('\gamma = .4');
% 
% ,'EH');
% 
% 
% 
% figure(2)
% subplot(2,3,1)
% plot(distribution_bi_pi_zero(:,1),distribution_bi_pi(:,1),'.');
% title('\gamma = .01');
% 
% subplot(2,3,2)
% plot(distribution_bi_pi_zero(:,2),distribution_bi_pi(:,2),'.');
% title('\gamma = .02');
% 
% subplot(2,3,3)
% plot(distribution_bi_pi_zero(:,3),distribution_bi_pi(:,3),'.');
% title('\gamma = .05');
% 
% subplot(2,3,4)
% plot(distribution_bi_pi_zero(:,4),distribution_bi_pi(:,4),'.');
% title('\gamma = .1');
% 
% subplot(2,3,5)
% plot(distribution_bi_pi_zero(:,5),distribution_bi_pi(:,5),'.');
% title('\gamma = .2');


