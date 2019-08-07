function price = pricefrominflation(pi, periods_simulations)
% PRICEFROMINFLATION generates the series for p_t from inflation simulation data
price(:,1) = pi(:,1);
for i=2:periods_simulations
    price(:,i) = price(:,i-1) + pi(:,i);
end
end