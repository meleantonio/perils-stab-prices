function [fe_pi, fe_gap] = forecasterror(gap_next, pi_next, b_pi, b_gap, gap)

% FORECASTERROR calculates the forecast error for x_t and pi_t. Vectorized for speed and optimal memory allocation.

fe_pi = pi_next - b_pi.*gap;
fe_gap = gap_next - b_gap.*gap;

    
end