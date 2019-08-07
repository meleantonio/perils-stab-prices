* ```p.alpha  ```: how much CB cares for output gap
* ```p.betta ```: discount factor
* ```p.kappa ```: parameter for the Phillips Curve
* ```p.sig ```: for IS curve
* ```p.gam .05; ```: gain parameter in case of constant gain (benchmark)
* ```p.nphi ```: number of gridpoints for testing the accuracy of the solution
* ```p.sigma ```: standard deviation of the shock
* ```p.rho = ```: Persistence of the shock
* ```p.omega ```: we take this parameter from Woodford (2003)
* ```p.nQuadr ```:number of quadrature points;
* ```p.QuadrPoints ```: parameter for calculating the expectations with quadrature
* ```p.QuadrWeights ```: parameter for calculating the expectations with quadrature
* ```p.b_x_comm ```: RE commitment b^x
* ```p.b_pi_comm ```: RE commitment b^pi
* ```p.c_x_comm ```: RE commitment c^gap
* ```p.c_pi_comm ```: RE commitment c^pi
* ```p.c_x_discr ```: RE discretion b^gap
* ```p.c_pi_discr ```: RE discretion b^pi
* ```p.rounds_approx ```: Number of rounds of approximation for the collocation
* ```p.Order_vector ```: Order of the polynomials (for Chebichev or splines) in a matrix (each column corresponds to a round of approximation)
* Extrema for x_{t-1}
    * ```p.gap_lag_min```
    * ```p.gap_lag_max```

* Extrema for b^{pi}_{t-1}
    * ```p.b_pi_min```
    * ```p.b_pi_max```

* Extrema for gamma_t
    * ```p.gamma_t_min```
    * ```p.gamma_t_max```

* ```p.sigeps = p.sigma/sqrt(1-p.rho^2);```: variance corrected for persistence
* Range for shock:
    * ```p.thmax```
    * ```p.thmin```

* Range on which we approximate the solution:
    * ```p.LowerBound```
    * ```p.UpperBound```
