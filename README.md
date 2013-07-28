black-scholes
=============

Solve Black Scholes using Crank-Nicolson Finite Difference method.

[![BlackScholesEquation](http://upload.wikimedia.org/math/5/e/f/5ef2fa747d3a5d684ae67bdc7236e6d4.png)](http://en.wikipedia.org/wiki/Black%E2%80%93Scholes)

This code numerically solves hyperbolic PDEs of the form: 
    
    Dt[u] + a Dx[u] + b Dy[u] + b Dxx[u] + u = F(t, x)
    
    where Dt[], Dx[], Dy[], and Dxx[] are the differential operators for t, x, and y

The solutions are animated in a window. Data is saved to text files when 's' is pressed.

Schemes
-------

The following explicit finite difference schemes are implemented:

* Forward-Time Forward-Space
* Forward-Time Backward-Space
* Forward-Time Central-Space
* Lax Fredrichs
* Leapfrog
* Equilibrium
* Lax Wendroff

The following implicit finite difference schemes are implemented:

*	Crank Nicolson
