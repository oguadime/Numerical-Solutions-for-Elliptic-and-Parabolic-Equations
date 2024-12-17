# Numerical Solutions for Elliptic and Parabolic Equations
This repository contains MATLAB implementations of numerical methods for solving elliptic and parabolic partial differential equations (PDEs). Methods include fixed-point iteration, fully implicit schemes, and IMEX (Implicit-Explicit) methods, with capabilities for adaptive time-stepping and error estimation.

# Decsription 
This repository includes MATLAB scripts that implement various numerical schemes for solving nonlinear PDEs and systems of ODEs. Key highlights of the repository include:
- Elliptic Equation Solvers: ELLIPTICFixedPoint1d.m; Solves 1D elliptic equations using a fixed-point iteration method with finite volume discretization​. ELLIPTIC1db.m; Demonstrates fixed-point solutions for 1D elliptic problems with nonlinear source terms and boundary conditions​.
- Fully Implicit Systems: Implements Newton's method for solving coupled nonlinear ODE systems with options for fixed and adaptive time-stepping. Scripts: FullyImplicitSystem1a.m, FullyImplicitSystem1b.m, FullyImplicitSystem1c.m, FullyImplicitSystem2a.m, FullyImplicitSystem2b.m​
- IMEX Scheme (IMEX.m): Combines implicit and explicit methods for solving parabolic equations with nonlinear source terms.


# Key Features
- Elliptic Solvers: Supports uniform and non-uniform grids. Handles nonlinear source terms and Dirichlet/Neumann boundary conditions. Implements finite volume discretization with transmissibilities.
- Fully Implicit Systems: Solves coupled nonlinear systems using Newton's method. Tracks solver performance metrics (e.g., average and maximum Newton iterations).
- IMEX Method: Combines implicit handling of diffusion terms with explicit treatment of nonlinearities for stability and efficiency.

# Usage
Each script is modular and can be executed directly in MATLAB. Key parameters (e.g., time step, grid resolution) can be adjusted via input arguments.
Example 1: Fixed-Point Solver for Elliptic Equations
```matlab
tau = 0.01;
eps = 0.1;
a_ = 1;
nxdx = 10; % Number of grid points
bcond = [0, 0, 0, 1]; % Dirichlet boundary conditions
ELLIPTICFixedPoint1d(tau, eps, a_, nxdx, bcond, 1, 1, 0);
```

Example 2: Fully Implicit System with Newton's Method
```matlab
alpha = 5;
c = 15;
tau = 0.025;
FullyImplicitSystem1a(alpha, c, tau);
```

Example 3: IMEX Scheme 
```matlab
IMEX(0, 1, 1, 1, 20, 1, 0.01);
```

## License
This project is licensed under the MIT License - see the LICENSE file for details.
```
Feel free to adjust any part of this README to better fit your specific needs or preferences.
