function [kappa0, kappa0_s] = calculate_suboff_curvature(x)
% CALCULATE_SUBOFF_CURVATURE Computes axial curvature and its derivative
% for the DARPA SUBOFF model head section based on Appendix F formulas.
%
% Inputs:
%   x - Vector of axial coordinates
%
% Outputs:
%   kappa0   - Axial curvature \kappa_0
%   kappa0_s - Derivative of axial curvature with respect to arc length d\kappa_0/ds
    RATIO=0.2;
    x=x/RATIO/0.3048;
    % Initialize outputs
    kappa0   = zeros(size(x));
    kappa0_s = zeros(size(x));

    % SUBOFF Constants
    A     = 1.126395101;
    B     = 0.442874707;
    Rmax  = 5/6;
    alpha = 1/2.1;

    % The head section is defined for 0 <= x <= 3.333333.
    % The parallel middle body (x > 3.333333) is a cylinder, so curvature is 0.
    head_idx = (x >= 0) & (x <= 3.333333);
    xh = x(head_idx);

    % Robustness: Prevent singularity at x = 0 where U(0) = 0 and U^(alpha-1) -> Inf
    xh(xh == 0) = eps;

    % Auxiliary variable z to simplify equations and prevent typos
    z = 0.3 * xh - 1;

    % 1. Calculate Auxiliary Functions U, V, W, J (F.2 - F.5)
    U = A .* xh .* z.^4 + B .* xh.^2 .* z.^3 + 1 - z.^4 .* (1.2 * xh + 1);
    
    V = 1.2 * A .* xh .* z.^3 + A .* z.^4 + 0.9 * B .* xh.^2 .* z.^2 ...
        + 2 * B .* xh .* z.^3 - 1.2 .* z.^3 .* (1.2 * xh + 1) - 1.2 .* z.^4;
        
    W = (2.4 * A + 2 * B - 2.88) .* z.^3 + (1.08 * A + 3.6 * B) .* xh .* z.^2 ...
        + 0.54 * B .* xh.^2 .* z - 1.08 .* z.^2 .* (1.2 * xh + 1);
        
    J = (3.24 * A + 5.4 * B - 3.888) .* z.^2 + (0.648 * A + 3.24 * B) .* xh .* z ...
        + 0.162 * B .* xh.^2 - 0.648 .* z .* (1.2 * xh + 1);

    % 2. Calculate derivatives of r0 with respect to x (F.7 - F.9)
    r0_x   = Rmax * alpha .* U.^(alpha - 1) .* V;
    
    r0_xx  = Rmax * alpha * (alpha - 1) .* U.^(alpha - 2) .* V.^2 ...
             + Rmax * alpha .* U.^(alpha - 1) .* W;
             
    r0_xxx = Rmax * alpha * (alpha - 1) * (alpha - 2) .* U.^(alpha - 3) .* V.^3 ...
             + 3 * Rmax * alpha * (alpha - 1) .* U.^(alpha - 2) .* V .* W ...
             + Rmax * alpha .* U.^(alpha - 1) .* J;

    % 3. Calculate Trignometric relations (F.10 - F.12)
    tan_phi = r0_x;
    cos_phi = (1 + tan_phi.^2).^(-1/2);
    % Using abs to ensure no complex numbers arise from floating-point inaccuracy
    sin_phi = sqrt(abs(1 - cos_phi.^2)); 

    % 4. Calculate derivatives of inclination angle phi with respect to x (F.13 - F.14)
    phi_x  = (1 ./ (1 + r0_x.^2)) .* r0_xx;
    
    phi_xx = -(2 ./ (1 + r0_x.^2).^2) .* r0_x .* r0_xx.^2 ...
             + (1 ./ (1 + r0_x.^2)) .* r0_xxx;

    % 5. Calculate curvature and its arc length derivative (F.15 - F.16)
    k0   = -phi_x .* cos_phi;
    k0_s = -phi_xx .* cos_phi.^2 + (phi_x.^2) .* sin_phi .* cos_phi;

    % Assign calculated head values back to the output vectors
    kappa0(head_idx)   = k0;
    kappa0_s(head_idx) = k0_s;
    
    kappa0 = kappa0/0.3048/RATIO;
    kappa0_s = kappa0_s/0.3048^2/RATIO^2;
end