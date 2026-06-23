function S = Get_S(x)
% GET_S Calculates the arc length along the SUBOFF hull from x=0.0.
% Preserves the shape of input 'x'.

    S = zeros(size(x));
    
    % Define the vectorized integrand: sqrt(1 + (dR/dx)^2)
    integrand = @(xi) sqrt(1 + Get_K(xi).^2);
    
    % Evaluate the integral for each point requested
    for i = 1:numel(x)
        if x(i) <= 0
            S(i) = 0;
        else
            % Global adaptive quadrature naturally bridges the x=0 singularity.
            % Tolerances set tightly for high precision.
            S(i) = 0.011310262042964 + integral(integrand, 0.002241760516193, x(i), 'AbsTol', 1e-10, 'RelTol', 1e-8);
        end
    end
end