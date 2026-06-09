function x = Get_x_from_R(R_target)
% GET_X_FROM_R Numerically calculates the x-coordinate for a given radius R.
% Restricted to the bow section of the SUBOFF hull where R(x) is one-to-one.

    RATIO = 0.2;
    R_MAX_SCALED = 5/6;
    
    % Calculate the actual maximum radius and max bow length in meters
    R_MAX_ACTUAL = R_MAX_SCALED * 0.3048 * RATIO; 
    x_max_bow = (10/3) * RATIO * 0.3048; 
    
    x = zeros(size(R_target));
    
    for i = 1:numel(R_target)
        Rt = R_target(i);
        
        % Check 1: Is the requested radius physically possible?
        if Rt < 0 || Rt > R_MAX_ACTUAL
            x(i) = NaN;
            warning('Radius %.4f is outside the bounds of the SUBOFF hull [0, %.4f].', Rt, R_MAX_ACTUAL);
            
        % Check 2: Are we at the absolute maximum radius?
        elseif abs(Rt - R_MAX_ACTUAL) < 1e-9
            x(i) = x_max_bow; % Return the start of the parallel midbody
            
        % Check 3: Are we at the very tip?
        elseif Rt == 0
            x(i) = 0;
            
        % Numerical Solver for the curved bow
        else
            % Define the objective function: f(x) = Get_R(x) - R_target = 0
            obj_fun = @(xi) Get_R(xi) - Rt;
            
            % Solve for x within the exact bounds of the bow section.
            % Bounding it guarantees rapid and flawless convergence.
            x(i) = fzero(obj_fun, [0, x_max_bow]);
        end
    end
end