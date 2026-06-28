function x = Get_x_from_S(S_target)
% GET_X_FROM_S Calculates the axial coordinate x for a given arc length S.
% Intelligently switches between numerical root-finding (for the bow) 
% and exact algebraic calculation (for the parallel mid-body).

    RATIO = 0.2;
    x_max_bow = (10/3) * RATIO * 0.3048;
    
    % Step 1: Pre-calculate the exact arc length where the curved bow ends.
    % We use our robust integral solver to find this exact boundary point.
    S_bow_max = Get_S(x_max_bow);
    
    x = zeros(size(S_target));
    
    for i = 1:numel(S_target)
        St = S_target(i);
        
        % Boundary Check
        if St < 0
            x(i) = NaN;
            warning('Arc length S cannot be negative.');
            
        elseif St == 0
            x(i) = 0;
            
        % Region 1: The Curved Bow Section (Numerical Solver)
        elseif St <= S_bow_max
            % Define objective function: f(x) = S(x) - S_target = 0
            obj_fun = @(xi) Get_S(xi) - St;
            x(i) = fzero(obj_fun, [0, x_max_bow]);
            
        % Region 2: The Parallel Middle Body (Analytical Exact Solution)
        else
            % Because slope K=0 here, arc length increases 1:1 with x.
            % We just take the max bow coordinate and add the remaining arc length.
            x(i) = x_max_bow + (St - S_bow_max);
        end
    end
end