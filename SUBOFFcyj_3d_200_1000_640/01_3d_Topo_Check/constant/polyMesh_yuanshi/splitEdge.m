function x = splitEdge(L, CN, lenRatios, cellRatios, expRatios)
    % Mimics OpenFOAM's simpleGrading multi-grading
    % 
    % Inputs:
    %   L          : Total length of the domain [0, L]
    %   CN         : Total number of cells
    %   lenRatios  : Vector of length ratios for each sub-region
    %   cellRatios : Vector of cell count ratios for each sub-region
    %   expRatios  : Vector of expansion ratios (end_cell/start_cell) for each sub-region
    %
    % Output:
    %   x          : Vector of nodal coordinates (size: 1 x CN+1)

    % Validate inputs
    numRegions = length(lenRatios);
    if length(cellRatios) ~= numRegions || length(expRatios) ~= numRegions
        error('Length ratios, cell ratios, and expansion ratios must have the same length.');
    end

    % Normalize length and cell ratios
    lenFrac = lenRatios / sum(lenRatios);
    cellFrac = cellRatios / sum(cellRatios);

    % Calculate actual length and cell count per sub-region
    L_sub = L * lenFrac;
    N_sub = round(CN * cellFrac);

    % Correct for potential rounding errors so total cells exactly match CN
    cellDiff = CN - sum(N_sub);
    if cellDiff ~= 0
        % Add or subtract the difference from the region with the most cells
        [~, maxIdx] = max(N_sub);
        N_sub(maxIdx) = N_sub(maxIdx) + cellDiff;
    end

  % Preallocate the node vector as a zeros matrix for performance
    % Total nodes is exactly Total Cells (CN) + 1
    x = zeros(1, CN + 1); 
    
    current_x = 0;
    idx = 1; % Tracks the current index in the preallocated array

    for i = 1:numRegions
        n = N_sub(i);
        l = L_sub(i);
        E = expRatios(i);

        if n == 0
            continue; % Skip if rounding pushed cell count to 0
        end

        if n == 1
            % Only 1 cell in this region
            nodes = current_x + l;
        elseif abs(E - 1.0) < 1e-7
            % Uniform grading (expansion ratio = 1)
            nodes = current_x + linspace(l/n, l, n);
        else
            % Geometric progression grading
            % Note: This mathematics intrinsically handles both E > 1 
            % (cells grow) and E < 1 (cells shrink). 
            % If E < 1, 'r' becomes < 1, making subsequent 'dx' values smaller.
            r = E^(1 / (n - 1));
            
            % Length of the first cell in this sub-region
            delta = l * (r - 1) / (r^n - 1);
            
            % Generate cell lengths
            dx = delta * (r .^ (0 : n-1));
            
            % Calculate node coordinates as a cumulative sum of cell lengths
            nodes = current_x + cumsum(dx);
        end

        % Insert the newly calculated nodes into the preallocated array
        x(idx+1 : idx+n) = nodes;
        
        % Update starting coordinate and index tracker for the next sub-region
        current_x = current_x + l;
        idx = idx + n;
    end
end