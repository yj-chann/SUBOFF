clc;clear;
load('Mesh.mat')
load('flow_field.mat')
NS=981;NN=440;
Us=flow_field.Us;
Ue=zeros(1,NS);
idx99=zeros(1,NS);
for i=1:NS
[Ue(i),idx99(i)]=max(Us(i,:));
end
delta99=N(idx99);


% Preallocate arrays for displacement and momentum thickness
delta_star = zeros(1, NS);
theta = zeros(1, NS);

for i = 1:NS
    % Determine the integration limit index for the current station
    edge_idx = idx99(i);
    
    % Extract the normal coordinate and velocity up to the boundary edge
    n = N(1:edge_idx);
    u = Us(i, 1:edge_idx);
    
    % Get local parameters for this station
    ue_local = Ue(i);
    cosphi_local = COSPHI(i);
    r0_local = r0(i);
    
    % Skip calculation if edge velocity is zero (e.g., exact stagnation point)
    if ue_local == 0
        delta_star(i) = 0;
        theta(i) = 0;
        continue;
    end
    
    % --- Calculate Integrands ---
    % Geometric expansion term: (1 + n*cos(phi)/r0)
    geom_term = 1 + (n .* cosphi_local) ./ r0_local;
    
    % Velocity terms
    vel_ratio = u ./ ue_local;
    vel_deficit = 1 - vel_ratio;
    
    % Integrands for I1 (delta_star) and I2 (theta)
    integrand_I1 = geom_term .* vel_deficit;
    integrand_I2 = geom_term .* vel_ratio .* vel_deficit;
    
    % --- Numerical Integration ---
    % Using the trapezoidal rule across the normal coordinates
    I1 = trapz(n, integrand_I1);
    I2 = trapz(n, integrand_I2);
    
    % --- Calculate Final Thicknesses ---
    % Add a small tolerance check to avoid division by zero if cos(phi) -> 0
    % (When cos(phi) -> 0, it reduces to the planar case where delta* = I1)
    if abs(cosphi_local) < 1e-8
        delta_star(i) = I1;
        theta(i) = I2;
    else
        % Explicit quadratic roots derived previously
        factor = r0_local / cosphi_local;
        root_term_coeff = (2 * cosphi_local) / r0_local;
        
        delta_star(i) = factor * (sqrt(1 + root_term_coeff * I1) - 1);
        theta(i)      = factor * (sqrt(1 + root_term_coeff * I2) - 1);
    end
end

% Optional: Plot the results to verify smoothness
figure;
plot(S, delta99, 'k-', 'LineWidth', 1.5);hold on;
plot(S, delta_star, 'b-', 'LineWidth', 1.5); 
plot(S, theta, 'r-', 'LineWidth', 1.5);
xlabel('Streamwise coordinate S');
ylabel('Thickness');
legend('\delta','\delta^* (Displacement)', '\theta (Momentum)');
title('Boundary Layer Thickness Development');

figure
plot(S,delta_star./theta)

DELTA.delta99=delta99;
DELTA.delta_star=delta_star;
DELTA.theta=theta;
DELTA.idx99=idx99;
DELTA.H=delta_star./theta;
DELTA.Ue=Ue;
save('DELTA.mat',"DELTA")


