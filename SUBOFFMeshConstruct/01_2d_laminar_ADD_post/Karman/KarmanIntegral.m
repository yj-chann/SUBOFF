clc; clear; close all;

nu = 2.87362146e-05;
Uref = 1.649194;     % New reference inflow velocity

% 1. Load Data
load('Mesh.mat');
load('flow_field.mat');
load('DELTA.mat');

% 2. Extract Variables
NS = length(S);      % Number of streamwise stations
N_coord = N;         % Normal coordinate vector

% From flow_field
Us = flow_field.Us;
Un = flow_field.Un; 

% From DELTA
Ue = DELTA.Ue;
delta_star = DELTA.delta_star;
theta = DELTA.theta;
idx99 = DELTA.idx99;

% 3. Calculate Streamwise Derivatives
% Using gradient for numerical differentiation w.r.t the streamwise coordinate S
dUe_ds = gradient(Ue, S);
dr0_ds = gradient(r0, S);
dtheta_ds = gradient(theta, S);

% Preallocate arrays for the terms and the final skin friction coefficient
T1 = zeros(1, NS);
T2 = zeros(1, NS);
T3 = zeros(1, NS);
T4 = zeros(1, NS);
Cf_karman = zeros(1, NS);

% 4. Compute the Karman Integral Equation term by term
for i = 1:NS
    % Prevent division by zero at the stagnation point or axis of symmetry
    if Ue(i) == 0 || r0(i) == 0
        Cf_karman(i) = NaN; % Undefined at exactly zero
        continue;
    end
    
    % --- Term 1 ---
    % (2/Ue) * (dUe/ds) * [ (delta* + 2*theta) + (cos(phi)/r0)*(0.5*delta*^2 + theta^2) ]
    bracket1 = (delta_star(i) + 2*theta(i)) + ...
               (COSPHI(i) / r0(i)) * (0.5 * delta_star(i)^2 + theta(i)^2);
    T1(i) = (2 / Ue(i)) * dUe_ds(i) * bracket1;
    
    % --- Term 2 ---
    % (2/r0) * (dr0/ds) * [ 0.5*kappa0*theta^2 + theta ]
    bracket2 = 0.5 * kappa0(i) * theta(i)^2 + theta(i);
    T2(i) = (2 / r0(i)) * dr0_ds(i) * bracket2;
    
    % --- Term 3 ---
    % 2 * [ 1 + (theta*cos(phi))/r0 ] * (dtheta/ds)
    T3(i) = 2 * (1 + (theta(i) * COSPHI(i)) / r0(i)) * dtheta_ds(i);
    
    % --- Term 4 (Integral Term) ---
    % (2 / (r0*Ue^2)) * \int (kappa0 * r * Us * Un) dn  (Laminar: no Reynolds stress)
    edge = idx99(i);
    n_local = N_coord(1:edge);
    us_local = Us(i, 1:edge);
    un_local = Un(i, 1:edge);
    
    % Local radius r = r0 + n*cos(phi)
    r_local = r0(i) + n_local .* COSPHI(i);
    
    % Integrand: kappa0 * r * (Us * Un)
    integrand = kappa0(i) .* r_local .* (us_local .* un_local);
    
    % Integrate using trapezoidal rule
    integral_T4 = trapz(n_local, integrand);
    T4(i) = (2 / (r0(i) * Ue(i)^2)) * integral_T4;
    
    % --- Local Cf Synthesis ---
    Cf_karman(i) = T1(i) + T2(i) + T3(i) - T4(i);
end

% --- Convert local non-dimensionalization to global Uref non-dimensionalization ---
% The Karman equation naturally yields Cf based on Ue^2. 
% We multiply by (Ue / Uref)^2 to map it to the Uref^2 basis.
scale_factor = (Ue ./ Uref).^2;

Cf_karman = Cf_karman .* scale_factor;
T1 = T1 .* scale_factor;
T2 = T2 .* scale_factor;
T3 = T3 .* scale_factor;
T4 = T4 .* scale_factor;

% --- Direct Wall Gradient Cf ---
% Updated to non-dimensionalize using Uref instead of local Ue
Cf = 2 * (nu * Us(:,1)' / N(1)) / (Uref^2);

% 5. Visualization
figure('Name', 'Karman Integral Cf Verification', 'Color', 'w');
plot(S, Cf_karman, 'k-', 'LineWidth', 1.5, 'DisplayName', 'C_f (Karman Integral)'); hold on
plot(S, Cf, 'b-', 'LineWidth', 1.5, 'DisplayName', 'C_f (Direct Wall Gradient)');
xlabel('Streamwise coordinate S', 'FontSize', 12);
ylabel('C_f (based on U_{ref})', 'FontSize', 12);
title('Calculated Skin Friction Coefficient (Non-dimensionalized by U_{ref})', 'FontSize', 14);
legend('Location', 'best');
grid on;

% Optional: Plot individual terms to see their relative contributions
figure('Name', 'Karman Integral Terms Breakdown', 'Color', 'w');
plot(S, T1, 'b-', 'LineWidth', 1.2, 'DisplayName', 'T_1 (dU_e/ds term)'); hold on;
plot(S, T2, 'g-', 'LineWidth', 1.2, 'DisplayName', 'T_2 (dr_0/ds term)');
plot(S, T3, 'r-', 'LineWidth', 1.2, 'DisplayName', 'T_3 (d\theta/ds term)');
plot(S, -T4, 'm-', 'LineWidth', 1.2, 'DisplayName', '-T_4 (Integral term)');
plot(S, Cf_karman, 'k--', 'LineWidth', 2, 'DisplayName', '\Sigma C_f');
xlabel('Streamwise coordinate S', 'FontSize', 12);
ylabel('Magnitude of Terms (based on U_{ref})', 'FontSize', 12);
title('Breakdown of Karman Integral Terms (Scaled to U_{ref})', 'FontSize', 14);
legend('Location', 'best');
grid on;