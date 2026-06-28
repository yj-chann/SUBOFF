clc; clear; close all;

% --- Global Plot Configurations ---
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultAxesFontName', 'Times New Roman');
set(groot, 'defaultAxesLineWidth', 1.2);

% --- Constants ---
nu = 1.31e-6;
Uref = 1.649194; 
L=4.356*0.2;

% --- Load Data ---
load('ReadData/Mesh.mat');
load('flow_field.mat');
load('DELTA.mat');

Us = flow_field.Us;
Un = flow_field.Un;
p = flow_field.p;
[NS, NN] = size(Us);

% 1. Calculate Skin Friction Coefficient (Cf)
% Since N(1) is the first cell center, distance from wall is N(1).
% Assuming no-slip condition (U = 0 at the wall).
dUdn_wall = Us(:, 1) ./ N(1);

% Calculate kinematic wall shear stress (tau_w / rho)
tau_w_kin = nu .* dUdn_wall;

% Calculate Cf. Density (rho) cancels out because we are using kinematic viscosity.
Cf = tau_w_kin ./ (0.5 * Uref^2);

% 2. Calculate Pressure Coefficient (Cp)
% Extract pressure at the first cell center (closest approximation to wall pressure).
p_wall = p(:, 1);

% Assuming freestream kinematic pressure is 0. Update if your reference is different.
p_inf = 0; 

% Calculate Cp. Density is omitted since p is kinematic.
Cp = (p_wall - p_inf) ./ (0.5 * Uref^2);

% 3. Calculate First Node Dimensionless Wall Distance (y_w^+)
% Calculate friction velocity u_tau = sqrt(tau_w / rho)
u_tau = sqrt(abs(tau_w_kin));

% Calculate y+ at the first cell center. 
yw_plus = (N(1) .* u_tau) ./ nu;

% 4. Plotting the Results
figure('Name', 'Boundary Layer Parameters', 'Position', [100, 100, 800, 800]);

% Plot Cf
subplot(3, 1, 1);
plot(S/L, Cf,'LineWidth',2);
xlabel('$s/L$');
ylabel('$C_f$');
ylim([0,0.012])
grid off;

% Plot Cp
subplot(3, 1, 2);
plot(S/L, Cp,'LineWidth',2);
xlabel('$s/L$');
ylabel('$C_p$');
grid off;

% Plot y_w^+
subplot(3, 1, 3);
plot(S/L, yw_plus,'LineWidth',2);
xlabel('$s/L$');
ylabel('$y_w^+$');
grid off;