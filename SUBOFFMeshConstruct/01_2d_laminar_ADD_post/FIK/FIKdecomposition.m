clc; clear; close all;
nu = 2.87362146e-05;
Uref = 1.649194;     
load('Mesh.mat');
load('flow_field.mat');
load('DELTA.mat');
Us = flow_field.Us;
Un = flow_field.Un;
p = flow_field.p;
[NS, NN] = size(Us);
idx99 = DELTA.idx99;
delta = DELTA.delta99;

% -------------------------------------------------------------------------
% 1. PRE-COMPUTE GRADIENTS
% -------------------------------------------------------------------------
% Note: gradient(F, hN, hS) applies hN to columns (dim 2) and hS to rows (dim 1)
% This maps perfectly to our matrices of size [NS, NN]
[dUs_dn, dUs_ds]         = gradient(Us, N, S);
[dUn_dn, dUn_ds]         = gradient(Un, N, S);
[~, d2Us_ds2]            = gradient(dUs_ds, N, S);
[dUsUs_dn, dUsUs_ds]     = gradient(Us.^2, N, S);
[dUsUn_dn, dUsUn_ds]     = gradient(Us.*Un, N, S);

% Assuming 'p' is kinematic pressure (P/rho). If it's absolute, divide by rho.
[~, dp_ds]               = gradient(p, N, S); 


% -------------------------------------------------------------------------
% 2. PRE-ALLOCATE FIK COMPONENT ARRAYS
% -------------------------------------------------------------------------
% Cartesian (Base) Components
Cf_nu_cart = zeros(1, NS);
Cf_PG_cart = zeros(1, NS);
Cf_C_cart  = zeros(1, NS);
Cf_D_cart  = zeros(1, NS);

% Geometric Components
Cf_nu_geo = zeros(1, NS);
Cf_PG_geo = zeros(1, NS);
Cf_C_geo  = zeros(1, NS);
Cf_D_geo  = zeros(1, NS);

% -------------------------------------------------------------------------
% 3. MAIN INTEGRATION LOOP OVER STREAMWISE STATIONS
% -------------------------------------------------------------------------
for i = 1:NS
    id_e = idx99(i);
    
    % Local parameters and coordinate extraction
    del  = delta(i);
    n    = N(1:id_e)';            % Normal coordinate (column vector)
    u    = Us(i, 1:id_e)';        % Streamwise velocity
    v    = Un(i, 1:id_e)';        % Normal velocity
    Ue   = u(end);                % Edge velocity
    
    % Local geometric variables
    k0   = kappa0(i);
    dk0  = kappa0_s(i);
    r    = r0(i);
    cosP = COSPHI(i);
    sinP = SINPHI(i);
    
    % Local Reynolds number based on momentum thickness/delta
    Re_del = Uref * del / nu;
    
    % Local gradients extracted as column vectors
    dpds    = dp_ds(i, 1:id_e)';
    dUs2ds  = dUsUs_ds(i, 1:id_e)';
    dUsUndn = dUsUn_dn(i, 1:id_e)';
    d2Usds2 = d2Us_ds2(i, 1:id_e)';
    dUsds   = dUs_ds(i, 1:id_e)';
    dUnds   = dUn_ds(i, 1:id_e)';
    
    % Linear weighting function
    W = 1 - n / del;
    
    % =====================================================================
    % CARTESIAN BASE COMPONENTS (Eq. 160)
    % =====================================================================
    
    % 1. Laminar Viscous Base Term
    Cf_nu_cart(i) = (2 / Re_del) * (Ue / Uref);
    
    % 2. Pressure Gradient Base Term
    integrand_PG = W .* (-dpds / Uref^2);
    Cf_PG_cart(i) = 2 * trapz(n, integrand_PG);
    
    % 3. Convective Base Term
    integrand_C = W .* (-dUs2ds - dUsUndn);
    Cf_C_cart(i) = (2 / Uref^2) * trapz(n, integrand_C);
    
    % 4. Other Derivatives (Diffusion) Base Term
    integrand_D = W .* (nu * d2Usds2);
    Cf_D_cart(i) = (2 / Uref^2) * trapz(n, integrand_D);


    % =====================================================================
    % GEOMETRIC CURVATURE COMPONENTS (Eqs. 161 - 165)
    % =====================================================================
    
    % 1. Laminar Viscous Geometric Component (Eq. 161)
    term1_nu_geo = (Ue / Uref) * (del * k0 + (1 + del * k0) * del * cosP / r);
    integrand_nu_geo = (k0 + (2 * k0 * n + 1) * cosP / r) .* (u / Uref);
    % Note: d(n/delta) = dn / delta, so we divide the integral by delta
    int_nu_geo = trapz(n, integrand_nu_geo); 
    Cf_nu_geo(i) = (2 / Re_del) * (term1_nu_geo - int_nu_geo);
    % Cf_nu_geo(i) = (2 / Re_del) * (term1_nu_geo);

    % 2. Pressure Gradient Geometric Component (Eq. 163)
    integrand_PG_geo = W .* (-dpds / Uref^2) .* (n * cosP / r);
    Cf_PG_geo(i) = 2 * trapz(n, integrand_PG_geo);
    
    % 3. Convective Geometric Component (Eq. 164)
    term1_Cgeo = dUs2ds .* (n * cosP / r);
    term2_Cgeo = dUsUndn .* (n * k0 + (1 + n * k0) .* n * cosP / r);
    term3_Cgeo = (1 + n * k0) * sinP / r .* u.^2;
    term4_Cgeo = (2 * k0 + (3 * k0 * n + 1) * cosP / r) .* (u .* v);
    integrand_C_geo = W .* (term1_Cgeo + term2_Cgeo + term3_Cgeo + term4_Cgeo);
    Cf_C_geo(i) = -(2 / Uref^2) * trapz(n, integrand_C_geo);
    
    % 4. Other Derivatives Geometric Component (Eq. 165)
    T1 = ((n * cosP - r * n * k0) ./ (r * (1 + n * k0))) .* d2Usds2;
    T2 = (sinP / r - (r + n * cosP) .* n ./ (r * (1 + n * k0).^2) * dk0) .* dUsds;
    T3 = (2 * k0 ./ (1 + n * k0)) .* (1 + n * cosP / r) .* dUnds;
    T4 = (k0^2 * (r + n * cosP) ./ (r * (1 + n * k0)) + (1 + n * k0) * sinP^2 ./ (r * (r + n * cosP))) .* u;
    T5 = ((r + n * cosP) ./ (r * (1 + n * k0).^2) * dk0 + sinP * (k0 * r - cosP) ./ (r * (r + n * cosP))) .* v;
    
    integrand_D_geo = W .* (nu .* (T1 + T2 + T3 - T4 + T5) );
    Cf_D_geo(i) = (2 / Uref^2) * trapz(n, integrand_D_geo);
    
end

% -------------------------------------------------------------------------
% 4. COMPILE RESULTS
% -------------------------------------------------------------------------
Cf_cartesian_total = Cf_nu_cart + Cf_PG_cart + Cf_C_cart + Cf_D_cart;
Cf_geometric_total = Cf_nu_geo + Cf_PG_geo + Cf_C_geo + Cf_D_geo;

Cf_Total = Cf_cartesian_total + Cf_geometric_total;
Cf = 2 * (nu * Us(:,1)' / N(1)) / (Uref^2);

Cf_nu=Cf_nu_cart+Cf_nu_geo;
Cf_PG=Cf_PG_cart+Cf_PG_geo;
Cf_C=Cf_C_cart+Cf_C_geo;
Cf_D=Cf_D_cart+Cf_D_geo;


% =========================================================================
% OPTIMIZED VISUALIZATION SETTINGS
% =========================================================================
% Set default interpreter to LaTeX and font to Times New Roman
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultAxesFontName', 'Times New Roman');
set(groot, 'defaultAxesLineWidth', 1.2);

% Calculate the sum of the terms to verify FIK identity
Cf_sum = Cf_nu + Cf_PG + Cf_C + Cf_D;


%% 

% -------------------------------------------------------------------------
% FIGURE 1: FIK Identity Verification and Decomposition
% -------------------------------------------------------------------------
L=4.356*0.2;
figure('Position', [400, 100, 900, 400], 'Name', 'FIK Decomposition');
plot(S/L, Cf, 'k', 'LineWidth', 2.5, 'DisplayName', '$C_f$'); hold on;
plot(S(1:25:end-5)/L, Cf_sum(1:25:end-5), 'r','LineStyle','none', 'LineWidth', 2, 'Marker','square',...
    'MarkerSize',8, 'DisplayName', '$\sum C_{f,i}$ (FIK)');
plot(S(1:25:end-5)/L, Cf_nu(1:25:end-5), 'b^-', 'LineWidth', 1.5, 'DisplayName', '$C_{\nu}$');
plot(S(1:25:end-5)/L, Cf_PG(1:25:end-5), 'g<-', 'LineWidth', 1.5, 'DisplayName', '$C_{PG}$');
plot(S(1:25:end-5)/L, Cf_C(1:25:end-5),'c>-', 'LineWidth', 1.5, 'DisplayName', '$C_{C}$');
plot(S(1:25:end-5)/L, Cf_D(1:25:end-5),'mv-', 'LineWidth', 1.5, 'DisplayName', '$C_{D}$');
legend('Location', 'northeast','FontSize',14);
set(gca, 'FontSize', 14);
xlabel('$s/L$','FontSize',18);
xlim([0,0.7])
ylim([-0.01,0.05])
grid off;


% -------------------------------------------------------------------------
% FIGURE 2: Ratio of the Terms to Total Cf
% -------------------------------------------------------------------------
L=4.356*0.2;
figure('Position', [400, 100, 900, 400]);
plot(S(5:25:end-5)/L, Cf_nu(5:25:end-5) ./ Cf(5:25:end-5),'b^-', 'LineWidth', 1.5, 'DisplayName', '$C_{\nu} / C_f$'); hold on;
plot(S(5:25:end-5)/L, Cf_PG(5:25:end-5) ./ Cf(5:25:end-5),'g<-', 'LineWidth', 1.5, 'DisplayName', '$C_{PG} / C_f$');
plot(S(5:25:end-5)/L, Cf_C(5:25:end-5) ./ Cf(5:25:end-5),'c>-', 'LineWidth', 1.5, 'DisplayName', '$C_{C} / C_f$');
plot(S(5:25:end-5)/L, Cf_D(5:25:end-5) ./ Cf(5:25:end-5),'mv-', 'LineWidth', 1.5, 'DisplayName', '$C_{D} / C_f$');

set(gca, 'FontSize', 14);
xlabel('$s/L$','FontSize',18);
legend('Location', 'southeast','FontSize',14);
xlim([0,0.7])
ylim([-1.5,1.5])
grid off;
%% 


% =========================================================================
% FIGURE 3: TOTAL CARTESIAN VS TOTAL GEOMETRIC EFFECT ON Cf
% =========================================================================
L = 4.356 * 0.2;
color_total = 'k';              
color_base  = [0.8 0 0];        
color_shade = [0.95 0.8 0.8];   
smooth_method = 'sgolay'; 

S_row = S(:)';
X_fill = [S_row, fliplr(S_row)] / L;

% --- SMOOTHING STEP ---
% Using the already compiled Cf_cartesian_total and Cf_Total arrays
Cf_cart_total_sm = smoothdata(Cf_cartesian_total, smooth_method);
Cf_Total_sm      = smoothdata(Cf_Total, smooth_method);

% --- SETUP MAIN FIGURE & AXES ---
fig5 = figure('Position', [450, 150, 900, 400], 'Name', 'Total Geometric Effect on Cf');
ax= axes('Position', [0.1, 0.15, 0.85, 0.8]);
hold(ax, 'on'); box(ax, 'on');

Y_fill_Total = [Cf_cart_total_sm(:)', fliplr(Cf_Total_sm(:)')];

% --- PLOT MAIN DATA ---
plot(ax, S/L, Cf_Total_sm, '-', 'Color', color_total, 'LineWidth', 2, ...
    'DisplayName', 'Total $C_f$');
plot(ax, S/L, Cf_cart_total_sm, '--', 'Color', color_base, 'LineWidth', 2, ...
    'DisplayName', 'Cartesian Term $C_f^{cart}$');
fill(ax, X_fill, Y_fill_Total, color_shade, 'FaceAlpha', 0.8, 'EdgeColor', 'none', ...
    'DisplayName', 'Geometric Effect $C_f^{geo}$');

yline(ax, 0, 'k-', 'LineWidth', 0.8, 'HandleVisibility', 'off');

% Format Main Axes
set(ax, 'FontSize', 14, 'LineWidth', 1.2);
xlabel(ax, '$s/L$', 'FontSize', 18, 'Interpreter', 'latex');
legend(ax, 'Location', 'northeast', 'FontSize', 14, 'Interpreter', 'latex');
xlim(ax, [0, 0.7]);
ylim(ax, [0, 0.07]);
grid(ax, 'off');

%% 

% =========================================================================
% FIGURE: PERCENTAGE OF CARTESIAN AND GEOMETRIC EFFECTS ON Cf
% =========================================================================
L = 4.356 * 0.2;
% --- CALCULATE PERCENTAGES ---
% Calculate the smoothed geometric total
Cf_geo_total_sm = Cf_Total_sm - Cf_cart_total_sm;

% Calculate percentage contributions
Pct_cart = (Cf_cart_total_sm ./ Cf_Total_sm);
Pct_geo  = (Cf_geo_total_sm ./ Cf_Total_sm);

% --- SETUP MAIN FIGURE & AXES ---
fig_pct = figure('Position', [450, 150, 900, 400], 'Name', 'Percentage Breakdown of Cf');
ax_pct = axes('Position', [0.1, 0.15, 0.85, 0.8]);
hold(ax_pct, 'on'); box(ax_pct, 'on');

% --- PLOT PERCENTAGE DATA ---
plot(ax_pct, S(15:25:end-5)/L, Pct_cart(15:25:end-5), 's-', 'Color', 'r', 'LineWidth', 2, ...
    'MarkerSize',8,'DisplayName', '$C_f^{cart} / C_f$');

% Using a distinct blue color for the geometric line to differentiate from the base red
plot(ax_pct, S(15:25:end-5)/L, Pct_geo(15:25:end-5), 's-', 'Color', 'b', 'LineWidth', 2, ...
    'MarkerSize',8,'DisplayName', '$C_f^{geo} / C_f$');

% Reference lines at 0% and 100%
yline(ax_pct, 0, 'k-', 'LineWidth', 0.8, 'HandleVisibility', 'off');
yline(ax_pct, 100, 'k--', 'LineWidth', 0.8, 'HandleVisibility', 'off', 'Color', [0.5 0.5 0.5]); 

% Format Main Axes
set(ax_pct, 'FontSize', 14, 'LineWidth', 1.2);
xlabel(ax_pct, '$s/L$', 'FontSize', 18, 'Interpreter', 'latex');
legend(ax_pct, 'Location', 'best', 'FontSize', 14, 'Interpreter', 'latex');
xlim(ax_pct, [0, 0.7]);
ylim(ax_pct, [-1, 2]); 
grid(ax_pct, 'off');

%% 
% =========================================================================
% PUBLICATION-QUALITY SETTINGS
% =========================================================================
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultAxesFontName', 'Times New Roman');
set(groot, 'defaultAxesLineWidth', 1.2);
set(groot, 'defaultAxesFontSize', 14);



% =========================================================================
% FIGURE 1: VISCOUS TERM WITH DUAL ZOOM INSETS
% =========================================================================
L = 4.356 * 0.2;

% Match colors exactly
color_total = 'k';              % Black for Total
color_base  = [0.8 0 0];        % Red for Cartesian Base
color_shade = [0.95 0.8 0.8];   % Pinkish for Geometric Shade

% --- SMOOTHING STEP ---
smooth_method = 'sgolay'; 
Cf_nu_cart_sm = smoothdata(Cf_nu_cart, smooth_method);
Cf_nu_sm      = smoothdata(Cf_nu, smooth_method);

% --- 1. SETUP MAIN FIGURE & AXES ---
fig1 = figure('Position', [400, 100, 900, 400]);
ax_main = axes('Position', [0.1, 0.15, 0.85, 0.8]); % Define explicit position
hold(ax_main, 'on'); box(ax_main, 'on');

% Prepare the X-coordinates for the fill polygon
S_row = S(:)';
X_fill = [S_row, fliplr(S_row)] / L;
Y_fill_nu = [Cf_nu_cart_sm(:)', fliplr(Cf_nu_sm(:)')];

% --- 2. PLOT MAIN DATA ---
plot(ax_main, S/L, Cf_nu_sm, '-', 'Color', color_total, 'LineWidth', 2, ...
    'DisplayName', 'Total $C_{\nu}$');
plot(ax_main, S/L, Cf_nu_cart_sm, '--', 'Color', color_base, 'LineWidth', 2, ...
    'DisplayName', 'Cartesian Term $C_{\nu}^{cart}$');
fill(ax_main, X_fill, Y_fill_nu, color_shade, 'FaceAlpha', 0.8, 'EdgeColor', 'none', ...
    'DisplayName', 'Geometric Effect $C_{\nu}^{geo}$');



% Format Main Axes
set(ax_main, 'FontSize', 14, 'LineWidth', 1.2);
xlabel(ax_main, '$s/L$', 'FontSize', 18);
legend(ax_main, 'Location', 'northeast', 'FontSize', 14);
xlim(ax_main, [0, 0.7]);
ylim(ax_main, [0.0, 0.02]);
grid(ax_main, 'off');

% --- Helper Functions for Coordinate Mapping ---
% Converts Data Coordinates to Normalized Figure Coordinates
xl = xlim(ax_main);
yl = ylim(ax_main);
pos = get(ax_main, 'Position');
x2norm = @(x) pos(1) + (x - xl(1)) / (xl(2) - xl(1)) * pos(3);
y2norm = @(y) pos(2) + (y - yl(1)) / (yl(2) - yl(1)) * pos(4);

% =========================================================================
% --- 3. FIRST ZOOM REGION (Peak area) ---
% =========================================================================
x_zoom1 = [0.015, 0.055]; 
y_zoom1 = [0.016, 0.019];

% Draw bounding box
plot(ax_main, [x_zoom1(1) x_zoom1(2) x_zoom1(2) x_zoom1(1) x_zoom1(1)], ...
              [y_zoom1(1) y_zoom1(1) y_zoom1(2) y_zoom1(2) y_zoom1(1)], ...
              'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');

% Create Inset Axes 1
inset_pos1 = [0.25, 0.55, 0.12, 0.35]; 
ax_inset1 = axes('Position', inset_pos1);
hold(ax_inset1, 'on'); box(ax_inset1, 'on');

fill(ax_inset1, X_fill, Y_fill_nu, color_shade, 'FaceAlpha', 0.8, 'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(ax_inset1, S/L, Cf_nu_cart_sm, '--', 'Color', color_base, 'LineWidth', 2, 'HandleVisibility', 'off');
plot(ax_inset1, S/L, Cf_nu_sm, '-', 'Color', color_total, 'LineWidth', 2, 'HandleVisibility', 'off');

xlim(ax_inset1, x_zoom1); ylim(ax_inset1, y_zoom1);
set(ax_inset1, 'XTick', [], 'YTick', [], 'LineWidth', 1.2);

% Connecting lines for Inset 1
annotation('line', [x2norm(x_zoom1(2)), inset_pos1(1)], [y2norm(y_zoom1(2)), inset_pos1(2) + inset_pos1(4)], 'LineStyle', '--', 'LineWidth', 1.5, 'Color', 'k');
annotation('line', [x2norm(x_zoom1(2)), inset_pos1(1)], [y2norm(y_zoom1(1)), inset_pos1(2)], 'LineStyle', '--', 'LineWidth', 1.5, 'Color', 'k');

% =========================================================================
% --- 4. SECOND ZOOM REGION (x = 0.48 to 0.52) ---
% =========================================================================
x_zoom2 = [0.48, 0.52];
% NOTE: Adjust y_zoom2 based on the actual min/max of your data in this specific range!
y_zoom2 = [0.0015, 0.0025]; 

% Draw bounding box
plot(ax_main, [x_zoom2(1) x_zoom2(2) x_zoom2(2) x_zoom2(1) x_zoom2(1)], ...
              [y_zoom2(1) y_zoom2(1) y_zoom2(2) y_zoom2(2) y_zoom2(1)], ...
              'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');

% Create Inset Axes 2 (Placed in the middle-top area)
inset_pos2 = [0.60, 0.40, 0.12, 0.2]; 
ax_inset2 = axes('Position', inset_pos2);
hold(ax_inset2, 'on'); box(ax_inset2, 'on');

fill(ax_inset2, X_fill, Y_fill_nu, color_shade, 'FaceAlpha', 0.8, 'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(ax_inset2, S/L, Cf_nu_cart_sm, '--', 'Color', color_base, 'LineWidth', 2, 'HandleVisibility', 'off');
plot(ax_inset2, S/L, Cf_nu_sm, '-', 'Color', color_total, 'LineWidth', 2, 'HandleVisibility', 'off');

xlim(ax_inset2, x_zoom2); ylim(ax_inset2, y_zoom2);
set(ax_inset2, 'XTick', [], 'YTick', [], 'LineWidth', 1.2);

% Connecting lines for Inset 2 (From UPPER EDGE of the dashed box to BOTTOM CORNERS of inset)
% Line 1: Top-Left of bounding box -> Bottom-Left of inset
annotation('line', [x2norm(x_zoom2(1)), inset_pos2(1)], ...
                   [y2norm(y_zoom2(2)), inset_pos2(2)], ...
                   'LineStyle', '--', 'LineWidth', 1.5, 'Color', 'k');

% Line 2: Top-Right of bounding box -> Bottom-Right of inset
annotation('line', [x2norm(x_zoom2(2)), inset_pos2(1) + inset_pos2(3)], ...
                   [y2norm(y_zoom2(2)), inset_pos2(2)], ...
                   'LineStyle', '--', 'LineWidth', 1.5, 'Color', 'k');


% =========================================================================
% FIGURE 2: PRESSURE GRADIENT TERM
% =========================================================================
L = 4.356 * 0.2;
color_total = 'k';              
color_base  = [0.8 0 0];        
color_shade = [0.95 0.8 0.8];   
smooth_method = 'sgolay'; 

S_row = S(:)';
X_fill = [S_row, fliplr(S_row)] / L;
% --- SMOOTHING STEP ---
Cf_PG_cart_sm = smoothdata(Cf_PG_cart, smooth_method);
Cf_PG_sm      = smoothdata(Cf_PG, smooth_method);

% --- 1. SETUP MAIN FIGURE & AXES ---
fig2 = figure('Position', [450, 150, 900, 400], 'Name', 'PG Term');
ax_main2 = axes('Position', [0.1, 0.15, 0.85, 0.8]);
hold(ax_main2, 'on'); box(ax_main2, 'on');

Y_fill_PG = [Cf_PG_cart_sm(:)', fliplr(Cf_PG_sm(:)')];

% --- 2. PLOT MAIN DATA ---
plot(ax_main2, S/L, Cf_PG_sm, '-', 'Color', color_total, 'LineWidth', 2, ...
    'DisplayName', 'Total $C_{PG}$');
plot(ax_main2, S/L, Cf_PG_cart_sm, '--', 'Color', color_base, 'LineWidth', 2, ...
    'DisplayName', 'Cartesian Term $C_{PG}^{cart}$');
fill(ax_main2, X_fill, Y_fill_PG, color_shade, 'FaceAlpha', 0.8, 'EdgeColor', 'none', ...
    'DisplayName', 'Geometric Effect $C_{PG}^{geo}$');


yline(ax_main2, 0, 'k-', 'LineWidth', 0.8, 'HandleVisibility', 'off');

% Format Main Axes
set(ax_main2, 'FontSize', 14, 'LineWidth', 1.2);
xlabel(ax_main2, '$s/L$', 'FontSize', 18);
legend(ax_main2, 'Location', 'northeast', 'FontSize', 14);
xlim(ax_main2, [0, 0.7]);
% ylim(ax_main2, [-0.01, 0.04]); % TWEAK THIS LIMIT BASED ON YOUR PG DATA
grid(ax_main2, 'off');

% --- 3. ZOOM REGION ---
xl2 = xlim(ax_main2); yl2 = ylim(ax_main2); pos2 = get(ax_main2, 'Position');
x2norm_2 = @(x) pos2(1) + (x - xl2(1)) / (xl2(2) - xl2(1)) * pos2(3);
y2norm_2 = @(y) pos2(2) + (y - yl2(1)) / (yl2(2) - yl2(1)) * pos2(4);

x_zoom_PG = [0.010, 0.038];  % TWEAK LOCATION
y_zoom_PG = [0.028, 0.035];  % TWEAK LOCATION

plot(ax_main2, [x_zoom_PG(1) x_zoom_PG(2) x_zoom_PG(2) x_zoom_PG(1) x_zoom_PG(1)], ...
               [y_zoom_PG(1) y_zoom_PG(1) y_zoom_PG(2) y_zoom_PG(2) y_zoom_PG(1)], ...
               'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');

inset_pos_PG = [0.25, 0.6, 0.1, 0.3]; 
ax_inset2 = axes('Position', inset_pos_PG);
hold(ax_inset2, 'on'); box(ax_inset2, 'on');

fill(ax_inset2, X_fill, Y_fill_PG, color_shade, 'FaceAlpha', 0.8, 'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(ax_inset2, S/L, Cf_PG_cart_sm, '--', 'Color', color_base, 'LineWidth', 2, 'HandleVisibility', 'off');
plot(ax_inset2, S/L, Cf_PG_sm, '-', 'Color', color_total, 'LineWidth', 2, 'HandleVisibility', 'off');

xlim(ax_inset2, x_zoom_PG); ylim(ax_inset2, y_zoom_PG);
set(ax_inset2, 'XTick', [], 'YTick', [], 'LineWidth', 1.2);

annotation('line', [x2norm_2(x_zoom_PG(2)), inset_pos_PG(1)], [y2norm_2(y_zoom_PG(2)), inset_pos_PG(2) + inset_pos_PG(4)], 'LineStyle', '--', 'LineWidth', 1.5, 'Color', 'k');
annotation('line', [x2norm_2(x_zoom_PG(2)), inset_pos_PG(1)], [y2norm_2(y_zoom_PG(1)), inset_pos_PG(2)], 'LineStyle', '--', 'LineWidth', 1.5, 'Color', 'k');

% =========================================================================
% --- 4. SECOND ZOOM REGION  ---
% =========================================================================
x_zoom2 = [0.22, 0.3];
y_zoom2 = [-0.006, -0.003]; 

% Draw bounding box
plot(ax_main2, [x_zoom2(1) x_zoom2(2) x_zoom2(2) x_zoom2(1) x_zoom2(1)], ...
              [y_zoom2(1) y_zoom2(1) y_zoom2(2) y_zoom2(2) y_zoom2(1)], ...
              'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');

% Create Inset Axes 3 (Placed in the middle-top area)
% (Renamed to ax_inset3 to avoid overwriting the first inset ax_inset2)
inset_pos2 = [0.50, 0.4, 0.18, 0.1]; 
ax_inset3 = axes('Position', inset_pos2);
hold(ax_inset3, 'on'); box(ax_inset3, 'on');

% FIXED: Changed _nu variables to _PG variables
fill(ax_inset3, X_fill, Y_fill_PG, color_shade, 'FaceAlpha', 0.8, 'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(ax_inset3, S/L, Cf_PG_sm, '-', 'Color', color_total, 'LineWidth', 2, 'HandleVisibility', 'off');
plot(ax_inset3, S/L, Cf_PG_cart_sm, '--', 'Color', color_base, 'LineWidth', 2, 'HandleVisibility', 'off');

xlim(ax_inset3, x_zoom2); ylim(ax_inset3, y_zoom2);
set(ax_inset3, 'XTick', [], 'YTick', [], 'LineWidth', 1.2);

% FIXED: Changed x2norm -> x2norm_2 and y2norm -> y2norm_2
% Connecting lines for Inset 2 (From UPPER EDGE of the dashed box to BOTTOM CORNERS of inset)
% Line 1: Top-Left of bounding box -> Bottom-Left of inset
annotation('line', [x2norm_2(x_zoom2(1)), inset_pos2(1)], ...
                   [y2norm_2(y_zoom2(2)), inset_pos2(2)], ...
                   'LineStyle', '--', 'LineWidth', 1.5, 'Color', 'k');

% Line 2: Top-Right of bounding box -> Bottom-Right of inset
annotation('line', [x2norm_2(x_zoom2(2)), inset_pos2(1) + inset_pos2(3)], ...
                   [y2norm_2(y_zoom2(2)), inset_pos2(2)], ...
                   'LineStyle', '--', 'LineWidth', 1.5, 'Color', 'k');



% =========================================================================
% FIGURE 3: CONVECTIVE TERM
% =========================================================================

% --- SMOOTHING STEP ---
L = 4.356 * 0.2;
color_total = 'k';              
color_base  = [0.8 0 0];        
color_shade = [0.95 0.8 0.8];   
smooth_method = 'sgolay'; 

Cf_C_cart_sm = smoothdata(Cf_C_cart, smooth_method);
Cf_C_sm      = smoothdata(Cf_C, smooth_method);

% --- 1. SETUP MAIN FIGURE & AXES ---
fig3 = figure('Position', [450, 150, 900, 400], 'Name', 'Convective Geometric Effect');
ax_main3 = axes('Position', [0.1, 0.15, 0.85, 0.8]);
hold(ax_main3, 'on'); box(ax_main3, 'on');

Y_fill_C = [Cf_C_cart_sm(:)', fliplr(Cf_C_sm(:)')];

% --- 2. PLOT MAIN DATA ---
plot(ax_main3, S/L, Cf_C_sm, '-', 'Color', color_total, 'LineWidth', 2, ...
    'DisplayName', 'Total $C_{C}$');
plot(ax_main3, S/L, Cf_C_cart_sm, '--', 'Color', color_base, 'LineWidth', 2, ...
    'DisplayName', 'Cartesian Term $C_{C}^{cart}$');
fill(ax_main3, X_fill, Y_fill_C, color_shade, 'FaceAlpha', 0.8, 'EdgeColor', 'none', ...
    'DisplayName', 'Geometric Effect $C_{C}^{geo}$');

yline(ax_main3, 0, 'k-', 'LineWidth', 0.8, 'HandleVisibility', 'off');

% Format Main Axes
set(ax_main3, 'FontSize', 14, 'LineWidth', 1.2);
xlabel(ax_main3, '$s/L$', 'FontSize', 18);
legend(ax_main3, 'Location', 'northeast', 'FontSize', 14);
xlim(ax_main3, [0, 0.7]);
grid(ax_main3, 'off');




% =========================================================================
% FIGURE 4: DIFFUSION TERM
% =========================================================================
L = 4.356 * 0.2;
color_total = 'k';              
color_base  = [0.8 0 0];        
color_shade = [0.95 0.8 0.8];   
smooth_method = 'sgolay'; 

S_row = S(:)';
X_fill = [S_row, fliplr(S_row)] / L;

% --- SMOOTHING STEP ---
Cf_D_cart_sm = smoothdata(Cf_D_cart, smooth_method);
Cf_D_sm      = smoothdata(Cf_D, smooth_method);

% --- 1. SETUP MAIN FIGURE & AXES ---
fig4 = figure('Position', [450, 150, 900, 400], 'Name', 'Diffusion Geometric Effect');
ax_main4 = axes('Position', [0.1, 0.15, 0.85, 0.8]);
hold(ax_main4, 'on'); box(ax_main4, 'on');

Y_fill_D = [Cf_D_cart_sm(:)', fliplr(Cf_D_sm(:)')];

% --- 2. PLOT MAIN DATA ---
plot(ax_main4, S/L, Cf_D_sm, '-', 'Color', color_total, 'LineWidth', 2, ...
    'DisplayName', 'Total $C_{D}$');
plot(ax_main4, S/L, Cf_D_cart_sm, '--', 'Color', color_base, 'LineWidth', 2, ...
    'DisplayName', 'Cartesian Term $C_{D}^{cart}$');
fill(ax_main4, X_fill, Y_fill_D, color_shade, 'FaceAlpha', 0.8, 'EdgeColor', 'none', ...
    'DisplayName', 'Geometric Effect $C_{D}^{geo}$');

yline(ax_main4, 0, 'k-', 'LineWidth', 0.8, 'HandleVisibility', 'off');

% Format Main Axes
set(ax_main4, 'FontSize', 14, 'LineWidth', 1.2);
xlabel(ax_main4, '$s/L$', 'FontSize', 18, 'Interpreter', 'latex');
legend(ax_main4, 'Location', 'northeast', 'FontSize', 14, 'Interpreter', 'latex');
xlim(ax_main4, [0, 0.7]);
ylim(ax_main4, [-0.0001, 0.0001]);
grid(ax_main4, 'off');



