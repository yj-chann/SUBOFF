clc; clear; close all;

nu = 1.31e-6;
Uref = 1.649194;     

load('Mesh.mat');
load('flow_field.mat');
load('DELTA.mat');

Us = flow_field.Us;
Un = flow_field.Un;

% Turbulent Reynolds Stresses
Rss = flow_field.Rss; % \overline{u_s'u_s'}
Rsn = flow_field.Rsn; % \overline{u_s'u_n'}
Rnn = flow_field.Rnn; % \overline{u_n'u_n'}
Rtt = flow_field.Rtt; % \overline{u_\theta'u_\theta'}

p = flow_field.p;
[NS, NN] = size(Us);

idx99 = DELTA.idx99;
delta = DELTA.delta99;

% -------------------------------------------------------------------------
% 1. PRE-COMPUTE GRADIENTS
% -------------------------------------------------------------------------
[dUs_dn, dUs_ds]         = gradient(Us, N, S);
[dUn_dn, dUn_ds]         = gradient(Un, N, S);
[~, d2Us_ds2]            = gradient(dUs_ds, N, S);
[dUsUs_dn, dUsUs_ds]     = gradient(Us.^2, N, S);
[dUsUn_dn, dUsUn_ds]     = gradient(Us.*Un, N, S);
[~, dp_ds]               = gradient(p, N, S); 

% Gradient for Turbulent Rss term needed in C_D
[~, dRss_ds]             = gradient(Rss, N, S);

% -------------------------------------------------------------------------
% 2. PRE-ALLOCATE FIK COMPONENT ARRAYS
% -------------------------------------------------------------------------
% Cartesian (Base) Components
Cf_nu_cart = zeros(1, NS);
Cf_PG_cart = zeros(1, NS);
Cf_C_cart  = zeros(1, NS);
Cf_D_cart  = zeros(1, NS);
Cf_RS_cart = zeros(1, NS); % Added for Turbulence

% Geometric Components
Cf_nu_geo  = zeros(1, NS);
Cf_PG_geo  = zeros(1, NS);
Cf_C_geo   = zeros(1, NS);
Cf_D_geo   = zeros(1, NS);
Cf_RS_geo  = zeros(1, NS); % Added for Turbulence

% -------------------------------------------------------------------------
% 3. MAIN INTEGRATION LOOP OVER STREAMWISE STATIONS
% -------------------------------------------------------------------------
for i = 1:NS
    id_e = idx99(i);
    
    % Local parameters and coordinate extraction
    del  = delta(i);
    n    = N(1:id_e)';            
    u    = Us(i, 1:id_e)';        
    v    = Un(i, 1:id_e)';        
    Ue   = u(end);                
    
    % Local geometric variables
    k0   = kappa0(i);
    dk0  = kappa0_s(i);
    r    = r0(i);
    cosP = COSPHI(i);
    sinP = SINPHI(i);
    
    % Local Reynolds number
    Re_del = Uref * del / nu;
    
    % Local gradients extracted as column vectors
    dpds    = dp_ds(i, 1:id_e)';
    dUs2ds  = dUsUs_ds(i, 1:id_e)';
    dUsUndn = dUsUn_dn(i, 1:id_e)';
    d2Usds2 = d2Us_ds2(i, 1:id_e)';
    dUsds   = dUs_ds(i, 1:id_e)';
    dUnds   = dUn_ds(i, 1:id_e)';
    
    % Local Reynolds stress gradients/variables
    Rsn_local = Rsn(i, 1:id_e)';
    Rss_local = Rss(i, 1:id_e)';
    Rtt_local = Rtt(i, 1:id_e)';
    dRssds    = dRss_ds(i, 1:id_e)';
    
    % Linear weighting function
    W = 1 - n / del;
    
    % =====================================================================
    % CARTESIAN BASE COMPONENTS
    % =====================================================================
    
    % 1. Laminar Viscous Base Term
    Cf_nu_cart(i) = (2 / Re_del) * (Ue / Uref);
    
    % 2. Pressure Gradient Base Term
    integrand_PG = W .* (-dpds / Uref^2);
    Cf_PG_cart(i) = 2 * trapz(n, integrand_PG);
    
    % 3. Convective Base Term
    integrand_C = W .* (-dUs2ds - dUsUndn);
    Cf_C_cart(i) = (2 / Uref^2) * trapz(n, integrand_C);
    
    % 4. Other Derivatives (Diffusion) Base Term (Updated with Rss derivative)
    integrand_D = W .* (nu * d2Usds2 - dRssds);
    Cf_D_cart(i) = (2 / Uref^2) * trapz(n, integrand_D);
    
    % 5. Reynolds Stress Base Term (Turbulent Addition)
    % Note: integral is respect to d(n/del) so we divide trapz by del
    integrand_RS = -Rsn_local / Uref^2;
    Cf_RS_cart(i) = (2 / del) * trapz(n, integrand_RS);

    % =====================================================================
    % GEOMETRIC CURVATURE COMPONENTS
    % =====================================================================
    
    % 1. Laminar Viscous Geometric Component
    term1_nu_geo = (Ue / Uref) * (del * k0 + (1 + del * k0) * del * cosP / r);
    integrand_nu_geo = (k0 + (2 * k0 * n + 1) * cosP / r) .* (u / Uref);
    int_nu_geo = trapz(n, integrand_nu_geo); 
    Cf_nu_geo(i) = (2 / Re_del) * (term1_nu_geo - int_nu_geo);
    
    % 2. Pressure Gradient Geometric Component
    integrand_PG_geo = W .* (-dpds / Uref^2) .* (n * cosP / r);
    Cf_PG_geo(i) = 2 * trapz(n, integrand_PG_geo);
    
    % 3. Convective Geometric Component
    term1_Cgeo = dUs2ds .* (n * cosP / r);
    term2_Cgeo = dUsUndn .* (n * k0 + (1 + n * k0) .* n * cosP / r);
    term3_Cgeo = (1 + n * k0) * sinP / r .* u.^2;
    term4_Cgeo = (2 * k0 + (3 * k0 * n + 1) * cosP / r) .* (u .* v);
    integrand_C_geo = W .* (term1_Cgeo + term2_Cgeo + term3_Cgeo + term4_Cgeo);
    Cf_C_geo(i) = -(2 / Uref^2) * trapz(n, integrand_C_geo);
    
    % 4. Other Derivatives Geometric Component (Updated with turbulent additions)
    T1 = ((n * cosP - r * n * k0) ./ (r * (1 + n * k0))) .* d2Usds2;
    T2 = (sinP / r - (r + n * cosP) .* n ./ (r * (1 + n * k0).^2) * dk0) .* dUsds;
    T3 = (2 * k0 ./ (1 + n * k0)) .* (1 + n * cosP / r) .* dUnds;
    T4 = (k0^2 * (r + n * cosP) ./ (r * (1 + n * k0)) + (1 + n * k0) * sinP^2 ./ (r * (r + n * cosP))) .* u;
    T5 = ((r + n * cosP) ./ (r * (1 + n * k0).^2) * dk0 + sinP * (k0 * r - cosP) ./ (r * (r + n * cosP))) .* v;
    
    % Turbulent specific terms inside C_D^geo
    T6 = -(n * cosP / r) .* dRssds;
    T7 = -(1 + n * k0) .* sinP / r .* (Rss_local - Rtt_local);
    
    integrand_D_geo = W .* (nu .* (T1 + T2 + T3 - T4 + T5) + T6 + T7);
    Cf_D_geo(i) = (2 / Uref^2) * trapz(n, integrand_D_geo);
    
    % 5. Reynolds Stress Geometric Component (Turbulent Addition)
    integrand_RS_geo = (-Rsn_local / Uref^2) .* (k0 * del + (1 + k0 * del) .* n * cosP / r);
    Cf_RS_geo(i) = 2 * trapz(n, integrand_RS_geo) / del;

end

% -------------------------------------------------------------------------
% 4. COMPILE RESULTS
% -------------------------------------------------------------------------
Cf_cartesian_total = Cf_nu_cart + Cf_PG_cart + Cf_C_cart + Cf_D_cart + Cf_RS_cart;
Cf_geometric_total = Cf_nu_geo + Cf_PG_geo + Cf_C_geo + Cf_D_geo + Cf_RS_geo;

Cf_Total = Cf_cartesian_total + Cf_geometric_total;
Cf = 2 * (nu * Us(:,1)' / N(1)) / (Uref^2);

Cf_nu = Cf_nu_cart + Cf_nu_geo;
Cf_PG = Cf_PG_cart + Cf_PG_geo;
Cf_C  = Cf_C_cart + Cf_C_geo;
Cf_D  = Cf_D_cart + Cf_D_geo;
Cf_RS = Cf_RS_cart + Cf_RS_geo; % Combined turbulent contribution


%% 

% =========================================================================
% OPTIMIZED VISUALIZATION SETTINGS
% =========================================================================
% Set default interpreter to LaTeX and font to Times New Roman
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultAxesFontName', 'Times New Roman');
set(groot, 'defaultAxesLineWidth', 1.2);



% -------------------------------------------------------------------------
% Ensure sum validation includes the new turbulent term
% -------------------------------------------------------------------------
Cf_sum = Cf_nu + Cf_PG + Cf_C + Cf_D + Cf_RS;

% -------------------------------------------------------------------------
% Apply Smoothing to FIK Terms (creating new _sm variables)
% -------------------------------------------------------------------------



Cf_nu_sm  = Cf_nu;
Cf_RS_sm  = Cf_RS;
Cf_PG_sm  = Cf_PG;
Cf_C_sm   = Cf_C;
Cf_D_sm   = Cf_D;

win_size = 50;
start_index=450;
Cf_nu_sm(start_index:end)  = smoothdata(Cf_nu(start_index:end), 'movmean', win_size);
Cf_RS_sm(start_index:end)  = smoothdata(Cf_RS(start_index:end), 'movmean', win_size);
Cf_PG_sm(start_index:end)  = smoothdata(Cf_PG(start_index:end), 'movmean', win_size);
Cf_C_sm(start_index:end)   = smoothdata(Cf_C(start_index:end), 'movmean', win_size);
Cf_D_sm(start_index:end)   = smoothdata(Cf_D(start_index:end), 'movmean', win_size);
% -------------------------------------------------------------------------
% FIGURE 1: FIK Identity Verification and Decomposition
% -------------------------------------------------------------------------
L = 4.356 * 0.2;
figure('Position', [400, 100, 900, 400], 'Name', 'FIK Decomposition');

plot(S/L, Cf, 'k', 'LineWidth', 2.5, 'DisplayName', '$C_f$'); hold on;
plot(S(1:25:end-5)/L, Cf_sum(1:25:end-5), 'r', 'LineStyle', 'none', 'LineWidth', 2, ...
    'Marker', 'square', 'MarkerSize', 8, 'DisplayName', '$\sum C_{f,i}$ (FIK)');
plot(S(1:25:end-5)/L, Cf_nu_sm(1:25:end-5), 'b^-', 'LineWidth', 1.5, 'DisplayName', '$C_{\nu}$');
plot(S(1:25:end-5)/L, Cf_RS_sm(1:25:end-5), 'd-', 'Color', [0.8500 0.3250 0.0980], ...
    'LineWidth', 1.5, 'DisplayName', '$C_{RS}$');
plot(S(1:25:end-5)/L, Cf_PG_sm(1:25:end-5), 'g<-', 'LineWidth', 1.5, 'DisplayName', '$C_{PG}$');
plot(S(1:25:end-5)/L, Cf_C_sm(1:25:end-5), 'c>-', 'LineWidth', 1.5, 'DisplayName', '$C_{C}$');
plot(S(1:25:end-5)/L, Cf_D_sm(1:25:end-5), 'mv-', 'LineWidth', 1.5, 'DisplayName', '$C_{D}$');

legend('Location', 'northeast', 'FontSize', 14, 'Interpreter', 'latex');
set(gca, 'FontSize', 14, 'TickLabelInterpreter', 'latex');
xlabel('$s/L$', 'FontSize', 18, 'Interpreter', 'latex');
xlim([0, 0.7]);
ylim([-0.005, 0.015]); 
grid off;


% -------------------------------------------------------------------------
% FIGURE 2: Ratio of the Terms to Total Cf
% -------------------------------------------------------------------------
figure('Position', [400, 100, 900, 400], 'Name', 'FIK Components Ratio');

plot(S(5:25:end-5)/L, Cf_nu_sm(5:25:end-5) ./ Cf(5:25:end-5), 'b^-', 'LineWidth', 1.5, 'DisplayName', '$C_{\nu} / C_f$'); hold on;
plot(S(5:25:end-5)/L, Cf_RS_sm(5:25:end-5) ./ Cf(5:25:end-5), 'd-', 'Color', [0.8500 0.3250 0.0980], ...
    'LineWidth', 1.5, 'DisplayName', '$C_{RS} / C_f$');
plot(S(5:25:end-5)/L, Cf_PG_sm(5:25:end-5) ./ Cf(5:25:end-5), 'g<-', 'LineWidth', 1.5, 'DisplayName', '$C_{PG} / C_f$');
plot(S(5:25:end-5)/L, Cf_C_sm(5:25:end-5) ./ Cf(5:25:end-5), 'c>-', 'LineWidth', 1.5, 'DisplayName', '$C_{C} / C_f$');
plot(S(5:25:end-5)/L, Cf_D_sm(5:25:end-5) ./ Cf(5:25:end-5), 'mv-', 'LineWidth', 1.5, 'DisplayName', '$C_{D} / C_f$');

set(gca, 'FontSize', 14, 'TickLabelInterpreter', 'latex');
xlabel('$s/L$', 'FontSize', 18, 'Interpreter', 'latex');
legend('Location', 'southeast', 'FontSize', 14, 'Interpreter', 'latex');
xlim([0, 0.7]);
ylim([-1, 1.5]);
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
% Note: Cf_cartesian_total and Cf_Total already contain the C_RS and updated C_D terms
Cf_cart_total_sm = smoothdata(Cf_cartesian_total, smooth_method);
Cf_Total_sm      = smoothdata(Cf_Total, smooth_method);

% --- SETUP MAIN FIGURE & AXES ---
figure('Position', [450, 150, 900, 400], 'Name', 'Total Geometric Effect on Cf');
ax = axes('Position', [0.1, 0.15, 0.85, 0.8]);
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
set(ax, 'FontSize', 14, 'LineWidth', 1.2, 'TickLabelInterpreter', 'latex');
xlabel(ax, '$s/L$', 'FontSize', 18, 'Interpreter', 'latex');
legend(ax, 'Location', 'northeast', 'FontSize', 14, 'Interpreter', 'latex');
xlim(ax, [0, 0.7]);
ylim(ax, [0, 0.02]);
grid(ax, 'off');



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

L = 4.356 * 0.2;
smooth_method = 'sgolay'; 
S_row = S(:)';
X_fill = [S_row, fliplr(S_row)] / L;
% =========================================================================
% FIGURE 1: CARTESIAN VS GEOMETRIC BREAKDOWN FOR C_nu 
% =========================================================================
fig1 = figure('Position', [400, 100, 900, 400]);
ax1 = axes('Parent', fig1); hold(ax1, 'on'); box(ax1, 'on');
Cf_nu_cart_sm = smoothdata(Cf_nu_cart, smooth_method);
Cf_nu_sm      = smoothdata(Cf_nu, smooth_method);
Y_fill_nu = [Cf_nu_cart_sm(:)', fliplr(Cf_nu_sm(:)')];

% Geometric fill (using the same logic as C_D and C_RS)
fill(ax1, X_fill, Y_fill_nu, [0.85 0.9 0.95], 'FaceAlpha', 0.8, 'EdgeColor', 'none', 'DisplayName', 'Geometric $C_{\nu}^{geo}$');
plot(ax1, S/L, Cf_nu_sm, '-', 'Color', 'b', 'LineWidth', 2, 'DisplayName', 'Total $C_{\nu}$');
plot(ax1, S/L, Cf_nu_cart_sm, '--', 'Color', [0 0 0.6], 'LineWidth', 2, 'DisplayName', 'Cartesian $C_{\nu}^{cart}$');

set(ax1, 'FontSize', 14, 'LineWidth', 1.2, 'TickLabelInterpreter', 'latex');
xlabel(ax1, '$s/L$', 'FontSize', 18, 'Interpreter', 'latex');
ylabel(ax1, '$C_{\nu}$', 'FontSize', 18, 'Interpreter', 'latex');
legend(ax1, 'Location', 'best', 'FontSize', 12, 'Interpreter', 'latex');
xlim(ax1, [0, 0.7]);

% =========================================================================
% FIGURE 2: CARTESIAN VS GEOMETRIC BREAKDOWN FOR C_PG
% =========================================================================
fig2 = figure('Position', [450, 100, 900, 400]);
ax2 = axes('Parent', fig2); hold(ax2, 'on'); box(ax2, 'on');
Cf_PG_cart_sm = smoothdata(Cf_PG_cart, smooth_method);
Cf_PG_sm      = smoothdata(Cf_PG, smooth_method);
Y_fill_PG = [Cf_PG_cart_sm(:)', fliplr(Cf_PG_sm(:)')];

fill(ax2, X_fill, Y_fill_PG, [0.85 0.95 0.85], 'FaceAlpha', 0.8, 'EdgeColor', 'none', 'DisplayName', 'Geometric $C_{PG}^{geo}$');
plot(ax2, S/L, Cf_PG_sm, '-', 'Color', [0 0.6 0], 'LineWidth', 2, 'DisplayName', 'Total $C_{PG}$');
plot(ax2, S/L, Cf_PG_cart_sm, '--', 'Color', [0 0.3 0], 'LineWidth', 2, 'DisplayName', 'Cartesian $C_{PG}^{cart}$');

yline(ax2, 0, 'k-', 'LineWidth', 0.8, 'HandleVisibility', 'off');
set(ax2, 'FontSize', 14, 'LineWidth', 1.2, 'TickLabelInterpreter', 'latex');
xlabel(ax2, '$s/L$', 'FontSize', 18, 'Interpreter', 'latex');
ylabel(ax2, '$C_{PG}$', 'FontSize', 18, 'Interpreter', 'latex');
legend(ax2, 'Location', 'best', 'FontSize', 12, 'Interpreter', 'latex');
xlim(ax2, [0, 0.7]);

% =========================================================================
% FIGURE 3: CARTESIAN VS GEOMETRIC BREAKDOWN FOR C_C
% =========================================================================
fig3 = figure('Position', [450, 150, 900, 400]);
ax3 = axes('Parent', fig3); hold(ax3, 'on'); box(ax3, 'on');
Cf_C_cart_sm = smoothdata(Cf_C_cart, smooth_method);
Cf_C_sm      = smoothdata(Cf_C, smooth_method);
Y_fill_C = [Cf_C_cart_sm(:)', fliplr(Cf_C_sm(:)')];

fill(ax3, X_fill, Y_fill_C, [0.85 0.95 0.95], 'FaceAlpha', 0.8, 'EdgeColor', 'none', 'DisplayName', 'Geometric $C_{C}^{geo}$');
plot(ax3, S/L, Cf_C_sm, '-', 'Color', [0 0.7 0.7], 'LineWidth', 2, 'DisplayName', 'Total $C_{C}$');
plot(ax3, S/L, Cf_C_cart_sm, '--', 'Color', [0 0.4 0.4], 'LineWidth', 2, 'DisplayName', 'Cartesian $C_{C}^{cart}$');

yline(ax3, 0, 'k-', 'LineWidth', 0.8, 'HandleVisibility', 'off');
set(ax3, 'FontSize', 14, 'LineWidth', 1.2, 'TickLabelInterpreter', 'latex');
xlabel(ax3, '$s/L$', 'FontSize', 18, 'Interpreter', 'latex');
ylabel(ax3, '$C_{C}$', 'FontSize', 18, 'Interpreter', 'latex');
legend(ax3, 'Location', 'best', 'FontSize', 12, 'Interpreter', 'latex');
xlim(ax3, [0, 0.7]);

% =========================================================================
% FIGURE 4: CARTESIAN VS GEOMETRIC BREAKDOWN FOR C_D
% =========================================================================
fig4 = figure('Position', [450, 150, 900, 400]);
ax1 = axes('Parent', fig4);
hold(ax1, 'on'); box(ax1, 'on');

Cf_D_cart_sm = smoothdata(Cf_D_cart, smooth_method);
Cf_D_total_sm = smoothdata(Cf_D, smooth_method);
Y_fill_D = [Cf_D_cart_sm(:)', fliplr(Cf_D_total_sm(:)')];

fill(ax1, X_fill, Y_fill_D, [0.9 0.8 0.95], 'FaceAlpha', 0.8, 'EdgeColor', 'none', ...
    'DisplayName', 'Geometric $C_D^{geo}$');
plot(ax1, S/L, Cf_D_total_sm, '-m', 'LineWidth', 2, 'DisplayName', 'Total $C_D$');
plot(ax1, S/L, Cf_D_cart_sm, '--', 'Color', [0.6 0 0.6], 'LineWidth', 2, ...
    'DisplayName', 'Cartesian $C_D^{cart}$');

yline(ax1, 0, 'k-', 'LineWidth', 0.8, 'HandleVisibility', 'off');
set(ax1, 'FontSize', 14, 'LineWidth', 1.2, 'TickLabelInterpreter', 'latex');
xlabel(ax1, '$s/L$', 'FontSize', 18, 'Interpreter', 'latex');
ylabel(ax1, '$C_D$', 'FontSize', 18, 'Interpreter', 'latex');
legend(ax1, 'Location', 'best', 'FontSize', 12, 'Interpreter', 'latex');
xlim(ax1, [0, 0.7]);

% =========================================================================
% FIGURE 5: CARTESIAN VS GEOMETRIC BREAKDOWN FOR C_RS
% =========================================================================
fig5 = figure('Position', [450, 150, 900, 400], 'Name', 'Reynolds Stress (C_RS) Breakdown');
ax2 = axes('Parent', fig5);
hold(ax2, 'on'); box(ax2, 'on');

Cf_RS_cart_sm = smoothdata(Cf_RS_cart, smooth_method);
Cf_RS_total_sm = smoothdata(Cf_RS, smooth_method);
Y_fill_RS = [Cf_RS_cart_sm(:)', fliplr(Cf_RS_total_sm(:)')];

fill(ax2, X_fill, Y_fill_RS, [0.95 0.85 0.75], 'FaceAlpha', 0.8, 'EdgeColor', 'none', ...
    'DisplayName', 'Geometric $C_{RS}^{geo}$');
plot(ax2, S/L, Cf_RS_total_sm, '-', 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 2, ...
    'DisplayName', 'Total $C_{RS}$');
plot(ax2, S/L, Cf_RS_cart_sm, '--', 'Color', [0.6 0.2 0], 'LineWidth', 2, ...
    'DisplayName', 'Cartesian $C_{RS}^{cart}$');

yline(ax2, 0, 'k-', 'LineWidth', 0.8, 'HandleVisibility', 'off');
set(ax2, 'FontSize', 14, 'LineWidth', 1.2, 'TickLabelInterpreter', 'latex');
xlabel(ax2, '$s/L$', 'FontSize', 18, 'Interpreter', 'latex');
ylabel(ax2, '$C_{RS}$', 'FontSize', 18, 'Interpreter', 'latex');
legend(ax2, 'Location', 'best', 'FontSize', 12, 'Interpreter', 'latex');
xlim(ax2, [0, 0.7]);