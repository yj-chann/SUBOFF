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
    integrand_D = W .* ((1 / Re_del) * d2Usds2);
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
    
    integrand_D_geo = W .* ( (1 / Re_del) .* (T1 + T2 + T3 - T4 + T5) );
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

figure
plot(S, Cf, 'k', 'LineWidth', 2); hold on;
plot(S, Cf_nu_cart, 'LineWidth', 2);
plot(S, Cf_PG_cart, 'LineWidth', 2);
plot(S, Cf_C_cart, 'LineWidth', 2);
plot(S, Cf_D_cart, 'LineWidth', 2);

figure
plot(S, Cf_nu_cart, 'LineWidth', 2);hold on;
plot(S, Cf_nu_geo, 'LineWidth', 2);

figure
plot(S, Cf_PG_cart, 'LineWidth', 2);hold on;
plot(S, Cf_PG_geo, 'LineWidth', 2);


figure
plot(S, Cf_C_cart, 'LineWidth', 2);hold on;
plot(S, Cf_C_geo, 'LineWidth', 2);
%% 

% Plotting Quick Check (Optional)
figure;

plot(S, Cf, 'k', 'LineWidth', 2); hold on;
plot(S, Cf_cartesian_total, 'b--', 'LineWidth', 1.5);
plot(S, Cf_geometric_total, 'r:', 'LineWidth', 1.5);
plot(S, Cf_cartesian_total+Cf_geometric_total, 'r:', 'LineWidth', 1.5);
xlabel('Streamwise Distance (S)');
ylabel('C_f');
legend('Total C_f', 'Cartesian Base C_f', 'Geometric Correction C_f','fik');
title('FIK Decomposition for Laminar Skin Friction');
grid on;


figure
plot(S,Cf);hold on;
plot(S,Cf_nu+Cf_PG+Cf_C)
% plot(S,Cf_nu)
% plot(S,Cf_PG)
% plot(S,Cf_C)
% plot(S,Cf_D)