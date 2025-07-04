% Plate Heat Exchanger PDE Analysis
% CO2 cooling with silicone oil coolant

clear; clc; close all;

%% Input Parameters
% Temperatures
T_CO2_in = 50;      % °C
T_CO2_out = -40;    % °C
T_coolant_in = -60; % °C % tweak this
T_coolant_out = -35; % °C

% Flow rates
m_CO2 = 1.67e-3;    % kg/s (converted from g/s)
V_coolant = 25/60/1000; % m³/s (converted from L/min)

% Plate geometry
plate_area = 0.0138;    % m² % tweak this
plate_ratio = 125/137;  % width/length ratio % tweak this
L_plate = sqrt(plate_area / plate_ratio); % Length
W_plate = L_plate * plate_ratio;          % Width
spacing = 2.4e-3;       % m (plate spacing)
N_plates = 50;         % Number of plates
t_plate = 1e-3;         % Plate thickness (assumed 1mm stainless steel)

% Material properties
% Stainless steel (304)
k_steel = 16.2;         % W/m·K
rho_steel = 8000;       % kg/m³
cp_steel = 500;         % J/kg·K

% CO2 properties (at average temperature 5°C, 1 bar)
rho_CO2 = 1.98;         % kg/m³
cp_CO2 = 844;           % J/kg·K
mu_CO2 = 1.37e-5;       % Pa·s
k_CO2 = 0.0146;         % W/m·K
Pr_CO2 = 0.79;          % Prandtl number

% Silicone oil properties (at average temperature -40°C)
rho_oil = 950;          % kg/m³
cp_oil = 1500;          % J/kg·K
mu_oil = 0.1;           % Pa·s (high viscosity at low temp)
k_oil = 0.15;           % W/m·K
Pr_oil = 1000;          % Prandtl number (high for oils)

fprintf('Plate dimensions: %.3f m x %.3f m\n', L_plate, W_plate);

%% Flow Analysis and Heat Transfer Coefficients

% Hydraulic diameter for rectangular channel
D_h = 2 * spacing;      % For parallel plates

% CO2 flow analysis
A_flow_CO2 = W_plate * spacing; % Flow area per channel
N_channels_CO2 = N_plates / 2;  % Half the plates for CO2
v_CO2 = m_CO2 / (rho_CO2 * A_flow_CO2 * N_channels_CO2);
Re_CO2 = rho_CO2 * v_CO2 * D_h / mu_CO2;

% Coolant flow analysis
m_oil = V_coolant * rho_oil;
A_flow_oil = A_flow_CO2;
N_channels_oil = N_plates / 2;  % Half the plates for coolant
v_oil = m_oil / (rho_oil * A_flow_oil * N_channels_oil);
Re_oil = rho_oil * v_oil * D_h / mu_oil;

fprintf('Reynolds numbers: Re_CO2 = %.1f, Re_oil = %.1f\n', Re_CO2, Re_oil);

% Nusselt number calculations (Dittus-Boelter for turbulent, Sieder-Tate for laminar)
if Re_CO2 > 2300
    Nu_CO2 = 0.023 * Re_CO2^0.8 * Pr_CO2^0.4; % Turbulent
else
    Nu_CO2 = 1.86 * (Re_CO2 * Pr_CO2 * D_h / L_plate)^(1/3); % Laminar VERIFY SIEDER TATE
    Nu_CO2 = max(Nu_CO2, 3.66); % Minimum for fully developed flow
end

if Re_oil > 2300
    Nu_oil = 0.023 * Re_oil^0.8 * Pr_oil^0.4; % Turbulent
else
    Nu_oil = 1.86 * (Re_oil * Pr_oil * D_h / L_plate)^(1/3); % Laminar
    Nu_oil = max(Nu_oil, 3.66);
end

% Heat transfer coefficients
h_CO2 = Nu_CO2 * k_CO2 / D_h;
h_oil = Nu_oil * k_oil / D_h;

fprintf('Heat transfer coefficients: h_CO2 = %.1f W/m²K, h_oil = %.1f W/m²K\n', h_CO2, h_oil);

% Overall heat transfer coefficient
U = 1 / (1/h_CO2 + t_plate/k_steel + 1/h_oil);
fprintf('Overall heat transfer coefficient: U = %.1f W/m²K\n', U);

%% Heat Balance Check
Q_CO2 = m_CO2 * cp_CO2 * (T_CO2_in - T_CO2_out);
Q_oil = m_oil * cp_oil * (T_coolant_out - T_coolant_in);
fprintf('Heat duty: Q_CO2 = %.1f W, Q_oil = %.1f W\n', Q_CO2, Q_oil);

% LMTD calculation
dT1 = T_CO2_in - T_coolant_out;
dT2 = T_CO2_out - T_coolant_in;
LMTD = (dT1 - dT2) / log(dT1 / dT2);
fprintf('LMTD = %.1f K\n', LMTD);

% Required area
A_required = Q_CO2 / (U * LMTD);
A_available = N_plates * plate_area;
fprintf('Required area: %.3f m², Available area: %.3f m²\n', A_required, A_available);
fprintf('Area safety factor: %.2f\n', A_available / A_required);

%% PDE Model Setup
% Create 2D geometry for a single plate
model = createpde('thermal');

% Create rectangular geometry
R1 = [3,4,0,L_plate,L_plate,0,0,0,W_plate,W_plate]';
g = decsg(R1);
geometryFromEdges(model,g);

% Generate mesh
generateMesh(model,'Hmax',min(L_plate,W_plate)/20);

figure(1);
pdemesh(model);
title('Mesh for Single Plate');
axis equal;

%% Material Properties and Boundary Conditions
thermalProperties(model,'ThermalConductivity',k_steel,...
                       'MassDensity',rho_steel,...
                       'SpecificHeat',cp_steel);

% Boundary conditions
% Edge 1: CO2 side (left edge, x=0)
thermalBC(model,'Edge',1,'ConvectionCoefficient',h_CO2,...
          'AmbientTemperature',(T_CO2_in + T_CO2_out)/2);

% Edge 3: Oil side (right edge, x=L_plate)  
thermalBC(model,'Edge',3,'ConvectionCoefficient',h_oil,...
          'AmbientTemperature',(T_coolant_in + T_coolant_out)/2);

% Edges 2 and 4: Adiabatic (top and bottom)
thermalBC(model,'Edge',[2,4],'HeatFlux',0);

% Initial condition
thermalIC(model,0); % Start at 0°C

%% Solve Steady-State PDE
result = solve(model);

% Extract temperature distribution
T_plate = result.Temperature;

%% Plot Results
figure(2);
pdeplot(model,'XYData',T_plate,'ColorMap','jet');
colorbar;
title('Temperature Distribution in Plate (°C)');
xlabel('Length (m)');
ylabel('Width (m)');
axis equal;

figure(3);
pdeplot(model,'XYData',T_plate,'ZData',T_plate,'ColorMap','jet');
title('3D Temperature Distribution');
xlabel('Length (m)');
ylabel('Width (m)');
zlabel('Temperature (°C)');

%% Temperature Profile Analysis
% Extract temperature along centerline
x_center = linspace(0, L_plate, 100);
y_center = W_plate/2 * ones(size(x_center));
T_centerline = interpolateTemperature(result, x_center, y_center);

figure(4);
plot(x_center, T_centerline, 'b-', 'LineWidth', 2);
xlabel('Distance along plate (m)');
ylabel('Temperature (°C)');
title('Temperature Profile along Plate Centerline');
grid on;

%% Uncertainty Analysis using Monte Carlo
fprintf('\n--- Uncertainty Analysis ---\n');

N_mc = 500; % Number of Monte Carlo simulations
uncertainties = struct();

% Define uncertainties (±5% for most parameters)
uncertainties.h_CO2_std = 0.05 * h_CO2;
uncertainties.h_oil_std = 0.05 * h_oil;
uncertainties.k_steel_std = 0.02 * k_steel;
uncertainties.T_CO2_in_std = 2; % ±2°C
uncertainties.T_coolant_in_std = 2; % ±2°C
uncertainties.m_CO2_std = 0.03 * m_CO2; % ±3%

% Monte Carlo simulation
U_mc = zeros(N_mc, 1);
Q_mc = zeros(N_mc, 1);
A_req_mc = zeros(N_mc, 1);

for i = 1:N_mc
    % Sample uncertain parameters
    h_CO2_i = h_CO2 + uncertainties.h_CO2_std * randn();
    h_oil_i = h_oil + uncertainties.h_oil_std * randn();
    k_steel_i = k_steel + uncertainties.k_steel_std * randn();
    T_CO2_in_i = T_CO2_in + uncertainties.T_CO2_in_std * randn();
    T_coolant_in_i = T_coolant_in + uncertainties.T_coolant_in_std * randn();
    m_CO2_i = m_CO2 + uncertainties.m_CO2_std * randn();
    
    % Calculate derived quantities
    U_i = 1 / (1/h_CO2_i + t_plate/k_steel_i + 1/h_oil_i);
    Q_i = m_CO2_i * cp_CO2 * (T_CO2_in_i - T_CO2_out);
    
    % LMTD with perturbed temperatures
    dT1_i = T_CO2_in_i - T_coolant_out;
    dT2_i = T_CO2_out - T_coolant_in_i;
    LMTD_i = (dT1_i - dT2_i) / log(dT1_i / dT2_i);
    
    A_req_i = Q_i / (U_i * LMTD_i);
    
    U_mc(i) = U_i;
    Q_mc(i) = Q_i;
    A_req_mc(i) = A_req_i;
end

% Statistical analysis
fprintf('Monte Carlo Results (N = %d):\n', N_mc);
fprintf('Overall HTC: %.1f ± %.1f W/m²K (95%% CI)\n', mean(U_mc), 1.96*std(U_mc));
fprintf('Heat duty: %.1f ± %.1f W (95%% CI)\n', mean(Q_mc), 1.96*std(Q_mc));
fprintf('Required area: %.3f ± %.3f m² (95%% CI)\n', mean(A_req_mc), 1.96*std(A_req_mc));

% Probability that current design is adequate
P_adequate = sum(A_req_mc <= A_available) / N_mc;
fprintf('Probability of adequate heat transfer: %.1f%%\n', P_adequate * 100);

% Design margin
margin = (A_available - mean(A_req_mc)) / mean(A_req_mc) * 100;
fprintf('Design margin: %.1f%%\n', margin);

%% Uncertainty Distribution Plots
figure(5);
subplot(2,2,1);
histogram(U_mc, 30, 'Normalization', 'probability');
xlabel('Overall HTC (W/m²K)');
ylabel('Probability');
title('Distribution of Overall Heat Transfer Coefficient');

subplot(2,2,2);
histogram(Q_mc, 30, 'Normalization', 'probability');
xlabel('Heat Duty (W)');
ylabel('Probability');
title('Distribution of Heat Duty');

subplot(2,2,3);
histogram(A_req_mc, 30, 'Normalization', 'probability');
xlabel('Required Area (m²)');
ylabel('Probability');
title('Distribution of Required Area');
hold on;
xline(A_available, 'r--', 'LineWidth', 2, 'DisplayName', 'Available Area');
legend();

subplot(2,2,4);
safety_factor = A_available ./ A_req_mc;
histogram(safety_factor, 30, 'Normalization', 'probability');
xlabel('Safety Factor');
ylabel('Probability');
title('Distribution of Safety Factor');
xline(1, 'r--', 'LineWidth', 2, 'DisplayName', 'Critical');
legend();

%% Heat Flux Analysis
% Calculate local heat flux on boundaries
[qx, qy] = evaluateHeatFlux(result);

figure(6);
pdeplot(model,'XYData',sqrt(qx.^2 + qy.^2),'ColorMap','hot');
colorbar;
title('Heat Flux Magnitude (W/m²)');
xlabel('Length (m)');
ylabel('Width (m)');
axis equal;

%% Summary Report
fprintf('\n=== DESIGN SUMMARY ===\n');
fprintf('Plate Heat Exchanger Analysis Results:\n');
fprintf('- Total plates: %d\n', N_plates);
fprintf('- Plate area each: %.4f m²\n', plate_area);
fprintf('- Total heat transfer area: %.3f m²\n', A_available);
fprintf('- Required area (nominal): %.3f m²\n', A_required);
fprintf('- Design is ADEQUATE with %.1f%% probability\n', P_adequate * 100);
fprintf('- Safety margin: %.1f%%\n', margin);
fprintf('- Max plate temperature: %.1f°C\n', max(T_plate));
fprintf('- Min plate temperature: %.1f°C\n', min(T_plate));

if P_adequate >= 0.95
    fprintf('✓ DESIGN RECOMMENDATION: Current design is robust\n');
elseif P_adequate >= 0.8
    fprintf('⚠ DESIGN RECOMMENDATION: Design is marginal, consider increasing area\n');
else
    fprintf('✗ DESIGN RECOMMENDATION: Design inadequate, increase plate count\n');
end