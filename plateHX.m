%% Plate Heat Exchanger Analysis with PDE Toolbox
clear; clc; close all;

%% Given Parameters
T_CO2_in = 50;      % °C
T_CO2_out = -40;    % °C
T_coolant_in = -45; % °C
T_coolant_out = -35; % °C
m_CO2 = 1.67e-3;    % kg/s
V_coolant = 25/60/1000; % m³/s (converted from L/min)

% Plate specifications
plate_area = 0.0551;    % m²
plate_ratio = [271, 532]; % width:height ratio
spacing = 2.4e-3;       % m
n_plates = 100;
plate_thickness = 1e-3; % m (assumed stainless steel thickness)

% Calculate plate dimensions
plate_width = sqrt(plate_area * plate_ratio(1) / plate_ratio(2));
plate_height = sqrt(plate_area * plate_ratio(2) / plate_ratio(1));

fprintf('Plate dimensions: %.3f m x %.3f m\n', plate_width, plate_height);
fprintf('Actual plate area: %.4f m²\n', plate_width * plate_height);

%% Fluid Properties (at average temperatures)
% CO2 properties at average temp (5°C)
T_avg_CO2 = (T_CO2_in + T_CO2_out) / 2;
rho_CO2 = 1.98;     % kg/m³
cp_CO2 = 0.85e3;    % J/kg·K
mu_CO2 = 14.8e-6;   % Pa·s
k_CO2 = 0.0166;     % W/m·K
Pr_CO2 = mu_CO2 * cp_CO2 / k_CO2;

% Silicone oil properties at average temp (-40°C)
T_avg_coolant = (T_coolant_in + T_coolant_out) / 2;
rho_coolant = 960;  % kg/m³
cp_coolant = 1.8e3; % J/kg·K
mu_coolant = 0.5;   % Pa·s (high viscosity at low temp)
k_coolant = 0.15;   % W/m·K
Pr_coolant = mu_coolant * cp_coolant / k_coolant;

% Stainless steel properties
k_steel = 16;       % W/m·K
rho_steel = 8000;   % kg/m³
cp_steel = 500;     % J/kg·K

%% Flow Analysis
% Hydraulic diameter
D_h = 2 * spacing;  % m

% Velocity calculations
A_flow_CO2 = plate_width * spacing * (n_plates - 1) / 2; % Half channels for CO2
A_flow_coolant = plate_width * spacing * (n_plates - 1) / 2; % Half channels for coolant

v_CO2 = m_CO2 / (rho_CO2 * A_flow_CO2);
v_coolant = V_coolant / A_flow_coolant;

% Reynolds numbers
Re_CO2 = rho_CO2 * v_CO2 * D_h / mu_CO2;
Re_coolant = rho_coolant * v_coolant * D_h / mu_coolant;

fprintf('\nFlow Analysis:\n');
fprintf('CO2 velocity: %.3f m/s, Re = %.0f\n', v_CO2, Re_CO2);
fprintf('Coolant velocity: %.3f m/s, Re = %.0f\n', v_coolant, Re_coolant);

%% Heat Transfer Coefficients
% Nusselt numbers for parallel plates (Gnielinski correlation)
f_CO2 = (0.79 * log(Re_CO2) - 1.64)^(-2);
Nu_CO2 = (f_CO2/8) * (Re_CO2 - 1000) * Pr_CO2 / ...
         (1 + 12.7 * sqrt(f_CO2/8) * (Pr_CO2^(2/3) - 1));

f_coolant = (0.79 * log(Re_coolant) - 1.64)^(-2);
Nu_coolant = (f_coolant/8) * (Re_coolant - 1000) * Pr_coolant / ...
             (1 + 12.7 * sqrt(f_coolant/8) * (Pr_coolant^(2/3) - 1));

h_CO2 = Nu_CO2 * k_CO2 / D_h;
h_coolant = Nu_coolant * k_coolant / D_h;

fprintf('\nHeat Transfer Coefficients:\n');
fprintf('h_CO2 = %.1f W/m²·K\n', h_CO2);
fprintf('h_coolant = %.1f W/m²·K\n', h_coolant);

%% Overall Heat Transfer Coefficient
R_CO2 = 1/h_CO2;
R_wall = plate_thickness / k_steel;
R_coolant = 1/h_coolant;
R_total = R_CO2 + R_wall + R_coolant;
U = 1/R_total;

fprintf('Overall heat transfer coefficient: %.1f W/m²·K\n', U);

%% Heat Transfer Calculation
Q_CO2 = m_CO2 * cp_CO2 * (T_CO2_in - T_CO2_out);
Q_coolant = rho_coolant * V_coolant * cp_coolant * (T_coolant_out - T_coolant_in);

fprintf('\nHeat Transfer Rates:\n');
fprintf('Q_CO2 = %.1f W\n', Q_CO2);
fprintf('Q_coolant = %.1f W\n', Q_coolant);

Q_avg = (Q_CO2 + Q_coolant) / 2;
fprintf('Average Q = %.1f W\n', Q_avg);

%% Required Area Calculation
LMTD = ((T_CO2_in - T_coolant_out) - (T_CO2_out - T_coolant_in)) / ...
       log((T_CO2_in - T_coolant_out) / (T_CO2_out - T_coolant_in));

A_required = Q_avg / (U * LMTD);
A_available = plate_area * n_plates;

fprintf('\nArea Analysis:\n');
fprintf('LMTD = %.2f K\n', LMTD);
fprintf('Required area: %.3f m²\n', A_required);
fprintf('Available area: %.3f m²\n', A_available);
fprintf('Safety factor: %.2f\n', A_available / A_required);

%% PDE Model Setup
% Create geometry for single plate
rect = [3, 4, 0, plate_width, plate_width, 0, 0, 0, plate_height, plate_height]';
g = decsg(rect);

% Create PDE model
model = createpde('thermal');
geometryFromEdges(model, g);

% Generate mesh
generateMesh(model, 'Hmax', min(plate_width, plate_height)/20);

% Material properties
thermalProperties(model, 'ThermalConductivity', k_steel, ...
                         'MassDensity', rho_steel, ...
                         'SpecificHeat', cp_steel);

%% Boundary Conditions for PDE
% Left edge: CO2 side convection
thermalBC(model, 'Edge', 4, 'ConvectionCoefficient', h_CO2, ...
          'AmbientTemperature', T_avg_CO2);

% Right edge: Coolant side convection  
thermalBC(model, 'Edge', 2, 'ConvectionCoefficient', h_coolant, ...
          'AmbientTemperature', T_avg_coolant);

% Top and bottom edges: insulated (symmetry)
thermalBC(model, 'Edge', [1, 3], 'HeatFlux', 0);

% Solve steady-state PDE
result = solve(model);

%% Plot Temperature Distribution
figure(1);
pdeplot(model, 'XYData', result.Temperature, 'Contour', 'on');
colorbar;
title('Temperature Distribution in Single Plate');
xlabel('Width (m)');
ylabel('Height (m)');

%% Uncertainty Analysis
n_simulations = 1000;
safety_factors = zeros(n_simulations, 1);

% Parameter uncertainties (±10%)
uncertainty = 0.1;

for i = 1:n_simulations
    % Add random variations
    h_CO2_var = h_CO2 * (1 + uncertainty * (2*rand - 1));
    h_coolant_var = h_coolant * (1 + uncertainty * (2*rand - 1));
    Q_var = Q_avg * (1 + uncertainty * (2*rand - 1));
    
    U_var = 1 / (1/h_CO2_var + R_wall + 1/h_coolant_var);
    A_req_var = Q_var / (U_var * LMTD);
    
    safety_factors(i) = A_available / A_req_var;
end

% Calculate probability of being adequate
prob_adequate = sum(safety_factors > 1) / n_simulations;

fprintf('\nUncertainty Analysis:\n');
fprintf('Probability of adequate heat transfer: %.2f%%\n', prob_adequate * 100);
fprintf('Mean safety factor: %.2f ± %.2f\n', mean(safety_factors), std(safety_factors));

%% Plot Histograms
figure(2);
subplot(2,1,1);
histogram(safety_factors, 50);
title('Distribution of Safety Factors');
xlabel('Safety Factor');
ylabel('Frequency');
xline(1, 'r--', 'LineWidth', 2);
legend('Distribution', 'Minimum Required', 'Location', 'best');

subplot(2,1,2);
histogram(1./safety_factors, 50);
title('Distribution of Area Utilization');
xlabel('Area Utilization Ratio');
ylabel('Frequency');
xline(1, 'r--', 'LineWidth', 2);

%% 3D Temperature Distribution Across All Plates
% Create 3D visualization
z_positions = (0:n_plates-1) * (spacing + plate_thickness);
[X, Y] = meshgrid(linspace(0, plate_width, 20), linspace(0, plate_height, 20));

% Temperature varies linearly from CO2 side to coolant side
T_distribution = zeros(size(X, 1), size(X, 2), n_plates);

for k = 1:n_plates
    % Alternate hot and cold sides for counter-flow
    if mod(k, 2) == 1
        T_left = T_avg_CO2;
        T_right = T_avg_coolant;
    else
        T_left = T_avg_coolant;
        T_right = T_avg_CO2;
    end
    
    for i = 1:size(X, 1)
        for j = 1:size(X, 2)
            x_norm = X(i, j) / plate_width;
            T_distribution(i, j, k) = T_left + (T_right - T_left) * x_norm;
        end
    end
end

%% Plot 3D Temperature Distribution
figure(3);
for k = 1:5:n_plates  % Plot every 5th plate for clarity
    Z = z_positions(k) * ones(size(X));
    surf(X, Y, Z, T_distribution(:, :, k), 'EdgeAlpha', 0.3);
    hold on;
end
colorbar;
title('3D Temperature Distribution Across Heat Exchanger');
xlabel('Width (m)');
ylabel('Height (m)');
zlabel('Position (m)');
view(45, 30);

%% Summary Plot
figure(4);
subplot(2,2,1);
plot(result.Temperature);
title('Temperature Profile - Single Plate');
xlabel('Node Number');
ylabel('Temperature (°C)');

subplot(2,2,2);
bar([A_required, A_available]);
title('Area Comparison');
ylabel('Area (m²)');
xticklabels({'Required', 'Available'});

subplot(2,2,3);
plot(z_positions(1:10), squeeze(mean(mean(T_distribution(:,:,1:10), 1), 2)));
title('Average Temperature vs Plate Position');
xlabel('Position (m)');
ylabel('Average Temperature (°C)');

subplot(2,2,4);
plot([h_CO2, h_coolant, U]);
title('Heat Transfer Coefficients');
ylabel('h (W/m²·K)');
xticklabels({'CO2', 'Coolant', 'Overall'});

%% Final Results Summary
fprintf('\n=== FINAL RESULTS ===\n');
if prob_adequate > 0.8
    status = 'ADEQUATE';
else
    status = 'INADEQUATE';
end

fprintf('Heat Exchanger is %s\n', status);
fprintf('Safety factor: %.2f (mean ± std)\n', mean(safety_factors));
fprintf('Confidence level: %.1f%%\n', prob_adequate * 100);
fprintf('Maximum temperature: %.1f°C\n', max(result.Temperature));
fprintf('Minimum temperature: %.1f°C\n', min(result.Temperature));

%% Save results
save('heat_exchanger_results.mat', 'result', 'safety_factors', 'T_distribution', ...
     'A_required', 'A_available', 'prob_adequate', 'U', 'h_CO2', 'h_coolant');