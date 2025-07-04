% Given parameters
m_dot_CO2 = 0.00167;       % CO2 mass flow rate [kg/s]
h_in_CO2 = 424.5e3;        % CO2 inlet enthalpy [J/kg] at 50°C, 10 bar
h_out_CO2 = 46.2e3;        % CO2 outlet enthalpy [J/kg] at -40°C, 10 bar
Q = m_dot_CO2 * (h_in_CO2 - h_out_CO2); % Heat duty [W]

V_dot_coolant = 0.0004167; % Coolant volume flow [m³/s]
rho_coolant = 950;         % Coolant density [kg/m³]
cp_coolant = 1500;         % Coolant specific heat [J/kg.K]
T_coolant_in = -60;        % Coolant inlet temperature [°C]
m_dot_coolant = V_dot_coolant * rho_coolant; % Coolant mass flow [kg/s]
delta_T_coolant = Q / (m_dot_coolant * cp_coolant); % Coolant temp rise [°C]
T_coolant_out = T_coolant_in + delta_T_coolant; % Coolant outlet [°C]

% LMTD calculation (counterflow)
dT1 = 50 - T_coolant_out;  % Temperature difference at hot end
dT2 = -40 - T_coolant_in;  % Temperature difference at cold end
LMTD = (dT1 - dT2) / log(dT1 / dT2); % Log mean temperature difference [°C]

% Overall heat transfer coefficient
U = 1500; % Estimated U [W/m².K]

% Total heat transfer area
A_total = Q / (U * LMTD); % [m²]
A_plate_max = 40e-4;      % Max plate area [m²] (40 cm²)
n_plates = ceil(A_total / A_plate_max) + 1; % Number of plates
A_plate = A_total / (n_plates - 1); % Area per plate [m²]

% Plate dimensions (rectangular)
W = 0.02; % Width [m]
L = A_plate / W; % Length [m]

% Spacing between plates
spacing = 0.003; % [m] (3 mm)

% PDE Model Setup
model = createpde();
rect = [0 L 0 W]; % [x_min, x_max, y_min, y_max]
g = decsg([3; 4; rect(1); rect(2); rect(2); rect(1); rect(3); rect(3); rect(4); rect(4)]);
geometryFromEdges(model, g);

% Material properties (stainless steel)
k_plate = 16; % Thermal conductivity [W/m.K]
t_plate = 0.0005; % Plate thickness [m]

% Heat transfer coefficients
h_CO2 = 3000; % CO2 side [W/m².K]
h_coolant = 3689; % Coolant side [W/m².K]
lambda = (h_CO2 + h_coolant) / (k_plate * t_plate);

% Define f coefficient function
fcoeff = @(location, state) (h_CO2 * (50 - 450 * location.x) + ...
    h_coolant * (-58.936 - 5.32 * location.x)) / (k_plate * t_plate);

% Specify PDE coefficients
specifyCoefficients(model, 'm', 0, 'd', 0, 'c', 1, 'a', -lambda, 'f', fcoeff);

% Boundary conditions (insulated edges)
applyBoundaryCondition(model, 'neumann', 'Edge', 1:4, 'g', 0, 'q', 0);

% Mesh generation
generateMesh(model, 'Hmax', 0.005);

% Solve PDE
result = solvepde(model);
T = result.NodalSolution;

% Plot temperature distribution
figure;
pdeplot(model, 'XYData', T, 'Contour', 'on', 'ColorMap', 'jet');
title('Temperature Distribution on Plate');
xlabel('Length [m]');
ylabel('Width [m]');
colorbar;

% Display results
fprintf('Total heat transfer area: %.4f m²\n', A_total);
fprintf('Number of plates: %d\n', n_plates);
fprintf('Area per plate: %.4f m² (%.1f cm²)\n', A_plate, A_plate * 10000);
fprintf('Plate dimensions: %.3f m x %.3f m\n', L, W);
fprintf('Spacing between plates: %.3f m\n', spacing);