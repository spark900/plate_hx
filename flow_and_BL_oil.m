% Input parameters
H = 2.4e-3;               % Plate spacing [m]
visc_coolant = 1.2e-5;    % Kinematic viscosity [m^2/s]
k_coolant = 0.1072;       % Thermal conductivity [W/m·K]
c_p_coolant = 1575;       % Specific heat [J/kg·K]
U = 0.6;                  % Flow velocity [m/s]
rho = 950;                % Assumed density [kg/m^3] (silicone oil)
L = 1.0;                  % Heat exchanger length [m] (adjust as needed)

% Derived properties
Dh = 2 * H;               % Hydraulic diameter [m]
Re = U * Dh / visc_coolant; % Reynolds number
mu = rho * visc_coolant;  % Dynamic viscosity [Pa·s]
Pr = (mu * c_p_coolant) / k_coolant; % Prandtl number

% Entrance lengths
L_hyd = 0.05 * Re * Dh;   % Hydrodynamic entrance length [m]
L_therm = 0.05 * Re * Pr * Dh; % Thermal entrance length [m]

% Calculate thermal boundary layer thickness (delta_t) and h
if L < L_hyd
    % Hydrodynamically developing region
    Re_x = U * L / visc_coolant;
    delta_h = 5.0 * L / sqrt(Re_x); % Hydrodynamic BL thickness
    delta_h = min(delta_h, H/2);     % Cap at half channel height
    delta_t = delta_h / (Pr^(1/3));  % Thermal BL thickness
    delta_t = min(delta_t, H/2);     % Cap at half channel height
    h = k_coolant / delta_t;         % Heat transfer coefficient
elseif L < L_therm
    % Thermally developing region (hydrodynamically developed)
    x_star = L / (Dh * Re * Pr);     % Dimensionless position
    Nu_x = 7.55 + 0.024 * (x_star)^(-0.64) * (Pr)^0.17; % Shah-London correlation
    h = Nu_x * k_coolant / Dh;       % Heat transfer coefficient
    delta_t = k_coolant / h;          % Thermal BL thickness
else
    % Fully developed region
    Nu_fd = 7.54;                    % Nusselt number for parallel plates
    h = Nu_fd * k_coolant / Dh;      % Heat transfer coefficient
    delta_t = k_coolant / h;          % Thermal BL thickness
end

% Nusselt number (based on Dh)
Nu = h * Dh / k_coolant;

% Display results
fprintf('Reynolds number (Re): %.2f (Laminar flow)\n', Re);
fprintf('Prandtl number (Pr): %.2f\n', Pr);
fprintf('Hydrodynamic entrance length: %.4f m\n', L_hyd);
fprintf('Thermal entrance length: %.4f m\n', L_therm);
fprintf('Thermal boundary layer thickness: %.6f m (%.2f mm)\n', delta_t, delta_t*1000);
fprintf('Heat transfer coefficient (h): %.2f W/m²K\n', h);
fprintf('Nusselt number (Nu): %.2f\n', Nu);

% Interpret heat transfer mode
if L < L_therm
    fprintf('\nFlow is thermally DEVELOPING at x = %.2f m.\n', L);
    fprintf('Convection dominates (Nu > 7.54).\n');
    fprintf('Boundary layers are thin; mixing enhances heat transfer.\n');
else
    fprintf('\nFlow is thermally FULLY DEVELOPED at x = %.2f m.\n', L);
    fprintf('Conduction dominates (Nu ≈ 7.54).\n');
    fprintf('Thick boundary layers limit heat transfer; conduction is the bottleneck.\n');
end