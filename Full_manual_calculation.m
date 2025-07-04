% T in K
T_CO_2_in = 273.15 + 50;
T_CO2_out = 273.15 - 40;
T_coolant_in = 273.15 - 50;
% T_coolant_out = 273.15 - 40; % It should be less than -30 Celsius at
% least.

% Flow rates in kg/s
m_dot_CO2 = 0.00167;
m_dot_coolant = 25 * 1000 * 1/60 * 0.87 * 1/1000;

% Material Parameters
k_foam = 0.0035; % W/mK