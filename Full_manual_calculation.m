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
k_steel = 15; % W/mK
visc_coolant = 1.2e-5; %m^2/s
k_coolant = 0.1072; % W/mk
c_p_coolant = 1.575e3; % J/kgK


% Logarithmic Temperature difference calculation: VERIFY !
deltaT1 = T_CO_2_in - T_coolant_in;
deltaT2 = T_CO2_out - T_coolant_in;
LMTD = (deltaT1 - deltaT2) / log(deltaT1 / deltaT2);

% Dittus-Boelter equation for Nusselt number
Re = (m_dot_coolant * k_coolant) / (visc_coolant * (T_CO_2_in - T_CO2_out)); % Reynolds number
Pr = (c_p_coolant * visc_coolant) / k_coolant; % Prandtl number
Nu = 0.023 * Re^0.8 * Pr^0.3; % Nusselt number for turbulent flow

% Calculate heat transfer coefficient
h = Nu * k_coolant / (T_CO_2_in - T_CO2_out); % Heat transfer coefficient

% Conduction through coolant
A = 1; % Assume a unit area for simplicity
Q_conduction = k_coolant * A * (T_CO_2_in - T_CO2_out); % Heat transfer by conduction