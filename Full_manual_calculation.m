% T in K
T_CO2_in = 273.15 + 50;
T_CO2_out = 273.15 - 40;
T_coolant_in = 273.15 - 50;
T_coolant_out = 273.15 - 30; % It should be less than -30 Celsius at
% least.

% Flow rates in kg/s
m_dot_CO2 = 0.00167;
m_dot_coolant = 25 * 1000 * 1/60 * 0.87 * 1/1000;
density_coolant = 0.87;
% Assuming the ZC6 plate heat exchanger and computing the heat transfer for
% one plate. The goal then is to find how many plates we need after knowing
% all the data for one plate. 
length = 0.271; % Breite eines Kanals zwischen zwei Platte (A) 
width = 0.0024; % Abstand zwischen Zwei Platten.
height = 0.532; % Höhe der Platten (B)
u_coolant = m_dot_coolant / (density_coolant * length * width);


% Material Parameters
k_foam = 0.0035; % W/mK
k_steel = 15; % W/mK
visc_coolant = 1.2e-5; %m^2/s
k_coolant = 0.1072; % W/mk
c_p_coolant = 1.575e3; % J/kgK


% T_coolant calculation
%{ 
The issue with this being that the mixing is very low for mineral oil. The
outer layer of mineral oil is barely moving and convection may be not
dominant with a warm layer of oil at the plates and a cold core in the
center.
Distance between two plates: 2.4 mm = 0.0024 m 
%}


T_coolant = T_coolant_in; % T_coolant_in must be found depending on 
% location in the heat exchanger.
% Logarithmic Temperature difference calculation
deltaT1 = T_coolant - T_CO2_in;
deltaT2 = T_coolant - T_CO2_out;
LMTD = (deltaT1 - deltaT2) / log(deltaT1 / deltaT2);

characteristic_len = (4 * (length * width)) / (2 * height + 2 * length);

Re = (u_coolant * characteristic_len) / visc_coolant; % Reynolds number
Pr = (c_p_coolant * visc_coolant) / k_coolant; % Prandtl number
Nu = 0.453 * sqrt(Re) * Pr^(1/3); % Nusselt relation for forced external
% convection for a flat plate. REDO THIS CONSIDERING THE CHARACTERISTIC
% LENGTH

% Calculate heat transfer coefficient
% REDO THIS WITH CHARACTERISTIC LENGTH
h = Nu * k_coolant / (T_coolant_out - T_coolant_in); % Heat transfer coefficient

% Conduction through coolant
Area_plate = length * height; % Area
T_coolant_core = T_coolant;
% In the entry point: 
% T_coolant_core = T_coolant_boundary = T_coolant_in
% Q_conduction = k_coolant * A * (T_coolant_core - T_coolant_boundary); % Heat transfer by conduction

% Boundary thickness of coolant: Assume it being 1/3rd of the interplate
% width on both sides.
boundary_coolant = 1/3 * width;
Conductive_Resistance_Coolant = boundary_coolant / ((Area_plate) * k_coolant);
Convective_Resistance_Coolant = 1 / (h * Area_plate);
% I assumed 0.5 mm to be the thickness of a steel plate. 
Conductive_Resistance_Steel = 0.0005 / ((Area_plate) * k_steel); 
