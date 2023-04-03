clc
clear


set(0,'DefaultTextFontSize',12);              % Settings for the plot
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultLegendInterpreter','Latex');
set(groot,'DefaultAxesTickLabelInterpreter','Latex');
set(0,'DefaultTextInterpreter','Latex');
set(0,'DefaultLegendFontSize',12);
T = 298.15;
P = [7:1:20];
data = nistdata('N2',T,P);
 
P_c = 7;                             % Pressure in the test chamber [bar]
d_p_ext = 12*1e-3;                   % Pipe external diameter [m]
t = 0.7*1e-3;                        % Thickness of the tube  [m]
d_p_int = d_p_ext - 2*t;             % Pipe internal diameter [m]
A_int = pi*(d_p_int/2)^2;            % Internal cross sectional area [m^2]
T_amb = 298.15;                      % Ambient temperature [K]
rho_N2 = data.Rho*data.Mw;           % Density of Nitrogen [kg/m^3] 
cp_N2 = data.Cp/data.Mw;             % Specific heat at constant pressure of Nitrogen [J/kgK]
cv_N2 = data.Cv/data.Mw;             % Specific heat at constant volume of Nitrogen [J/kgK]
gamma = cp_N2./cv_N2;                % Ratio of specific heats [-]
mu_N2 = data.mu;                     % Viscosity of Nitrogen [Pa*s]
L_test = 30*1e-3;                    % Dimension of the square test chamber [m]
rho = 1.832;                         % Density in full scale facility [kg/m^3]
u = 22.12;                           % Velocity in full scale facility [m/s]
mom_flux = rho*u^2;                  % Momentum flux in full scale facility [kg/ms^2]
m_dot_N2 = sqrt(mom_flux./rho_N2).*rho_N2*L_test^2; % Nitrogen mass flow rate [kg/s]
v_N2 = m_dot_N2./(A_int*rho_N2);     % Velocity of the fluid inside the pipes [m/s]
Re_N2 = 