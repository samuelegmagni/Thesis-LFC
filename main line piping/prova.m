clear
clc

set(0,'DefaultTextFontSize',12);              % Settings for the plot
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultLegendInterpreter','Latex');
set(groot,'DefaultAxesTickLabelInterpreter','Latex');
set(0,'DefaultTextInterpreter','Latex');
set(0,'DefaultLegendFontSize',12);

T13 = 600;
P13 = 1.5;                           % Temperature and pressure at the exit of the injector

T = (floor(T13)-10):0.5:(ceil(T13));
P = (floor(P13)):0.1:(ceil(P13)+10);
data = nistdata('N2',T,P);

R = 8314/28;

rho_N2 = data.Rho*data.Mw;           % Density of Nitrogen [kg/m^3] 
cp_N2 = data.Cp/data.Mw;             % Specific heat at constant pressure of Nitrogen [J/kgK]
cv_N2 = data.Cv/data.Mw;             % Specific heat at constant volume of Nitrogen [J/kgK]
gamma_N2 = cp_N2./cv_N2;             % Ratio of specific heats [-]
mu_N2 = data.mu;  

rho13 = rho_N2(find(T==round(T13)),find(abs(P - round(P13,1)) < 0.001));     % Density downstream the check valve [kg/m^3]
gamma13 = gamma_N2(find(T==round(T13)),find(abs(P - round(P13,1)) < 0.001)); % Ratio of specific heats the check valve  [-]
gamma12 = gamma_N2(find(T==round(T13)),find(abs(P - round(P13,1)) < 0.001));

M13 = 1; 
m_dot_new = 78*1e-3;
d12 = 16.05*1e-3;
A12 = 0.25*pi*d12^2;
C_d = 0.61;
c13 = sqrt(gamma13*R*T13);
v_inj = c13

deltaP_inj = (rho13*v_inj^2)

A_inj = m_dot_new/(C_d*sqrt(2*deltaP_inj*rho13));

d_inj = sqrt((4*A_inj)/pi)

z = @(x) A12/x - (M13/M12)