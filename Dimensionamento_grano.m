clear
clc

% deltaHreag=-2129.73; %KW/kg dato dal CEA
% % per i prodotti tengo conto di CO,CO2,HCl,H20,N2 e H2 hanno 0
% deltaHprod= 0.32848*(-110.53/(28.01*1e-3))+0.10562*(-393.52/(44.01*1e-3))+0.24448*(-92.31/(36.458*1e-3))+0.20013*(-241.83/(18.02*1e-3)); %KW/kg
% deltaHreaz= deltaHprod-deltaHreag; %KW/kg
% Qcomb=-deltaHreaz*1e3; % W/kg

m_dot_N2 = 60*1e-3;          % kg/s
T_fin = 650;                 % K
T_amb = 298.15;              % K
cp_N2 = 1050;                % J/kgK
q_N2 = m_dot_N2*cp_N2*(T_fin - T_amb);   % W
cp_g = 2363;                 % J/kgK
T_fl = 2305;                 % K

f = @(x) m_dot_N2*cp_N2*(T_fin - T_amb) + x*cp_g*(T_fl) - (m_dot_N2 + x) * ( (cp_N2*m_dot_N2)/(m_dot_N2 + x) + (cp_g*x)/(m_dot_N2 + x) )*(T_fin);
m_dot_p = 1.2*fzero(f,0.5) ;


P = 1;     % bar
r_b = (1.152*P^0.768)*1e-3;  % m/s

m_AP = 0.8;     % kg
m_HTPB = 0.2;     % kg
rho_AP = 1950;   % kg/m^3
rho_HTPB = 920;    % kg/m^3
V_AP = m_AP/rho_AP; %m^3
V_HTPB = m_HTPB/rho_HTPB;  %m^3

rho_mix = (m_AP + m_HTPB)/(V_AP + V_HTPB);  % kg/m^3
t_b = 20;            % s

A_b = (m_dot_p/(r_b*rho_mix))*1e6;    % mm^2
d_b = 2*sqrt(A_b/pi);                 % mm

thickness = r_b*t_b*1e3;              % mm
m_prop=A_b*thickness*1e-9*rho_mix;
