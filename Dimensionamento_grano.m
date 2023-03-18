clear
clc

P = 1;     % bar
r_b = (1.152*P^0.768)*1e-3;  % m/s

m_AP = 0.8;     % kg
m_HTPB = 0.2;     % kg
rho_AP = 1950;   % kg/m^3
rho_HTPB = 920;    % kg/m^3
V_AP = m_AP/rho_AP; %m^3
V_HTPB = m_HTPB/rho_HTPB;  %m^3

rho_mix = (m_AP + m_HTPB)/(V_AP + V_HTPB);  % kg/m^3

m_dot_p = 0.0039   % kg/s
t_b = 20;            % s

A_b = (m_dot_p/(r_b*rho_mix))*1e6    % mm^2
d_b = 2*sqrt(A_b/pi)                 % mm

thickness = r_b*t_b*1e3              % mm

m_dot_N2 = 50*1e-3;          % kg/s
T_fin = 650;                 % K
T_amb = 298.15;              % K
cp_N2 = 1050;                % J/kgK
q_N2 = m_dot_N2*cp_N2*(T_fin - T_amb);   % W
cp_g = 2363;                 % J/kgK
T_fl = 2305;                 % K
m_dot_p = q_N2/(cp_g*(T_fl - T_amb));   % kg/s