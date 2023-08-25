clear
clc

% all data calculated at 4 bar 

m_dot_N2 = 60*1e-3;          % kg/s
T_fin = 1.2*600;             % K
T_amb = 295;              % K
cp_N2_in = 1046;                % J/kgK
cp_N2_fin = 1104;          % J/kgK
cp_g = 2159;                 % J/kgK
T_fl = 0.8*2325;             % K

f = @(x) m_dot_N2*cp_N2_in*T_amb + x*cp_g*(T_fl) - (m_dot_N2 + x) * ( (cp_N2_fin*m_dot_N2 + cp_g*x)/(m_dot_N2 + x) )*(T_fin);
m_dot_p = fzero(f,0.5) 


P = 4;     % bar
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
d_b = 2*sqrt(A_b/pi)                 % mm

thickness = r_b*t_b*1e3              % mm
m_prop=A_b*thickness*1e-9*rho_mix