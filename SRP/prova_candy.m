clear
clc


% NASA CEA data assumed at 4 bar 
m_dot_N2 = 60*1e-3;                     % kg/s
T_fin = 1.2*600;                        % K
T_beforeinj = 295;                      % K
cp_N2 = 1046;                         % J/kgK at 295 K and 12 bar
cp_N2_fin = 1104;                     % J/kgK at 600 K and 12 bar
cp_g = 5059.6;                        % J/kgK
T_fl = 0.7*1580;                     % K

f = @(x) m_dot_N2*cp_N2*T_beforeinj + x*cp_g*(T_fl) - (m_dot_N2 + x) * ( (cp_N2_fin*m_dot_N2 + cp_g*x)/(m_dot_N2 + x) )*(T_fin);
m_dot_p = fzero(f,0.5);


P = 4*14.504;                        % psi
r_b = (0.0665*P^0.319)/39.37;         % m/s

m_KNO3 = 0.65;                        % kg
m_Sucrose = 0.35;                     % kg
rho_KNO3 = 2110;                      % kg/m^3
rho_Sucrose = 1590;                   % kg/m^3
V_KNO3 = m_KNO3/rho_KNO3;             % m^3
V_Sucrose = m_Sucrose/rho_Sucrose;    % m^3

rho_mix = (m_KNO3 + m_Sucrose)/(V_KNO3 + V_Sucrose);  % kg/m^3
t_b = 20;                                             % s

A_b = (m_dot_p/(r_b*rho_mix))*1e6;                    % mm^2
d_b = 2*sqrt(A_b/pi);                                 % mm

thickness = r_b*t_b*1e3;                              % mm
m_prop=A_b*thickness*1e-9*rho_mix;                    % kg