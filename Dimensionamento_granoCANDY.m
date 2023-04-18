clear
clc

m_dot_N2 = 60*1e-3;          % kg/s
T_fin = 1.2*650;             % K
T_amb = 298.15;              % K
cp_N2 = 1052.08;             % J/kgK
cp_N2_fin = 1119.02;         % J/kgK
cp_g = 3830.4;               % J/kgK
T_fl = 0.7*1635.62;          % K

f = @(x) m_dot_N2*cp_N2*T_amb + x*cp_g*(T_fl) - (m_dot_N2 + x) * ( (cp_N2_fin*m_dot_N2 + cp_g*x)/(m_dot_N2 + x) )*(T_fin);
m_dot_p = fzero(f,0.5);


P = 10*14.504;     % psi
r_b = (0.0665*P^0.319)/39.37;  % m/s

m_KNO3 = 0.65;     % kg
m_Sucrose = 0.35;     % kg
rho_KNO3 = 2110;   % kg/m^3
rho_Sucrose = 1590;    % kg/m^3
V_KNO3 = m_KNO3/rho_KNO3; %m^3
V_Sucrose = m_Sucrose/rho_Sucrose;  %m^3

rho_mix = (m_KNO3 + m_Sucrose)/(V_KNO3 + V_Sucrose);  % kg/m^3
t_b = 20;            % s

A_b = (m_dot_p/(r_b*rho_mix))*1e6;    % mm^2
d_b = 2*sqrt(A_b/pi);                % mm

thickness = r_b*t_b*1e3;              % mm
m_prop=A_b*thickness*1e-9*rho_mix;

%%

set(0,'DefaultTextFontSize',12);              % Settings for the plot
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultLegendInterpreter','Latex');
set(groot,'DefaultAxesTickLabelInterpreter','Latex');
set(0,'DefaultTextInterpreter','Latex');
set(0,'DefaultLegendFontSize',12);

T = 1.2*[300:5:650];
P = [1:0.2:10];

T_amb = 298.15;
cp_g = 3830.4;                 % J/kgK
T_fl = 0.7*1635.62; 

data = nistdata('N2',T,P); 

L_test = 30*1e-3;
rho = 1.832;
u = 22.12;
mom_flux = rho*u^2;      % mom_flux_paper = 896.89;

u_exp = sqrt(mom_flux./(data.Rho*data.Mw));

m_dot_N2 = data.Rho*data.Mw.*u_exp*L_test^2;

cp_N2 = data.Cp/data.Mw;

m_dot_p = zeros(length(T),length(P));

for i = 1 : length(T)

    for m = 1 : length(P)

        f = @(x) m_dot_N2(i,m)*cp_N2(1,m)*T_amb + x*cp_g*(T_fl) - (m_dot_N2(i,m) + x) * ( (cp_N2(end,m)*m_dot_N2(i,m) + (cp_g*x))/(m_dot_N2(i,m) + x) )*(T(i));
        z = fzero(f,0.05) ;

        m_dot_p(i,m) = z; 

    end

end

figure()
contourf(data.P*1e-5,data.T,m_dot_N2*1e3); 
title('Slab 30x30: momentum flux analogy')
xlabel('Pressure $[bar]$')
ylabel('Temperature $[K]$')
c = colorbar;
c.Label.Interpreter = 'latex';
c.Label.String = 'Nitrogen mass flow rate $\dot{m}_{N_2}$ [g/s]';

figure()
contourf(data.P*1e-5,data.T,m_dot_p*1e3); 
title('Slab 30x30: momentum flux analogy')
xlabel('Pressure $[bar]$')
ylabel('Temperature $[K]$')
c = colorbar;
c.Label.Interpreter = 'latex';
c.Label.String = 'Propellant mass flow rate $\dot{m}_{prop}$ [g/s]';


figure()
contourf(data.P*1e-5,data.T,d_b); 
title('Slab 30x30: momentum flux analogy')
xlabel('Pressure $[bar]$')
ylabel('Temperature $[K]$')
c = colorbar;
c.Label.Interpreter = 'latex';
c.Label.String = 'SRM charge diameter $d_b$ [mm]';

figure()
contourf(data.P*1e-5,data.T,m_prop); 
title('Slab 30x30: momentum flux analogy')
xlabel('Pressure $[bar]$')
ylabel('Temperature $[K]$')
c = colorbar;
c.Label.Interpreter = 'latex';
c.Label.String = 'SRM charge mass $m_{prop}$ [kg]';