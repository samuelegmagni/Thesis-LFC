clear
clc

% deltaHreag=-2129.73; %KW/kg dato dal CEA
% % per i prodotti tengo conto di CO,CO2,HCl,H20,N2 e H2 hanno 0
% deltaHprod= 0.32848*(-110.53/(28.01*1e-3))+0.10562*(-393.52/(44.01*1e-3))+0.24448*(-92.31/(36.458*1e-3))+0.20013*(-241.83/(18.02*1e-3)); %KW/kg
% deltaHreaz= deltaHprod-deltaHreag; %KW/kg
% Qcomb=-deltaHreaz*1e3; % W/kg


m_dot_N2 = 60*1e-3;          % kg/s
T_fin = 1.2*650;             % K
T_amb = 298.15;              % K
cp_N2 = 1050;                % J/kgK
q_N2 = m_dot_N2*cp_N2*(T_fin - T_amb);   % W
cp_g = 2363;                 % J/kgK
T_fl = 0.8*2305;             % K

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

%%

clear
clc
set(0,'DefaultTextFontSize',12);              % Settings for the plot
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultLegendInterpreter','Latex');
set(groot,'DefaultAxesTickLabelInterpreter','Latex');
set(0,'DefaultTextInterpreter','Latex');
set(0,'DefaultLegendFontSize',12);

T = 1.2*[300:5:650];
P = [1:0.2:10];

T_amb = 298.15;
cp_g = 2363;                 % J/kgK
T_fl = 0.8*2305; 

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

        f = @(x) m_dot_N2(i,m)*cp_N2(i,m)*(T(i) - T_amb) + x*cp_g*(T_fl) - (m_dot_N2(i,m) + x) * ( (cp_N2(i,m)*m_dot_N2(i,m))/(m_dot_N2(i,m) + x) + (cp_g*x)/(m_dot_N2(i,m) + x) )*(T(i));
        z = fsolve(f,0.05) ;
        m_dot_p(i,m) = z; 

    end

end

figure()
contourf(data.P*1e-5,data.T,m_dot_N2*1e3); 
title('Slab 30x30: mass flow rate momentum flux analogy')
xlabel('Pressure $[bar]$')
ylabel('Temperature $[K]$')
c = colorbar;
c.Label.String = 'm_{dot,N_2} [g/s]';

figure()
contourf(data.P*1e-5,data.T,m_dot_p*1e3); 
title('Slab 30x30: mass flow rate momentum flux analogy')
xlabel('Pressure $[bar]$')
ylabel('Temperature $[K]$')
c = colorbar;
c.Label.String = 'm_{dot,p} [g/s]';

P_prop = 1;     % bar
r_b = (1.152*P_prop^0.768)*1e-3;  % m/s

m_AP = 0.8;     % kg
m_HTPB = 0.2;     % kg
rho_AP = 1950;   % kg/m^3
rho_HTPB = 920;    % kg/m^3
V_AP = m_AP/rho_AP; %m^3
V_HTPB = m_HTPB/rho_HTPB;  %m^3

rho_mix = (m_AP + m_HTPB)/(V_AP + V_HTPB);  % kg/m^3
t_b = 20;            % s

A_b = (m_dot_p./(r_b*rho_mix))*1e6;    % mm^2
d_b = 2*sqrt(A_b/pi);                 % mm

thickness = r_b*t_b*1e3;              % mm
m_prop=A_b*thickness*1e-9*rho_mix;

figure()
contourf(data.P*1e-5,data.T,d_b); 
title('Slab 30x30: mass flow rate momentum flux analogy')
xlabel('Pressure $[bar]$')
ylabel('Temperature $[K]$')
c = colorbar;
c.Label.String = 'd_b [mm]';

figure()
contourf(data.P*1e-5,data.T,m_prop); 
title('Slab 30x30: mass flow rate momentum flux analogy')
xlabel('Pressure $[bar]$')
ylabel('Temperature $[K]$')
c = colorbar;
c.Label.String = 'm_{prop} [kg]';
