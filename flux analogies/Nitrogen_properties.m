clear
clc

format long

set(0,'DefaultTextFontSize',12);              % Settings for the plot
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultLegendInterpreter','Latex');
set(0,'DefaultAxesTickLabelInterpreter','Latex');
set(0,'DefaultTextInterpreter','Latex');
set(0,'DefaultLegendFontSize',12);

T = [300:5:850];
P = [1:5:200];

data = nistdata('N2',T,P); 
%% Density

figure()
contourf(data.P*1e-5,data.T,data.Rho*data.Mw); 
title('Density')
xlabel('Pressure $[bar]$')
ylabel('Temperature $[K]$')
c = colorbar;
c.Label.String = '\rho_{N_2} [kg/m^3]';
%% Heat capacity at constant pressure

figure()
contourf(data.P*1e-5,data.T,data.Cp/data.Mw); 
title('Heat capacity at constant pressure')
xlabel('Pressure $[bar]$')
ylabel('Temperature $[K]$')
c = colorbar;
c.Label.String = 'c_{p,N_2} [J/Kg*K]';

%% Conductivity 

figure()
contourf(data.P*1e-5,data.T,data.k); 
title('Conductivity')
xlabel('Pressure $[bar]$')
ylabel('Temperature $[K]$')
c = colorbar;
c.Label.String = '\lambda_{N_2} [W/m*K]';

%% Viscosity 

figure()
contourf(data.P*1e-5,data.T,data.mu); 
title('Viscosity')
xlabel('Pressure $[bar]$')
ylabel('Temperature $[K]$')
c = colorbar;
c.Label.String = '\mu_{N_2} [Pa*s]';

