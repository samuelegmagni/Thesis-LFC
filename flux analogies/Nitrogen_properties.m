clear
clc

format long

set(0,'DefaultTextFontSize',12);              % Settings for the plot
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultLegendInterpreter','Latex');
set(0,'DefaultAxesTickLabelInterpreter','Latex');
set(0,'DefaultTextInterpreter','Latex');
set(0,'DefaultLegendFontSize',12);

load('nitrogenThermoPhysicalProp.mat')

rho_N2 = rho;           % Density of Nitrogen [kg/m^3] 
cp_N2 = cp;             % Specific heat at constant pressure of Nitrogen [J/kgK]
cv_N2 = cv;             % Specific heat at constant volume of Nitrogen [J/kgK]
gamma_N2 = cp_N2./cv_N2;             % Ratio of specific heats [-]
mu_N2 = mu;                     % Viscosity of Nitrogen [Pa*s]
P = P(1,:)/1e5;
T = T(:,1)';

clear C; clear cp; clear Cp; clear cv; clear Cv; clear h; clear H; clear JT; clear mu; clear omega; clear pc; clear rho;
clear Rho; clear s; clear S; clear species; clear Tc; clear u;
clear U; clear V
%% Density

figure()
contourf(P,T,rho_N2); 
title('Nitrogen density')
xlabel('Pressure, P $[bar]$')
ylabel('Temperature, T $[K]$')
c = colorbar;
c.Label.String = 'Density, \rho [kg/m^3]';
%% Heat capacity at constant pressure

figure()
contourf(P,T,cp_N2); 
title('Nitrogen heat capacity')
xlabel('Pressure, P $[bar]$')
ylabel('Temperature, T $[K]$')
c = colorbar;
c.Label.String = 'Specific heat, c_p [J/Kg*K]';

%% Conductivity 

figure()
contourf(P,T,k); 
title('Nitrogen conductivity')
xlabel('Pressure, P $[bar]$')
ylabel('Temperature, T $[K]$')
c = colorbar;
c.Label.String = 'Thermal conductivity, \lambda [W/m*K]';

%% Viscosity 

figure()
contourf(P,T,mu_N2); 
title('Nitrogen viscosity')
xlabel('Pressure, P $[bar]$')
ylabel('Temperature, T $[K]$')
c = colorbar;
c.Label.String = 'Viscosity, \mu [Pa*s]';

