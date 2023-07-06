clc
clear


set(0,'DefaultTextFontSize',12);              % Settings for the plot
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultLegendInterpreter','Latex');
set(groot,'DefaultAxesTickLabelInterpreter','Latex');
set(0,'DefaultTextInterpreter','Latex');
set(0,'DefaultLegendFontSize',12);

%% After pressure regulator (point 1)

T_chamber = 600;                                       % Temperature downstream the pressure regulator [K]
P_chamber = 2;  

load('nitrogenThermoPhysicalProp.mat')

rho_N2 = rho;           % Density of Nitrogen [kg/m^3] 
cp_N2 = cp;             % Specific heat at constant pressure of Nitrogen [J/kgK]
cv_N2 = cv;             % Specific heat at constant volume of Nitrogen [J/kgK]
gamma_N2 = cp_N2./cv_N2;             % Ratio of specific heats [-]
mu_N2 = mu;                     % Viscosity of Nitrogen [Pa*s]
P = P(1,:)/1e5;
T = T(:,1)';

clear C; clear cp; clear Cp; clear cv; clear Cv; clear h; clear H; clear JT;
clear k; clear mu; clear Mw; clear omega; clear pc; clear rho;
clear Rho; clear s; clear S; clear species; clear Tc; clear u;
clear U; clear V

m_dot_N2 = 20*1e-3;                  % Nitrogen mass flow rate [kg/s]
R_N2 = 8314/28;                         % Specific ideal gas constant [J/kgK]

R_SRP = 8314/36.32;
m_dot_SRP = 18*1e-3;
gamma_SRP = 1.1122;
rho_SRP = 0.92189;

R_mix = (m_dot_SRP*R_SRP + m_dot_N2*R_N2)/(m_dot_N2 + m_dot_SRP);

m_dot_mix = m_dot_N2 + m_dot_SRP;

rho_chamber_N2 = rho_N2(find(abs(T - round(T_chamber,1))==min(abs(T - round(T_chamber,1)))) ,find( abs(P - round(P_chamber,1))==min(abs(P - round(P_chamber,1)))) );
rho_chamber = (m_dot_SRP*rho_SRP + m_dot_N2*rho_chamber_N2)/(m_dot_N2 + m_dot_SRP);

gamma_chamber_N2 = gamma_N2(find(abs(T - round(T_chamber,1))==min(abs(T - round(T_chamber,1)))) ,find( abs(P - round(P_chamber,1))==min(abs(P - round(P_chamber,1)))) );
gamma_chamber = (m_dot_SRP*gamma_SRP + m_dot_N2*gamma_chamber_N2)/(m_dot_N2 + m_dot_SRP);

gamma30_N2 = gamma_N2(find(abs(T - round(T_chamber,1))==min(abs(T - round(T_chamber,1)))) ,find( abs(P - round(P_chamber,1))==min(abs(P - round(P_chamber,1)))) );
gamma30 = (m_dot_SRP*gamma_SRP + m_dot_N2*gamma30_N2)/(m_dot_N2 + m_dot_SRP);


d_inj = 11.9*1e-3;    % [m]
A_throat_int = 0.25*pi*d_inj^2;

M_throat = 1;

P_tot = P_chamber*((1 + 0.5*(gamma_chamber - 1)*M_throat^2)^(gamma_chamber/(gamma_chamber - 1)));

d29_30_ext = 19.05*1e-3;
t29_30 = 1.5*1e-3;
d29_30_int = d29_30_ext - 2*t29_30;
A29_30 = pi*(d29_30_int/2)^2;

eps = 0.015*1e-3;
eps29_30_rel = eps/d29_30_int;               % Relative roughness of stainless steel [-]

z = @(x) A_throat_int/A29_30 - (x/M_throat)*sqrt( ((1 + 0.5*(gamma_chamber - 1)*M_throat^2)/(1 + 0.5*(gamma_chamber - 1)*x^2))^((gamma_chamber + 1)/(gamma_chamber - 1)) );

M30 = fsolve(z,0.08)

iter = 0;
err = 1;

while err > 1e-3
    iter = iter + 1;

    P30 = P_tot/((1 + 0.5*(gamma30 - 1)*M30^2)^(gamma30/(gamma30 - 1)));

    T_tot = T_chamber*(1 + 0.5*(gamma_chamber - 1)*M_throat^2);

    T30 = T_tot/(1 + 0.5*(gamma30 - 1)*M30^2);

    gamma30_N2_new = gamma_N2(find(abs(T - round(T30,1))==min(abs(T - round(T30,1)))) ,find( abs(P - round(P30,1))==min(abs(P - round(P30,1)))) ); 
    gamma30_new = (m_dot_SRP*gamma_SRP + m_dot_N2*gamma30_N2_new)/(m_dot_N2 + m_dot_SRP);

    err = abs(gamma30 - gamma30_new);

    gamma30 = gamma30_new;

end

clear gamma30_new
clear gamma30_N2_new
