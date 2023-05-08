clear
clc

format long

set(0,'DefaultTextFontSize',12);              % Settings for the plot
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultLegendInterpreter','Latex');
set(0,'DefaultAxesTickLabelInterpreter','Latex');
set(0,'DefaultTextInterpreter','Latex');
set(0,'DefaultLegendFontSize',12);

%%
T_vect = [280:5:300];
P_vect = [1:0.5:10];

data = nistdata('N2',T_vect,P_vect);

%% GE-C parameters

R_CC = 79.2 * 1e-3;
R_t = 14.42 * 1e-3;
eps_c = 30.16;
eps_e = 125;
L_CC = 22*1e-2;
m_dot = 0.798;
rho = 1.832;
T = 2787;
u = 22.12;
M_m = 21.23; %g/mol
c_p = 2519;  %J/kgK
mu = 9.76*1e-5;
lambda = 0.3372;

A_tot = pi*2*R_CC*L_CC + 2*pi*R_CC^2;

%% Full scale test section: Re similarity

D_test = 2*R_CC;

Re = rho*u*2*R_CC/mu;     %Re_paper = 65797

u_exp = Re*(data.mu./( data.Rho*data.Mw*D_test ));

m_dot_exp = data.Rho*data.Mw.*u_exp*0.25*pi*D_test^2;

adjust_vect = ones(length(T_vect),length(P_vect));

T_amb = 298;

q_dot_exp = m_dot_exp.*( (data.Cp/data.Mw).*T_vect' - ( (data.Cp(1,1)/data.Mw)*T_amb*adjust_vect ) );

figure()
contourf(data.P*1e-5,data.T,m_dot_exp); 
shading flat;
title('Mass flow rate $Re_D$ similarity')
xlabel('Pressure $[bar]$')
ylabel('Temperature $[K]$')
c = colorbar;
c.Label.Interpreter = 'latex';
c.Label.String = '$\dot{m}_{N_2}$ [kg/s]';

figure()
contourf(data.P*1e-5,data.T,1e-3*q_dot_exp); 
shading flat;
title('Energy $Re_D$ similarity')
xlabel('Pressure $[bar]$')
ylabel('Temperature $[K]$')
c = colorbar;
c.Label.Interpreter = 'latex';
c.Label.String = '$\dot{q}_{N_2}$ [kW]';

%%
clear
clc

d_ext = 12*1e-3;
t = 1.5*1e-3;                        % Thickness of the tube  [m]
d1 = d_ext - 2*t;
A1 = 0.25*pi*d1^2;
m_dot_N2 = 200*1e-3;       % [g/s]
P1 = 56;                   % [bar]
T1 = 295;                  % [K]
R = 8314/28;               % [J/KgK]
eps = 0.015*1e-3;                    % Absolute roughness of stainless steel [m]
eps_rel = eps/d1;               % Relative roughness of stainless steel [-]


T = (floor(T1)-30):0.5:(ceil(T1));
P = (floor(P1)-12):0.1:(ceil(P1));

data = nistdata('N2',T,P);

rho_N2 = data.Rho*data.Mw;           % Density of Nitrogen [kg/m^3] 
cp_N2 = data.Cp/data.Mw;             % Specific heat at constant pressure of Nitrogen [J/kgK]
cv_N2 = data.Cv/data.Mw;             % Specific heat at constant volume of Nitrogen [J/kgK]
gamma_N2 = cp_N2./cv_N2;             % Ratio of specific heats [-]
mu_N2 = data.mu;                     % Viscosity of Nitrogen [Pa*s]


rho1 = rho_N2(find(T==round(T1)),find(abs(P - round(P1,1)) < 0.001));
mu1 = mu_N2(find(T==round(T1)),find(abs(P - round(P1,1)) < 0.001));
gamma1 = gamma_N2(find(T==round(T1)),find(abs(P - round(P1,1)) < 0.001));

v1 = m_dot_N2/(rho1*A1);
c1 = sqrt(gamma1*R*T1);
M1 = v1/c1; 
Re1 = (rho1*v1*d1)/mu1;              % Alto quindi effetti inerziali dominanti rispetto a quelli viscosi

d2 = 7.5*1e-3;                        % Throat diameter [m]
A2 = 0.25*pi*d2^2;                  % Throat cross sectional area [m^2]
beta = sqrt(A2/A1);
C = 0.995;                          % Disharge coefficient
k = gamma1;
%%
z = @(x) m_dot_N2 - sqrt( ((x/P1)^(2/k)) * (k/(k-1)) * ((1 - (x/P1)^((k-1)/k))/(1 - x/P1)) * ((1 - beta^4)/(1 - beta^4*((x/P1)^(2/k)))))*C*A2*sqrt(2*rho1*(P1 - x)/(1 - (d2/d1)^4));
P2 = fsolve(z,0)
%%
T2 = T1;
rho2 = rho_N2(find(T==round(T1)),find(abs(P - round(P2,1)) < 0.001));
mu2 = mu_N2(find(T==round(T1)),find(abs(P - round(P2,1)) < 0.001));
gamma2 = gamma_N2(find(T==round(T2)),find(abs(P - round(P2,1)) < 0.001));
gamma3 = gamma_N2(find(T==round(T2)),find(abs(P - round(P2,1)) < 0.001));
L = 0.05;                            % Length of the throat [m]
v2 = m_dot_N2/(rho2*A2);
c2 = sqrt(gamma1*R*T1);
M2 = v2/c2;
Re2 = (rho2*v2*d2)/mu2;

if Re2 < 2300

        lambda = 64/Re2;

    else 

        z = @(x) 1/sqrt(x) + 2*log10(2.51/(Re2*sqrt(x)) + eps_rel/3.71);   % Colebrook-White correlation
        lambda = fsolve(z,0.0004);

end

iter = 0;
err = 1;

while err > 1e-3

    iter = iter + 1;

    g_M2 = (1 - M2^2)/(gamma2*M2^2) + ((gamma2 + 1)/(2*gamma2))*log(((gamma2 + 1)*M2^2)/(2 + (gamma2 - 1)*M2^2) );
    g_M3 = g_M2 - lambda*(L/d2);
    
    y = @(x) g_M3 - (1 - x^2)/(gamma3*x^2) + ((gamma3 + 1)/(2*gamma3))*log(((gamma3 + 1)*x^2)/(2 + (gamma3 - 1)*x^2) );
    M3 = fsolve(y,0.006);
    
    T_star = T2/(0.5*(gamma2 + 1)/(1 + (gamma2 - 1)*0.5*M2^2));
    T3 = T_star*(0.5*(gamma3 + 1)/(1 + (gamma3 - 1)*0.5*M3^2));
    
    P_star = P2/((1/M2)*sqrt(0.5*(gamma2 + 1)/(1 + (gamma2 - 1)*0.5*M2^2)));
    P3 = P_star*((1/M3)*sqrt(0.5*(gamma3 + 1)/(1 + (gamma3 - 1)*0.5*M3^2)));
    
    rho_star = rho2/((1/M2)*sqrt( 2*(1 + (gamma2 - 1)*0.5*M2^2)/(gamma2 + 1)));
    rho3 = rho_star*((1/M3)*sqrt( 2*(1 + (gamma3 - 1)*0.5*M3^2)/(gamma3 + 1)));
    
    c3 = sqrt(gamma3*R*T3);
    v3 = M3*c3;

    gamma3_new = gamma_N2(find(T==round(T3)),find(abs(P - round(P3,1)) < 0.001));

    err = abs(gamma3 - gamma3_new);

    gamma3 = gamma3_new;

end

clear gamma3_new

