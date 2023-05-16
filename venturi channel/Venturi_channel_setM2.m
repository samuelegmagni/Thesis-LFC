%% Before convergent
clear
clc

d_ext = 12*1e-3;
t = 1.5*1e-3;                        % Thickness of the tube  [m]
d1 = d_ext - 2*t;
A1 = 0.25*pi*d1^2;
m_dot_N2 = 140*1e-3;       % [g/s]
P1 = 39;                   % [bar]
T1 = 295;                  % [K]
R = 8314/28;               % [J/KgK]
eps = 0.015*1e-3;               % Absolute roughness of stainless steel [m]
eps_rel = eps/d1;               % Relative roughness of stainless steel [-]


T = (floor(T1)-3):0.5:(ceil(T1));
P = (floor(P1)-2):0.1:(ceil(P1));

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

%% Convergent part

T_tot = T1*(1 + ((gamma1 - 1)/2)*M1^2);
P_tot2 = P1*(1 + ((gamma1 - 1)/2)*M1^2)^(gamma1/(gamma1 - 1));
M2 = 1;
z = @(x) x/A1 - (M1/M2)*sqrt( ((1 + 0.5*(gamma1 - 1)*M2^2)/(1 + 0.5*(gamma1 - 1)*M1^2))^((gamma1 + 1)/(gamma1 - 1)) );
A2 = fsolve(z,0.8);
d2 = sqrt((4*A2)/pi);
T2 = T_tot/(1 + ((gamma1 - 1)/2)*M2^2);
P2 = P_tot2/(1 + ((gamma1 - 1)/2)*M2^2)^(gamma1/(gamma1 - 1));
c2 = sqrt(gamma1*R*T2);
v2 = c2;
rho2 = (rho1*v1*A1)/(A2*v2);

%%
T = (floor(T2)):0.5:(ceil(T2)+15);
P = (floor(P2)):0.1:(ceil(P2)+7);

data = nistdata('N2',T,P);

rho_N2 = data.Rho*data.Mw;           % Density of Nitrogen [kg/m^3] 
cp_N2 = data.Cp/data.Mw;             % Specific heat at constant pressure of Nitrogen [J/kgK]
cv_N2 = data.Cv/data.Mw;             % Specific heat at constant volume of Nitrogen [J/kgK]
gamma_N2 = cp_N2./cv_N2;             % Ratio of specific heats [-]
mu_N2 = data.mu;                     % Viscosity of Nitrogen [Pa*s]

v2 = (rho1*A1*v1/(A2*rho2));

mu2 = mu_N2(find(T==round(T2)),find(abs(P - round(P2,1)) < 0.001));
gamma2 = gamma_N2(find(T==round(T2)),find(abs(P - round(P2,1)) < 0.001));
gamma3 = gamma_N2(find(T==round(T2)),find(abs(P - round(P2,1)) < 0.001));
Re2 = (rho2*v2*d2)/mu2;

if Re2 < 2300

        lambda = 64/Re2;

    else 

        z = @(x) 1/sqrt(x) + 2*log10(2.51/(Re2*sqrt(x)) + eps_rel/3.71);   % Colebrook-White correlation
        lambda = fsolve(z,0.0004);

end

%%

L = 0.08;
iter = 0;
err = 1;

while err > 1e-3

    iter = iter + 1;

    g_M2 = (1 - M2^2)/(gamma2*M2^2) + ((gamma2 + 1)/(2*gamma2))*log(((gamma2 + 1)*M2^2)/(2 + (gamma2 - 1)*M2^2) );
    g_M3 = abs(g_M2 - (lambda/d2)*L);

    y = @(x) g_M3 - (1 - x^2)/(gamma3*x^2) + ((gamma3 + 1)/(2*gamma3))*log(((gamma3 + 1)*x^2)/(2 + (gamma3 - 1)*x^2) );
    M3 = fsolve(y,0.6);

    T_star = T2/(0.5*(gamma2 + 1)/(1 + (gamma2 - 1)*0.5*M2^2));
    T3 = T_star*(0.5*(gamma3 + 1)/(1 + (gamma3 - 1)*0.5*M3^2));
  
    P_star = P2/((1/M2)*sqrt(0.5*(gamma2 + 1)/(1 + (gamma2 - 1)*0.5*M2^2)));
    P3 = P_star*((1/M3)*sqrt(0.5*(gamma3 + 1)/(1 + (gamma3 - 1)*0.5*M3^2)));

    P_tot_star = (M2*P_tot2)/( ((1 + 0.5*(gamma2 - 1)*M2^2)/(0.5*(gamma2 + 1)))^( (0.5*(gamma2 + 1))/(gamma2 - 1) ));
    P_tot3 = (P_tot_star/M3)*( ((1 + 0.5*(gamma3 - 1)*M3^2)/(0.5*(gamma3 + 1)))^( (0.5*(gamma3 + 1))/(gamma3 - 1) ));

    rho_star = rho2/((1/M2)*sqrt( 2*(1 + (gamma2 - 1)*0.5*M2^2)/(gamma2 + 1)));
    rho3 = rho_star*((1/M3)*sqrt( 2*(1 + (gamma3 - 1)*0.5*M3^2)/(gamma3 + 1)));

    c3 = sqrt(gamma3*R*T3);
    v3 = M3*c3;

    gamma3_new = gamma_N2(find(T==round(T3)),find(abs(P - round(P3,1)) < 0.001));

    err = abs(gamma3 - gamma3_new);

    gamma3 = gamma3_new;

end

clear gamma3_new

A3 = A2;
A4 = A1;

z = @(x) A4/A3 - (M3/x)*sqrt( ((1 + 0.5*(gamma3 - 1)*x^2)/(1 + 0.5*(gamma3 - 1)*M3^2))^((gamma3 + 1)/(gamma3 - 1)) );
M4 = fsolve(z,0.8);

P4 = P_tot3/(1 + ((gamma3 - 1)/2)*M4^2)^(gamma3/(gamma3 - 1));
T4 = T_tot/(1 + ((gamma3 - 1)/2)*M4^2);
c4 = sqrt(gamma3*R*T3);
v4 = M4*c4;
rho4 = (rho3*v3*A3)/(A4*v4);