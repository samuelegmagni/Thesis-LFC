clear
clc

set(0,'DefaultTextFontSize',12);              % Settings for the plot
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultLegendInterpreter','Latex');
set(groot,'DefaultAxesTickLabelInterpreter','Latex');
set(0,'DefaultTextInterpreter','Latex');
set(0,'DefaultLegendFontSize',12);

T13 = 600;
P13 = 1.5;                           % Temperature and pressure at the exit of the injector

T = (floor(T13)):0.5:(ceil(T13)+80);
P = (floor(P13)):0.1:(ceil(P13)+5);
data = nistdata('N2',T,P);

eps = 0.015*1e-3;                    % Absolute roughness of stainless steel [m]
eps_rel = eps/d_p_int;               % Relative roughness of stainless steel [-]

R_N2 = 8314/28;
R_SRP = 8314/36.32;
d_p_int = 16.05*1e-3;
A12 = 0.25*pi*d_p_int^2;

m_dot_SRP = 18*1e-3;
m_dot_N2 = 60*1e-3;
m_dot12 = m_dot_N2 + m_dot_SRP;

rho_N2 = data.Rho*data.Mw;           % Density of Nitrogen [kg/m^3] 
cp_N2 = data.Cp/data.Mw;             % Specific heat at constant pressure of Nitrogen [J/kgK]
cv_N2 = data.Cv/data.Mw;             % Specific heat at constant volume of Nitrogen [J/kgK]
gamma_N2 = cp_N2./cv_N2;             % Ratio of specific heats [-]
mu_N2 = data.mu;  

gamma_N2_13 = gamma_N2(find(T==round(T13)),find(abs(P - round(P13,1)) < 0.001)); % Ratio of specific heats the check valve  [-]

gamma_SRP = 1.1122;
gamma13 = (m_dot_SRP*gamma_SRP + m_dot_N2*gamma_N2_13)/(m_dot_N2 + m_dot_SRP);

rho_SRP = 0.92189;
rho_N2_13 = rho_N2(find(T==round(T13)),find(abs(P - round(P13,1)) < 0.001));
rho13 = (m_dot_SRP*rho_SRP + m_dot_N2*rho_N2_13)/(m_dot_N2 + m_dot_SRP);

R12 = (m_dot_SRP*R_SRP + m_dot_N2*R_N2)/(m_dot_N2 + m_dot_SRP);
c13 = sqrt(gamma13*R12*T13); 
v_inj = c13;

C_d = 0.61;
deltaP_inj = (0.5*rho13*v_inj^2)/C_d^2;
P12 = P13 + deltaP_inj*1e-5;

A_inj = m_dot12/(v_inj*rho13);
d_inj = sqrt((4/pi)*A_inj);

M13 = 1;
T_tot13 = T13*(1 + 0.5*(gamma13 - 1)*M13^2);

gamma_N2_12 = gamma_N2(find(T==round(T13)),find(abs(P - round(P12,1)) < 0.001));
gamma12 = (m_dot_SRP*gamma_SRP + m_dot_N2*gamma_N2_12)/(m_dot_N2 + m_dot_SRP);
 
rho_N2_12 = rho_N2(find(T==round(T13)),find(abs(P - round(P12,1)) < 0.001));
rho12 = (m_dot_SRP*rho_SRP + m_dot_N2*rho_N2_12)/(m_dot_N2 + m_dot_SRP);

iter = 0;
err = 1;

T12 = T13;

while err > 1e-3

    iter = iter + 1;

    c12 = sqrt(gamma12*R12*T12);

    v12 = m_dot12/(A_inj*rho12);

    M12 = v12/c12;

    T12 = T_tot13/(1 + 0.5*(gamma12 - 1)*M12^2);

    gamma_N2_12_new = gamma_N2(find(T==round(T12)),find(abs(P - round(P12,1)) < 0.001));
    gamma12_new = (m_dot_SRP*gamma_SRP + m_dot_N2*gamma_N2_12)/(m_dot_N2 + m_dot_SRP);

    rho_N2_12_new = rho_N2(find(T==round(T12)),find(abs(P - round(P12,1)) < 0.001));
    rho12_new = (m_dot_SRP*rho_SRP + m_dot_N2*rho_N2_12)/(m_dot_N2 + m_dot_SRP);

    c12_new = sqrt(gamma12_new*R12*T12);

    v12_new = m_dot12/(A_inj*rho12_new);

    M12_new = v12_new/c12_new;

    T12_new = T_tot13/(1 + 0.5*(gamma12_new - 1)*M12_new^2);

    err = abs(T12 - T12_new);

    c12 = c12_new;
    v12 = v12_new;
    M12 = M12_new;
    T12 = T12_new;

end

clear c12_new
clear v12_new
clear M12_new
clear T12_new
clear gamma12_new
clear rho12_new

%% Mixing chamber pressure's evolution

P_tank = P12*1e5 + 0.5*rho12*v12^2; 

T11 = 300;

T = (floor(T11)-10):0.5:(ceil(T11));
P = (floor(P12)):0.1:(ceil(P12)+5);
data = nistdata('N2',T,P);

rho_N2 = data.Rho*data.Mw;           % Density of Nitrogen [kg/m^3] 
cp_N2 = data.Cp/data.Mw;             % Specific heat at constant pressure of Nitrogen [J/kgK]
cv_N2 = data.Cv/data.Mw;             % Specific heat at constant volume of Nitrogen [J/kgK]
gamma_N2 = cp_N2./cv_N2;             % Ratio of specific heats [-]
mu_N2 = data.mu;

rho11 = rho_N2(find(T==round(T11)),find(abs(P - round(P12,1)) < 0.001));

v11 = m_dot_N2/(A12*rho11);

P11 = (P_tank + rho11*v11^2)*1e-5;

err = 1;
iter = 0;

while err > 1e-3
   
    rho11_new = rho_N2(find(T==round(T11)),find(abs(P - round(P11,1)) < 0.001));

    v11_new = m_dot_N2/(A12*rho11_new);

    P11_new = (P_tank + rho11_new*v11_new^2)*1e-5; 

    err = abs(rho11 - rho11_new);

    rho11 = rho11_new;
    v11 = v11_new;
    P11 = P11_new;

end

clear rho11_new
clear v11_new
clear P11_new

gamma11 = gamma_N2(find(T==round(T11)),find(abs(P - round(P11,1)) < 0.001));

c11 = sqrt(gamma11*R_N2*T11);
M11 = v11/c11;

%% Before mixing chamber 


T = (floor(T11)):0.5:(ceil(T11)+30);
P = (floor(P11)):0.1:(ceil(P11)+10);
data = nistdata('N2',T,P);

rho_N2 = data.Rho*data.Mw;           % Density of Nitrogen [kg/m^3] 
cp_N2 = data.Cp/data.Mw;             % Specific heat at constant pressure of Nitrogen [J/kgK]
cv_N2 = data.Cv/data.Mw;             % Specific heat at constant volume of Nitrogen [J/kgK]
gamma_N2 = cp_N2./cv_N2;             % Ratio of specific heats [-]
mu_N2 = data.mu;  

rho11 = rho_N2(find(T==round(T11)),find(abs(P - round(P11,1)) < 0.001));     
gamma11 = gamma_N2(find(T==round(T11)),find(abs(P - round(P11,1)) < 0.001)); 
gamma10 = gamma_N2(find(T==round(T11)),find(abs(P - round(P11,1)) < 0.001));
mu11 = mu_N2(find(T==round(T11)),find(abs(P - round(P11,1)) < 0.001)); 
Re11 = (rho11*v11*d_p_int)/mu11;              

if Re11 < 2300

        lambda = 64/Re11;

    else 

        z = @(x) 1/sqrt(x) + 2*log10(2.51/(Re11*sqrt(x)) + eps_rel/3.71);   % Colebrook-White correlation
        lambda = fsolve(z,0.0004);

end

L = 0.5;

g_M11 = (1 - M11^2)/(gamma11*M11^2) + ((gamma11 + 1)/(2*gamma11))*log(((gamma11 + 1)*M11^2)/(2 + (gamma11 - 1)*M11^2) );
g_M10 = g_M11 + lambda*(L/d_p_int);

y = @(x) g_M10 - (1 - x^2)/(gamma10*x^2) + ((gamma10 + 1)/(2*gamma10))*log(((gamma10 + 1)*x^2)/(2 + (gamma10 - 1)*x^2) );
M10 = fsolve(y,0.006);

T_star = T11/(0.5*(gamma11 + 1)/(1 + (gamma11 - 1)*0.5*M11^2));
T10 = T_star*(0.5*(gamma10 + 1)/(1 + (gamma10 - 1)*0.5*M10^2));

P_star = P11/((1/M11)*sqrt(0.5*(gamma11 + 1)/(1 + (gamma11 - 1)*0.5*M11^2)));
P10 = P_star*((1/M10)*sqrt(0.5*(gamma10 + 1)/(1 + (gamma10 - 1)*0.5*M10^2)));

iter = 0;
err = 1;

while err > 1e-3

    iter = iter + 1;

    gamma10_new = gamma_N2(find(T==round(T10)),find(abs(P - round(P10,1)) < 0.001));

    y_new = @(x) g_M10 - (1 - x^2)/(gamma10_new*x^2) + ((gamma10_new + 1)/(2*gamma10_new))*log(((gamma10_new + 1)*x^2)/(2 + (gamma10_new - 1)*x^2) );
    M10_new = fsolve(y,0.006);

    T10_new = T_star*(0.5*(gamma10_new + 1)/(1 + (gamma10_new - 1)*0.5*M10_new^2));

    P10_new = P_star*((1/M10)*sqrt(0.5*(gamma10 + 1)/(1 + (gamma10_new - 1)*0.5*M10_new^2)));  

    err = abs(gamma10 - gamma10_new);

    gamma10 = gamma10_new;
    M10 = M10_new;
    y = y_new;
    T10 = T10_new;
    P10 = P10_new;

end

clear gamma10_new
clear M10_new
clear y_new
clear T10_new
clear P10_new


    rho_star = rho11/((1/M11)*sqrt( 2*(1 + (gamma11 - 1)*0.5*M11^2)/(gamma11 + 1)));
    rho10 = rho_star*((1/M10)*sqrt( 2*(1 + (gamma10 - 1)*0.5*M10^2)/(gamma10 + 1)));

    c10 = sqrt(gamma10*R_N2*T10);
    v10 = c10*M10;