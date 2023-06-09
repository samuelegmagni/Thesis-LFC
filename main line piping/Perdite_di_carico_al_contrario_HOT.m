clear
clc

set(0,'DefaultTextFontSize',12);              % Settings for the plot
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultLegendInterpreter','Latex');
set(groot,'DefaultAxesTickLabelInterpreter','Latex');
set(0,'DefaultTextInterpreter','Latex');
set(0,'DefaultLegendFontSize',12);

T12 = 600;
P13 = 1.5;                           % Temperature and pressure at the exit of the injector

T = (floor(T12)-100):0.5:(ceil(T12));
P = (floor(P13)):0.1:(ceil(P13)+5);
data = nistdata('N2',T,P);

R_N2 = 8314/28;
R_SRP = 8314/36.32;
d12 = 16.05*1e-3;
A12 = 0.25*pi*d12^2;

m_dot_SRP = 18*1e-3;
m_dot_N2 = 60*1e-3;
m_dot12 = m_dot_N2 + m_dot_SRP;

rho_N2 = data.Rho*data.Mw;           % Density of Nitrogen [kg/m^3] 
cp_N2 = data.Cp/data.Mw;             % Specific heat at constant pressure of Nitrogen [J/kgK]
cv_N2 = data.Cv/data.Mw;             % Specific heat at constant volume of Nitrogen [J/kgK]
gamma_N2 = cp_N2./cv_N2;             % Ratio of specific heats [-]
mu_N2 = data.mu;  

gamma_N2_13 = gamma_N2(find(T==round(T12)),find(abs(P - round(P13,1)) < 0.001)); % Ratio of specific heats the check valve  [-]
gamma_N2_12 = gamma_N2(find(T==round(T12)),find(abs(P - round(P13,1)) < 0.001));

gamma_SRP = 1.1122;
gamma12 = (m_dot_SRP*gamma_SRP + m_dot_N2*gamma_N2_12)/(m_dot_N2 + m_dot_SRP);
gamma13 = (m_dot_SRP*gamma_SRP + m_dot_N2*gamma_N2_13)/(m_dot_N2 + m_dot_SRP);

%%

M13 = 1; 
% d_inj = 12*1e-3;
% A_inj = 0.25*pi*d_inj^2;

P_tot13 = P13*(1 + 0.5*(gamma13 - 1)*M13^2)^(gamma13/(gamma13 - 1));

iter = 0;
err = 1;

while err > 1e-3

    iter = iter + 1;

    M12 = 0.4;
    Area_ratio = (1/M12)*((2/(gamma12 + 1))*(1 + 0.5*(gamma12 - 1)*M12^2))^((gamma12 + 1)/(2*(gamma12 - 1)));
    A_inj = A12/Area_ratio;
    d_inj = sqrt((4*A_inj)/pi);
    % z = @(x) A12/A_inj - (M13/x)*((1 + 0.5*(gamma12 - 1)*x^2)/(1 + 0.5*(gamma13 - 1)*M13^2))^((gamma13 + 1)/(2*(gamma13 -1)));
    % M12 = fsolve(z,0.8);
    
    T_tot12 = T12*(1 + 0.5*(gamma12 - 1)*M12^2);
    P12 = P_tot13/(1 + 0.5*(gamma12 - 1)*M12^2)^(gamma12/(gamma12 - 1));
    T13 = T_tot12/(1 + 0.5*(gamma13 - 1)*M13^2);

    gamma_N2_12_new = gamma_N2(find(T==round(T12)),find(abs(P - round(P12,1)) < 0.001));
    gamma_N2_13 = gamma_N2(find(T==round(T13)),find(abs(P - round(P13,1)) < 0.001));


    gamma12 = (m_dot_SRP*gamma_SRP + m_dot_N2*gamma_N2_12_new)/(m_dot_N2 + m_dot_SRP);
    gamma13 = (m_dot_SRP*gamma_SRP + m_dot_N2*gamma_N2_13)/(m_dot_N2 + m_dot_SRP);

    err = abs(gamma_N2_12 - gamma_N2_12_new);

    gamma_N2_12 = gamma_N2_12_new;

end

clear gamma_N2_12_new

cp_N2_12 = cp_N2(find(T==round(T12)),find(abs(P - round(P12,1)) < 0.001));
rho_N2_12 = rho_N2(find(T==round(T12)),find(abs(P - round(P12,1)) < 0.001));
rho_N2_13 = rho_N2(find(T==round(T13)),find(abs(P - round(P13,1)) < 0.001));
cp_SRP = 5515.2;
rho_SRP = 0.92189;

cp_12 = (m_dot_SRP*cp_SRP + m_dot_N2*cp_N2_12)/(m_dot_N2 + m_dot_SRP);
rho12 = (m_dot_SRP*rho_SRP + m_dot_N2*rho_N2_12)/(m_dot_N2 + m_dot_SRP);
rho13 = (m_dot_SRP*rho_SRP + m_dot_N2*rho_N2_13)/(m_dot_N2 + m_dot_SRP);
R12 = (m_dot_SRP*R_SRP + m_dot_N2*R_N2)/(m_dot_N2 + m_dot_SRP);

v_inj = m_dot12/(A_inj*rho13)
c13 = sqrt(gamma13*R12*T13)

v12 = m_dot12/(A12*rho12)
P_tank = P12*1e5 + 0.5*rho12*v12^2; 

T_SRP = 1600;
T11 = (m_dot12*T12 - m_dot_SRP*T_SRP)/(m_dot_N2);

T = (floor(T11)-10):0.5:(ceil(T11));
P = (floor(P12)):0.1:(ceil(P12)+5);
data = nistdata('N2',T,P);

rho_N2 = data.Rho*data.Mw;           % Density of Nitrogen [kg/m^3] 
cp_N2 = data.Cp/data.Mw;             % Specific heat at constant pressure of Nitrogen [J/kgK]
cv_N2 = data.Cv/data.Mw;             % Specific heat at constant volume of Nitrogen [J/kgK]
gamma_N2 = cp_N2./cv_N2;             % Ratio of specific heats [-]
mu_N2 = data.mu;  

rho11 = rho_N2(find(T==round(T11)),find(abs(P - round(P12,1)) < 0.001));

iter = 0;
err = 1;

while err > 1e-3

    iter = iter + 1;

    v11 = m_dot_N2/(A12*rho11);
    P11 = (P_tank + rho11*v11^2)*1e-5;
   
    rho11_new = rho_N2(find(T==round(T11)),find(abs(P - round(P11,1)) < 0.001));

    err = abs(rho11 - rho11_new);

    rho11 = rho11_new;

end

clear rho11_new
