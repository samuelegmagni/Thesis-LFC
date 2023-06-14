%% Pressure drops in N2 line with 3/4 inch diameter tubes

clc
clear


set(0,'DefaultTextFontSize',12);              % Settings for the plot
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultLegendInterpreter','Latex');
set(groot,'DefaultAxesTickLabelInterpreter','Latex');
set(0,'DefaultTextInterpreter','Latex');
set(0,'DefaultLegendFontSize',12);

%% After pressure regulator (point 1)

T1 = 298;                                       % Temperature downstream the pressure regulator [K]
P1 = 30.2;  

T = (floor(T1)-3):0.5:(ceil(T1));
P = (floor(P1)-13):0.1:(ceil(P1));
data = nistdata('N2',T,P);

rho_N2 = data.Rho*data.Mw;           % Density of Nitrogen [kg/m^3] 
cp_N2 = data.Cp/data.Mw;             % Specific heat at constant pressure of Nitrogen [J/kgK]
cv_N2 = data.Cv/data.Mw;             % Specific heat at constant volume of Nitrogen [J/kgK]
gamma_N2 = cp_N2./cv_N2;             % Ratio of specific heats [-]
mu_N2 = data.mu;                     % Viscosity of Nitrogen [Pa*s]
m_dot_N2 = 65*1e-3;                  % Nitrogen mass flow rate [kg/s]
R = 8314/28;                         % Specific ideal gas constant [J/kgK]

d1_ext = 6.35*1e-3;                   % Pipe external diameter [m]
t1 = 1*1e-3;                        % Thickness of the tube  [m]
d1_int = d1_ext - 2*t1;             % Pipe internal diameter [m]
A1_int = pi*(d1_int/2)^2;            % Internal cross sectional area [m^2]

rho1 = rho_N2(find(T==round(T1)),find(abs(P - round(P1,1)) < 0.001)); % Density downstream the pressure regulator [kg/m^3]
v1 = m_dot_N2/(A1_int*rho1);                                             % Velocity downstream the pressure regulator [m/s]
P2 = 1e-5*(P1*1e5 - 1*rho1*v1^2);                                  % Pressure downstream the pipe bending after the pressure regulator [bar]
gamma1 = gamma_N2(find(T==round(T1)),find(abs(P - round(P1,1)) < 0.001));   % Ratio of specific heats downstream the pipe bending after the pressure regulator [-]
mu1 = mu_N2(find(T==round(T1)),find(abs(P - round(P1,1)) < 0.001));         % Viscosity downstream the pipe bending after the pressure regulator [Pa*s]
c1 = (gamma1*R*T1)^0.5;                         % Sound speed downstream the pipe bending after pressure regulator [m/s]
M1 = v1/c1;                                     % Mach number downstream the pipe bending after pressure regulator [-]
Re1 = (rho1*v1*d1_int)/mu1; 

%% 12 mm tube between blue and green fittings (point 2 and 3)
d2_3_ext = 12*1e-3;
t2_3 = 1.5*1e-3;
d2_3_int = d2_3_ext - 2*t2_3;
A2_3 = pi*(d2_3_int/2)^2;

eps = 0.015*1e-3;                    % Absolute roughness of stainless steel [m]
eps2_3_rel = eps/d2_3_int;               % Relative roughness of stainless steel [-]

T2 = T1;
rho2 = rho_N2(find(T==round(T2)),find(abs(P - round(P2,1)) < 0.001));       % Density downstream the pipe bending after the pressure regulator [kg/m^3]
gamma2 = gamma_N2(find(T==round(T2)),find(abs(P - round(P2,1)) < 0.001));   % Ratio of specific heats downstream the pipe bending after the pressure regulator [-]
mu2 = mu_N2(find(T==round(T2)),find(abs(P - round(P2,1)) < 0.001));         % Viscosity downstream the pipe bending after the pressure regulator [Pa*s]
v2 = m_dot_N2/(A2_3*rho2);                     % Gas velocity downstream the pipe bending after pressure regulator [m/s]
c2 = (gamma2*R*T2)^0.5;                         % Sound speed downstream the pipe bending after pressure regulator [m/s]
M2 = v2/c2;                                     % Mach number downstream the pipe bending after pressure regulator [-]
Re2 = (rho2*v2*d2_3_int)/mu2;

L2_3 = 5*1e-2;
gamma3 = gamma_N2(find(T==round(T2)),find(abs(P - round(P2,1)) < 0.001));   % Ratio of specific heats downstream the pipe bending after the pressure regulator [-]

if Re2 < 2300

        lambda = 64/Re2;

    else 

        z = @(x) 1/sqrt(x) + 2*log10(2.51/(Re2*sqrt(x)) + eps2_3_rel/3.71);   % Colebrook-White correlation
        lambda = fsolve(z,0.0004);

end

iter = 0;
err = 1;

while err > 1e-3

    iter = iter + 1;

    g_M2 = (1 - M2^2)/(gamma2*M2^2) + ((gamma2 + 1)/(2*gamma2))*log(((gamma2 + 1)*M2^2)/(2 + (gamma2 - 1)*M2^2) );
    g_M3 = g_M2 - lambda*(L2_3/d2_3_int);
    
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

%% After green fitting (12 mm -> 3/4") (point 4)

P4 = 1e-5*(P3*1e5 - 1*rho3*v3^2);

%% 3/4" mm tube before blue elbow (point 4 and 5)

d4_5_ext = 19.05*1e-3;
t4_5 = 1.5*1e-3;
d4_5_int = d4_5_ext - 2*t4_5;
A4_5 = pi*(d4_5_int/2)^2;

eps = 0.015*1e-3;                    % Absolute roughness of stainless steel [m]
eps4_5_rel = eps/d4_5_int;               % Relative roughness of stainless steel [-]

T4 = T3;
rho4 = rho_N2(find(T==round(T4)),find(abs(P - round(P4,1)) < 0.001));       % Density downstream the pipe bending after the pressure regulator [kg/m^3]
gamma4 = gamma_N2(find(T==round(T4)),find(abs(P - round(P4,1)) < 0.001));   % Ratio of specific heats downstream the pipe bending after the pressure regulator [-]
mu4 = mu_N2(find(T==round(T4)),find(abs(P - round(P4,1)) < 0.001));         % Viscosity downstream the pipe bending after the pressure regulator [Pa*s]
v4 = m_dot_N2/(A4_5*rho4);                     % Gas velocity downstream the pipe bending after pressure regulator [m/s]
c4 = (gamma4*R*T4)^0.5;                         % Sound speed downstream the pipe bending after pressure regulator [m/s]
M4 = v4/c4;                                     % Mach number downstream the pipe bending after pressure regulator [-]
Re4 = (rho4*v4*d4_5_int)/mu4;

L4_5 = 5*1e-2;
gamma5 = gamma_N2(find(T==round(T4)),find(abs(P - round(P4,1)) < 0.001));   % Ratio of specific heats downstream the pipe bending after the pressure regulator [-]

if Re4 < 2300

        lambda = 64/Re4;

    else 

        z = @(x) 1/sqrt(x) + 2*log10(2.51/(Re4*sqrt(x)) + eps4_5_rel/3.71);   % Colebrook-White correlation
        lambda = fsolve(z,0.0004);

end

iter = 0;
err = 1;

while err > 1e-3

    iter = iter + 1;

    g_M4 = (1 - M4^2)/(gamma4*M4^2) + ((gamma4 + 1)/(2*gamma4))*log(((gamma4 + 1)*M4^2)/(2 + (gamma4 - 1)*M4^2) );
    g_M5 = g_M4 - lambda*(L4_5/d4_5_int);
    
    y = @(x) g_M5 - (1 - x^2)/(gamma5*x^2) + ((gamma5 + 1)/(2*gamma5))*log(((gamma5 + 1)*x^2)/(2 + (gamma5 - 1)*x^2) );
    M5 = fsolve(y,0.006);
    
    T_star = T4/(0.5*(gamma4 + 1)/(1 + (gamma4 - 1)*0.5*M4^2));
    T5 = T_star*(0.5*(gamma5 + 1)/(1 + (gamma5 - 1)*0.5*M5^2));
    
    P_star = P4/((1/M4)*sqrt(0.5*(gamma4 + 1)/(1 + (gamma4 - 1)*0.5*M4^2)));
    P5 = P_star*((1/M5)*sqrt(0.5*(gamma5 + 1)/(1 + (gamma5 - 1)*0.5*M5^2)));
    
    rho_star = rho4/((1/M4)*sqrt( 2*(1 + (gamma4 - 1)*0.5*M4^2)/(gamma4 + 1)));
    rho5 = rho_star*((1/M5)*sqrt( 2*(1 + (gamma5 - 1)*0.5*M5^2)/(gamma5 + 1)));
    
    c5 = sqrt(gamma5*R*T5);
    v5 = M5*c5;

    gamma5_new = gamma_N2(find(T==round(T5)),find(abs(P - round(P5,1)) < 0.001));

    err = abs(gamma5 - gamma5_new);

    gamma5 = gamma5_new;

end

clear gamma5_new

%% After blue elbow (point 6)

P6 = 1e-5*(P5*1e5 - 0.7*rho5*v5^2);

%% 3/4" mm tube after blue elbow (point 6 and 7)

d6_7_ext = 19.05*1e-3;
t6_7 = 1.5*1e-3;
d6_7_int = d6_7_ext - 2*t6_7;
A6_7 = pi*(d6_7_int/2)^2;

eps = 0.015*1e-3;                    % Absolute roughness of stainless steel [m]
eps6_7_rel = eps/d6_7_int;               % Relative roughness of stainless steel [-]

T6 = T5;

T = (floor(T6)-3):0.5:(ceil(T6));
P = (floor(P6)-5):0.1:(ceil(P6));
data = nistdata('N2',T,P);

rho_N2 = data.Rho*data.Mw;           % Density of Nitrogen [kg/m^3] 
cp_N2 = data.Cp/data.Mw;             % Specific heat at constant pressure of Nitrogen [J/kgK]
cv_N2 = data.Cv/data.Mw;             % Specific heat at constant volume of Nitrogen [J/kgK]
gamma_N2 = cp_N2./cv_N2;             % Ratio of specific heats [-]
mu_N2 = data.mu; 

L6_7 = 5*1e-2;
gamma6 = gamma_N2(find(T==round(T6)),find(abs(P - round(P6,1)) < 0.001));   % Ratio of specific heats downstream the pipe bending after the pressure regulator [-]
gamma7 = gamma_N2(find(T==round(T6)),find(abs(P - round(P6,1)) < 0.001));   % Ratio of specific heats downstream the pipe bending after the pressure regulator [-]
mu6 = mu_N2(find(T==round(T6)),find(abs(P - round(P6,1)) < 0.001));   % Ratio of specific heats downstream the pipe bending after the pressure regulator [-]
rho6 = rho_N2(find(T==round(T6)),find(abs(P - round(P6,1)) < 0.001));   % Ratio of specific heats downstream the pipe bending after the pressure regulator [-]

v6 = m_dot_N2/(A6_7*rho6);                     % Gas velocity downstream the pipe bending after pressure regulator [m/s]
c6 = (gamma6*R*T6)^0.5;                         % Sound speed downstream the pipe bending after pressure regulator [m/s]
M6 = v6/c6;       

Re6 = (rho6*v6*d6_7_int)/mu6;

if Re6 < 2300

        lambda = 64/Re6;

    else 

        z = @(x) 1/sqrt(x) + 2*log10(2.51/(Re6*sqrt(x)) + eps6_7_rel/3.71);   % Colebrook-White correlation
        lambda = fsolve(z,0.0004);

end

iter = 0;
err = 1;

while err > 1e-3

    iter = iter + 1;

    g_M6 = (1 - M6^2)/(gamma6*M6^2) + ((gamma6 + 1)/(2*gamma6))*log(((gamma6 + 1)*M6^2)/(2 + (gamma6 - 1)*M6^2) );
    g_M7 = g_M6 - lambda*(L6_7/d6_7_int);
    
    y = @(x) g_M7 - (1 - x^2)/(gamma7*x^2) + ((gamma7 + 1)/(2*gamma7))*log(((gamma7 + 1)*x^2)/(2 + (gamma7 - 1)*x^2) );
    M7 = fsolve(y,0.006);
    
    T_star = T6/(0.5*(gamma6 + 1)/(1 + (gamma6 - 1)*0.5*M6^2));
    T7 = T_star*(0.5*(gamma7 + 1)/(1 + (gamma7 - 1)*0.5*M7^2));
    
    P_star = P6/((1/M6)*sqrt(0.5*(gamma6 + 1)/(1 + (gamma6 - 1)*0.5*M6^2)));
    P7 = P_star*((1/M7)*sqrt(0.5*(gamma7 + 1)/(1 + (gamma7 - 1)*0.5*M7^2)));
    
    rho_star = rho6/((1/M6)*sqrt( 2*(1 + (gamma6 - 1)*0.5*M6^2)/(gamma6 + 1)));
    rho7 = rho_star*((1/M7)*sqrt( 2*(1 + (gamma7 - 1)*0.5*M7^2)/(gamma7 + 1)));
    
    c7 = sqrt(gamma7*R*T7);
    v7 = M7*c7;

    gamma7_new = gamma_N2(find(T==round(T7)),find(abs(P - round(P7,1)) < 0.001));

    err = abs(gamma7 - gamma7_new);

    gamma7 = gamma7_new;

end

clear gamma7_new

%% After first T-fitting (point 8)

P8 = 1e-5*(P7*1e5 - 1.3*rho7*v7^2);

%% Between T-fitting and manual ball valve (point 8 and 9)

d8_9_ext = 19.05*1e-3;
t8_9 = 1.5*1e-3;
d8_9_int = d8_9_ext - 2*t8_9;
A8_9 = pi*(d8_9_int/2)^2;

eps = 0.015*1e-3;                    % Absolute roughness of stainless steel [m]
eps8_9_rel = eps/d8_9_int;               % Relative roughness of stainless steel [-]

T8 = T7;

T = (floor(T8)-3):0.5:(ceil(T8));
P = (floor(P8)-5):0.1:(ceil(P8));
data = nistdata('N2',T,P);

rho_N2 = data.Rho*data.Mw;           % Density of Nitrogen [kg/m^3] 
cp_N2 = data.Cp/data.Mw;             % Specific heat at constant pressure of Nitrogen [J/kgK]
cv_N2 = data.Cv/data.Mw;             % Specific heat at constant volume of Nitrogen [J/kgK]
gamma_N2 = cp_N2./cv_N2;             % Ratio of specific heats [-]
mu_N2 = data.mu; 

L8_9 = 5*1e-2;
gamma8 = gamma_N2(find(T==round(T8)),find(abs(P - round(P8,1)) < 0.001));   % Ratio of specific heats downstream the pipe bending after the pressure regulator [-]
gamma9 = gamma_N2(find(T==round(T8)),find(abs(P - round(P8,1)) < 0.001));   % Ratio of specific heats downstream the pipe bending after the pressure regulator [-]
mu8 = mu_N2(find(T==round(T8)),find(abs(P - round(P8,1)) < 0.001));   % Ratio of specific heats downstream the pipe bending after the pressure regulator [-]
rho8 = rho_N2(find(T==round(T8)),find(abs(P - round(P8,1)) < 0.001));   % Ratio of specific heats downstream the pipe bending after the pressure regulator [-]

v8 = m_dot_N2/(A8_9*rho8);                     % Gas velocity downstream the pipe bending after pressure regulator [m/s]
c8 = (gamma8*R*T8)^0.5;                         % Sound speed downstream the pipe bending after pressure regulator [m/s]
M8 = v8/c8;       

Re8 = (rho8*v8*d8_9_int)/mu8;

if Re8 < 2300

        lambda = 64/Re8;

    else 

        z = @(x) 1/sqrt(x) + 2*log10(2.51/(Re8*sqrt(x)) + eps8_9_rel/3.71);   % Colebrook-White correlation
        lambda = fsolve(z,0.0004);

end

iter = 0;
err = 1;

while err > 1e-3

    iter = iter + 1;

    g_M8 = (1 - M8^2)/(gamma8*M8^2) + ((gamma8 + 1)/(2*gamma8))*log(((gamma8 + 1)*M8^2)/(2 + (gamma8 - 1)*M8^2) );
    g_M9 = g_M8 - lambda*(L8_9/d8_9_int);
    
    y = @(x) g_M9 - (1 - x^2)/(gamma9*x^2) + ((gamma9 + 1)/(2*gamma9))*log(((gamma9 + 1)*x^2)/(2 + (gamma9 - 1)*x^2) );
    M9 = fsolve(y,0.006);
    
    T_star = T8/(0.5*(gamma8 + 1)/(1 + (gamma8 - 1)*0.5*M8^2));
    T9 = T_star*(0.5*(gamma9 + 1)/(1 + (gamma9 - 1)*0.5*M9^2));
    
    P_star = P8/((1/M8)*sqrt(0.5*(gamma8 + 1)/(1 + (gamma8 - 1)*0.5*M8^2)));
    P9 = P_star*((1/M9)*sqrt(0.5*(gamma9 + 1)/(1 + (gamma9 - 1)*0.5*M9^2)));
    
    rho_star = rho8/((1/M8)*sqrt( 2*(1 + (gamma8 - 1)*0.5*M8^2)/(gamma8 + 1)));
    rho9 = rho_star*((1/M9)*sqrt( 2*(1 + (gamma9 - 1)*0.5*M9^2)/(gamma9 + 1)));
    
    c9 = sqrt(gamma9*R*T9);
    v9 = M9*c9;

    gamma9_new = gamma_N2(find(T==round(T9)),find(abs(P - round(P9,1)) < 0.001));

    err = abs(gamma9 - gamma9_new);

    gamma9 = gamma9_new;

end

clear gamma9_new

%% After manual ball valve (point 10)

G_g = rho9/1000;                      % Nitrogen specific gravity [-]
q_N2 = (m_dot_N2/rho9)*1000;           % Nitrogen volumetric flow rate [L/s]
% q_N2_std = (P4*q_N2*T4*60)/(1*273.15); % Nitrogen volumetric flow rate at std conditions [std L/min]
C_V = 3.8;

P10 = P9 - (G_g*(q_N2*60)^2)/(14.42*C_V)^2;
T10 = T9;

%% After first manual ball valve (point 10) and before second T-fitting (point 11) 

d10_11_ext = 19.05*1e-3;
t10_11 = 1.5*1e-3;
d10_11_int = d10_11_ext - 2*t10_11;
A10_11 = pi*(d10_11_int/2)^2;

eps10_11_rel = eps/d10_11_int;               % Relative roughness of stainless steel [-]


T = (floor(T10)-5):0.5:(ceil(T10));
P = (floor(P10)-5):0.1:(ceil(P10));
data = nistdata('N2',T,P);

rho_N2 = data.Rho*data.Mw;           % Density of Nitrogen [kg/m^3] 
cp_N2 = data.Cp/data.Mw;             % Specific heat at constant pressure of Nitrogen [J/kgK]
cv_N2 = data.Cv/data.Mw;             % Specific heat at constant volume of Nitrogen [J/kgK]
gamma_N2 = cp_N2./cv_N2;             % Ratio of specific heats [-]
mu_N2 = data.mu;                     % Viscosity of Nitrogen [Pa*s]

L10_11 = 5*1e-2;
rho10 = rho_N2(find(T==round(T10)),find(abs(P - round(P10,1)) < 0.001));     % Density downstream the manual ball valve [kg/m^3]
gamma10 = gamma_N2(find(T==round(T10)),find(abs(P - round(P10,1)) < 0.001)); % Ratio of specific heats downstream the manual ball valve  [-]
gamma11 = gamma_N2(find(T==round(T10)),find(abs(P - round(P10,1)) < 0.001));
mu10 = mu_N2(find(T==round(T10)),find(abs(P - round(P10,1)) < 0.001));       % Viscosity downstream the manual ball valve [Pa*s]
v10 = m_dot_N2/(A10_11*rho10);                     % Gas velocity downstream the manual ball valve [m/s]
c10 = (gamma10*R*T10)^0.5;                         % Sound speed downstream the manual ball valve [m/s]
M10 = v10/c10;                                     % Mach number downstream the manual ball valve [-]
Re10 = (rho10*v10*d10_11_int)/mu10;                    % Reynolds number downstream the manual ball valve [-]

if Re10 < 2300

        lambda = 64/Re10;

    else 

        z = @(x) 1/sqrt(x) + 2*log10(2.51/(Re10*sqrt(x)) + eps10_11_rel/3.71);   % Colebrook-White correlation
        lambda = fsolve(z,0.0004);

end

iter = 0;
err = 1;

while err > 1e-3

    iter = iter + 1;

    g_M10 = (1 - M10^2)/(gamma10*M10^2) + ((gamma10 + 1)/(2*gamma10))*log(((gamma10 + 1)*M10^2)/(2 + (gamma10 - 1)*M10^2) );
    g_M11 = g_M10 - lambda*(L10_11/d10_11_int);

    y = @(x) g_M11 - (1 - x^2)/(gamma11*x^2) + ((gamma11 + 1)/(2*gamma11))*log(((gamma11 + 1)*x^2)/(2 + (gamma11 - 1)*x^2) );
    M11 = fsolve(y,0.006);
    
    T_star = T10/(0.5*(gamma10 + 1)/(1 + (gamma10 - 1)*0.5*M10^2));
    T11 = T_star*(0.5*(gamma11 + 1)/(1 + (gamma11 - 1)*0.5*M11^2));
    
    P_star = P10/((1/M10)*sqrt(0.5*(gamma10 + 1)/(1 + (gamma10 - 1)*0.5*M10^2)));
    P11 = P_star*((1/M11)*sqrt(0.5*(gamma11 + 1)/(1 + (gamma11 - 1)*0.5*M11^2)));

    rho_star = rho10/((1/M10)*sqrt( 2*(1 + (gamma10 - 1)*0.5*M10^2)/(gamma10 + 1)));
    rho11 = rho_star*((1/M11)*sqrt( 2*(1 + (gamma11 - 1)*0.5*M11^2)/(gamma11 + 1)));
    
    c11 = sqrt(gamma11*R*T11);
    v11 = c11*M11;

    gamma11_new = gamma_N2(find(T==round(T11)),find(abs(P - round(P11,1)) < 0.001));

    err = abs(gamma11 - gamma11_new);

    gamma11 = gamma11_new;

end

clear gamma11_new

%% After second T-fitting (point 12)

P12 = 1e-5*(P11*1e5 - 1.3*rho11*v11^2);
T12 = T11;

%% Between second and third T-fitting (point 12 and 13)

d12_13_ext = 19.05*1e-3;
t12_13 = 1.5*1e-3;
d12_13_int = d12_13_ext - 2*t12_13;
A12_13 = pi*(d12_13_int/2)^2;

eps12_13_rel = eps/d12_13_int;               % Relative roughness of stainless steel [-]


T = (floor(T12)-5):0.5:(ceil(T12));
P = (floor(P12)-5):0.1:(ceil(P12));
data = nistdata('N2',T,P);

rho_N2 = data.Rho*data.Mw;           % Density of Nitrogen [kg/m^3] 
cp_N2 = data.Cp/data.Mw;             % Specific heat at constant pressure of Nitrogen [J/kgK]
cv_N2 = data.Cv/data.Mw;             % Specific heat at constant volume of Nitrogen [J/kgK]
gamma_N2 = cp_N2./cv_N2;             % Ratio of specific heats [-]
mu_N2 = data.mu;                     % Viscosity of Nitrogen [Pa*s]

L12_13 = 5*1e-2;
rho12 = rho_N2(find(T==round(T12)),find(abs(P - round(P12,1)) < 0.001));     % Density downstream the manual ball valve [kg/m^3]
gamma12 = gamma_N2(find(T==round(T12)),find(abs(P - round(P12,1)) < 0.001)); % Ratio of specific heats downstream the manual ball valve  [-]
gamma13 = gamma_N2(find(T==round(T12)),find(abs(P - round(P12,1)) < 0.001));
mu12 = mu_N2(find(T==round(T12)),find(abs(P - round(P12,1)) < 0.001));       % Viscosity downstream the manual ball valve [Pa*s]
v12 = m_dot_N2/(A12_13*rho12);                     % Gas velocity downstream the manual ball valve [m/s]
c12 = (gamma12*R*T12)^0.5;                         % Sound speed downstream the manual ball valve [m/s]
M12 = v12/c12;                                     % Mach number downstream the manual ball valve [-]
Re12 = (rho12*v12*d12_13_int)/mu12;                    % Reynolds number downstream the manual ball valve [-]

if Re12 < 2300

        lambda = 64/Re12;

    else 

        z = @(x) 1/sqrt(x) + 2*log10(2.51/(Re12*sqrt(x)) + eps12_13_rel/3.71);   % Colebrook-White correlation
        lambda = fsolve(z,0.0004);

end

iter = 0;
err = 1;

while err > 1e-3

    iter = iter + 1;

    g_M12 = (1 - M12^2)/(gamma12*M12^2) + ((gamma12 + 1)/(2*gamma12))*log(((gamma12 + 1)*M12^2)/(2 + (gamma12 - 1)*M12^2) );
    g_M13 = g_M12 - lambda*(L12_13/d12_13_int);

    y = @(x) g_M13 - (1 - x^2)/(gamma13*x^2) + ((gamma13 + 1)/(2*gamma13))*log(((gamma13 + 1)*x^2)/(2 + (gamma13 - 1)*x^2) );
    M13 = fsolve(y,0.006);
    
    T_star = T12/(0.5*(gamma12 + 1)/(1 + (gamma12 - 1)*0.5*M12^2));
    T13 = T_star*(0.5*(gamma13 + 1)/(1 + (gamma13 - 1)*0.5*M13^2));
    
    P_star = P12/((1/M12)*sqrt(0.5*(gamma12 + 1)/(1 + (gamma12 - 1)*0.5*M12^2)));
    P13 = P_star*((1/M13)*sqrt(0.5*(gamma13 + 1)/(1 + (gamma13 - 1)*0.5*M13^2)));

    rho_star = rho12/((1/M12)*sqrt( 2*(1 + (gamma12 - 1)*0.5*M12^2)/(gamma12 + 1)));
    rho13 = rho_star*((1/M13)*sqrt( 2*(1 + (gamma13 - 1)*0.5*M13^2)/(gamma13 + 1)));
    
    c13 = sqrt(gamma13*R*T13);
    v13 = c13*M13;

    gamma13_new = gamma_N2(find(T==round(T13)),find(abs(P - round(P13,1)) < 0.001));

    err = abs(gamma13 - gamma13_new);

    gamma13 = gamma13_new;

end

clear gamma13_new

%% After third T-fitting (point 14)

P14 = 1e-5*(P13*1e5 - 1.3*rho13*v13^2);
T14 = T13;


%% After third T-fitting and before green fitting (point 14 and 15)

d14_15_ext = 19.05*1e-3;
t14_15 = 1.5*1e-3;
d14_15_int = d14_15_ext - 2*t14_15;
A14_15 = pi*(d14_15_int/2)^2;

eps14_15_rel = eps/d14_15_int;               % Relative roughness of stainless steel [-]


T = (floor(T14)-5):0.5:(ceil(T14));
P = (floor(P14)-5):0.1:(ceil(P14));
data = nistdata('N2',T,P);

rho_N2 = data.Rho*data.Mw;           % Density of Nitrogen [kg/m^3] 
cp_N2 = data.Cp/data.Mw;             % Specific heat at constant pressure of Nitrogen [J/kgK]
cv_N2 = data.Cv/data.Mw;             % Specific heat at constant volume of Nitrogen [J/kgK]
gamma_N2 = cp_N2./cv_N2;             % Ratio of specific heats [-]
mu_N2 = data.mu;                     % Viscosity of Nitrogen [Pa*s]

L14_15 = 5*1e-2;
rho14 = rho_N2(find(T==round(T14)),find(abs(P - round(P14,1)) < 0.001));     % Density downstream the manual ball valve [kg/m^3]
gamma14 = gamma_N2(find(T==round(T14)),find(abs(P - round(P14,1)) < 0.001)); % Ratio of specific heats downstream the manual ball valve  [-]
gamma15 = gamma_N2(find(T==round(T14)),find(abs(P - round(P14,1)) < 0.001));
mu14 = mu_N2(find(T==round(T14)),find(abs(P - round(P14,1)) < 0.001));       % Viscosity downstream the manual ball valve [Pa*s]
v14 = m_dot_N2/(A14_15*rho14);                     % Gas velocity downstream the manual ball valve [m/s]
c14 = (gamma12*R*T14)^0.5;                         % Sound speed downstream the manual ball valve [m/s]
M14 = v14/c14;                                     % Mach number downstream the manual ball valve [-]
Re14 = (rho14*v14*d14_15_int)/mu14;                    % Reynolds number downstream the manual ball valve [-]

if Re14 < 2300

        lambda = 64/Re14;

    else 

        z = @(x) 1/sqrt(x) + 2*log10(2.51/(Re14*sqrt(x)) + eps14_15_rel/3.71);   % Colebrook-White correlation
        lambda = fsolve(z,0.0004);

end

iter = 0;
err = 1;

while err > 1e-3

    iter = iter + 1;

    g_M14 = (1 - M14^2)/(gamma14*M14^2) + ((gamma14 + 1)/(2*gamma14))*log(((gamma14 + 1)*M14^2)/(2 + (gamma14 - 1)*M14^2) );
    g_M15 = g_M14 - lambda*(L14_15/d14_15_int);

    y = @(x) g_M15 - (1 - x^2)/(gamma15*x^2) + ((gamma15 + 1)/(2*gamma15))*log(((gamma15 + 1)*x^2)/(2 + (gamma15 - 1)*x^2) );
    M15 = fsolve(y,0.006);
    
    T_star = T14/(0.5*(gamma14 + 1)/(1 + (gamma14 - 1)*0.5*M14^2));
    T15 = T_star*(0.5*(gamma15 + 1)/(1 + (gamma15 - 1)*0.5*M15^2));
    
    P_star = P14/((1/M14)*sqrt(0.5*(gamma14 + 1)/(1 + (gamma14 - 1)*0.5*M14^2)));
    P15 = P_star*((1/M15)*sqrt(0.5*(gamma15 + 1)/(1 + (gamma15 - 1)*0.5*M15^2)));

    rho_star = rho14/((1/M14)*sqrt( 2*(1 + (gamma14 - 1)*0.5*M14^2)/(gamma14 + 1)));
    rho15 = rho_star*((1/M15)*sqrt( 2*(1 + (gamma15 - 1)*0.5*M15^2)/(gamma15 + 1)));
    
    c15 = sqrt(gamma15*R*T15);
    v15 = c15*M15;

    gamma15_new = gamma_N2(find(T==round(T15)),find(abs(P - round(P15,1)) < 0.001));

    err = abs(gamma15 - gamma15_new);

    gamma15 = gamma15_new;

end

clear gamma15_new


%% After green fitting (3/4" -> 1/2") (point 16)

P16 = 1e-5*(P15*1e5 - 0.5*rho15*v15^2);
T16 = T15;


%% After green T-fitting and before pink fitting (point 16 and 17)

d16_17_ext = 12.07*1e-3;
t16_17 = 1.5*1e-3;
d16_17_int = d16_17_ext - 2*t16_17;
A16_17 = pi*(d16_17_int/2)^2;

eps16_17_rel = eps/d16_17_int;               % Relative roughness of stainless steel [-]


T = (floor(T16)-5):0.5:(ceil(T16));
P = (floor(P16)-5):0.1:(ceil(P16));
data = nistdata('N2',T,P);

rho_N2 = data.Rho*data.Mw;           % Density of Nitrogen [kg/m^3] 
cp_N2 = data.Cp/data.Mw;             % Specific heat at constant pressure of Nitrogen [J/kgK]
cv_N2 = data.Cv/data.Mw;             % Specific heat at constant volume of Nitrogen [J/kgK]
gamma_N2 = cp_N2./cv_N2;             % Ratio of specific heats [-]
mu_N2 = data.mu;                     % Viscosity of Nitrogen [Pa*s]

L16_17 = 5*1e-2;
rho16 = rho_N2(find(T==round(T16)),find(abs(P - round(P16,1)) < 0.001));     % Density downstream the manual ball valve [kg/m^3]
gamma16 = gamma_N2(find(T==round(T16)),find(abs(P - round(P16,1)) < 0.001)); % Ratio of specific heats downstream the manual ball valve  [-]
gamma17 = gamma_N2(find(T==round(T16)),find(abs(P - round(P16,1)) < 0.001));
mu16 = mu_N2(find(T==round(T16)),find(abs(P - round(P16,1)) < 0.001));       % Viscosity downstream the manual ball valve [Pa*s]
v16 = m_dot_N2/(A16_17*rho16);                     % Gas velocity downstream the manual ball valve [m/s]
c16 = (gamma16*R*T16)^0.5;                         % Sound speed downstream the manual ball valve [m/s]
M16 = v16/c16;                                     % Mach number downstream the manual ball valve [-]
Re16 = (rho16*v16*d16_17_int)/mu16;                    % Reynolds number downstream the manual ball valve [-]

if Re16 < 2300

        lambda = 64/Re16;

    else 

        z = @(x) 1/sqrt(x) + 2*log10(2.51/(Re16*sqrt(x)) + eps16_17_rel/3.71);   % Colebrook-White correlation
        lambda = fsolve(z,0.0004);

end

iter = 0;
err = 1;

while err > 1e-3

    iter = iter + 1;

    g_M16 = (1 - M16^2)/(gamma16*M16^2) + ((gamma16 + 1)/(2*gamma16))*log(((gamma16 + 1)*M16^2)/(2 + (gamma16 - 1)*M16^2) );
    g_M17 = g_M16 - lambda*(L16_17/d16_17_int);

    y = @(x) g_M16 - (1 - x^2)/(gamma16*x^2) + ((gamma16 + 1)/(2*gamma16))*log(((gamma16 + 1)*x^2)/(2 + (gamma16 - 1)*x^2) );
    M17 = fsolve(y,0.006);
    
    T_star = T16/(0.5*(gamma16 + 1)/(1 + (gamma16 - 1)*0.5*M16^2));
    T17 = T_star*(0.5*(gamma17 + 1)/(1 + (gamma17 - 1)*0.5*M17^2));
    
    P_star = P16/((1/M16)*sqrt(0.5*(gamma16 + 1)/(1 + (gamma16 - 1)*0.5*M16^2)));
    P17 = P_star*((1/M17)*sqrt(0.5*(gamma17 + 1)/(1 + (gamma17 - 1)*0.5*M17^2)));

    rho_star = rho16/((1/M16)*sqrt( 2*(1 + (gamma16 - 1)*0.5*M16^2)/(gamma16 + 1)));
    rho17 = rho_star*((1/M17)*sqrt( 2*(1 + (gamma17 - 1)*0.5*M17^2)/(gamma17 + 1)));
    
    c17 = sqrt(gamma17*R*T17);
    v17 = c17*M17;

    gamma17_new = gamma_N2(find(T==round(T17)),find(abs(P - round(P17,1)) < 0.001));

    err = abs(gamma17 - gamma17_new);

    gamma17 = gamma17_new;

end

clear gamma17_new


%% After pink fitting ( 1/2" -> 12mm) (point 18)

P18 = 1e-5*(P17*1e5 - 0.5*rho17*v17^2);
T18 = T17;


rho18 = rho_N2(find(T==round(T18)),find(abs(P - round(P18,1)) < 0.001));       % Density downstream the pipe bending after the pressure regulator [kg/m^3]
gamma18 = gamma_N2(find(T==round(T18)),find(abs(P - round(P18,1)) < 0.001));   % Ratio of specific heats downstream the pipe bending after the pressure regulator [-]
mu18 = mu_N2(find(T==round(T18)),find(abs(P - round(P18,1)) < 0.001));         % Viscosity downstream the pipe bending after the pressure regulator [Pa*s]
v18 = m_dot_N2/(A16_17*rho18);                     % Gas velocity downstream the pipe bending after pressure regulator [m/s]
c18 = (gamma18*R*T18)^0.5;                         % Sound speed downstream the pipe bending after pressure regulator [m/s]
M18 = v18/c18;    

%% Before MFM (point 18) and after MFM (point 19)
d18_19_ext = 12*1e-3;
t18_19 = 1.5*1e-3;
d18_19_int = d18_19_ext - 2*t18_19;
A18_19 = pi*(d18_19_int/2)^2;

G_g = rho18/1000;                       % Nitrogen specific gravity [-]
q_N2 = (m_dot_N2/rho18)*1000;           % Nitrogen volumetric flow rate [L/s]
C_V = 0.73;                             % Flow coefficient needle valve

P19 = P18 - (G_g*(q_N2*60)^2)/(14.42*C_V)^2;                                            % Pressure downstream the mass flow meter (needle valve approx) [bar]
T19 = T18;     

T = (floor(T19)-5):0.5:(ceil(T19));
P = (floor(P19)-4):0.1:(ceil(P19));
data = nistdata('N2',T,P);

rho_N2 = data.Rho*data.Mw;           % Density of Nitrogen [kg/m^3] 
cp_N2 = data.Cp/data.Mw;             % Specific heat at constant pressure of Nitrogen [J/kgK]
cv_N2 = data.Cv/data.Mw;             % Specific heat at constant volume of Nitrogen [J/kgK]
gamma_N2 = cp_N2./cv_N2;             % Ratio of specific heats [-]
mu_N2 = data.mu;                     % Viscosity of Nitrogen [Pa*s]


rho19 = rho_N2(find(T==round(T19)),find(abs(P - round(P19,1)) < 0.001));       % Density downstream the pipe bending after the pressure regulator [kg/m^3]
gamma19 = gamma_N2(find(T==round(T19)),find(abs(P - round(P19,1)) < 0.001));   % Ratio of specific heats downstream the pipe bending after the pressure regulator [-]
mu19 = mu_N2(find(T==round(T19)),find(abs(P - round(P19,1)) < 0.001));         % Viscosity downstream the pipe bending after the pressure regulator [Pa*s]
v19 = m_dot_N2/(A18_19*rho19);                     % Gas velocity downstream the pipe bending after pressure regulator [m/s]
c19 = (gamma19*R*T19)^0.5;                         % Sound speed downstream the pipe bending after pressure regulator [m/s]
M19 = v19/c19;   

%% After pink fitting ( 1/2" -> 12mm) (point 20)

P20 = 1e-5*(P19*1e5 - 1*rho19*v19^2);
T20 = T19;


%% After pink fitting and before green fitting (point 20 and 21)

d20_21_ext = 12.07*1e-3;
t20_21 = 1.5*1e-3;
d20_21_int = d20_21_ext - 2*t20_21;
A20_21 = pi*(d20_21_int/2)^2;

eps20_21_rel = eps/d20_21_int;               % Relative roughness of stainless steel [-]


T = (floor(T20)-5):0.5:(ceil(T20));
P = (floor(P20)-4):0.1:(ceil(P20));
data = nistdata('N2',T,P);

rho_N2 = data.Rho*data.Mw;           % Density of Nitrogen [kg/m^3] 
cp_N2 = data.Cp/data.Mw;             % Specific heat at constant pressure of Nitrogen [J/kgK]
cv_N2 = data.Cv/data.Mw;             % Specific heat at constant volume of Nitrogen [J/kgK]
gamma_N2 = cp_N2./cv_N2;             % Ratio of specific heats [-]
mu_N2 = data.mu;                     % Viscosity of Nitrogen [Pa*s]

L20_21 = 5*1e-2;
rho20 = rho_N2(find(T==round(T20)),find(abs(P - round(P20,1)) < 0.001));       % Density downstream the pipe bending after the pressure regulator [kg/m^3]
gamma20 = gamma_N2(find(T==round(T20)),find(abs(P - round(P20,1)) < 0.001));   % Ratio of specific heats downstream the pipe bending after the pressure regulator [-]
gamma21 = gamma_N2(find(T==round(T20)),find(abs(P - round(P20,1)) < 0.001));
mu20 = mu_N2(find(T==round(T20)),find(abs(P - round(P20,1)) < 0.001));         % Viscosity downstream the pipe bending after the pressure regulator [Pa*s]
v20 = m_dot_N2/(A20_21*rho20);                     % Gas velocity downstream the pipe bending after pressure regulator [m/s]
c20 = (gamma20*R*T20)^0.5;                         % Sound speed downstream the pipe bending after pressure regulator [m/s]
M20 = v20/c20;           
                  
Re20 = (rho20*v20*d20_21_int)/mu20;                    % Reynolds number downstream the manual ball valve [-]

if Re20 < 2300

        lambda = 64/Re20;

    else 

        z = @(x) 1/sqrt(x) + 2*log10(2.51/(Re20*sqrt(x)) + eps20_21_rel/3.71);   % Colebrook-White correlation
        lambda = fsolve(z,0.0004);

end

iter = 0;
err = 1;

while err > 1e-3

    iter = iter + 1;

    g_M20 = (1 - M20^2)/(gamma20*M20^2) + ((gamma20 + 1)/(2*gamma20))*log(((gamma20 + 1)*M20^2)/(2 + (gamma20 - 1)*M20^2) );
    g_M21 = g_M20- lambda*(L20_21/d20_21_int);

    y = @(x) g_M21 - (1 - x^2)/(gamma21*x^2) + ((gamma21 + 1)/(2*gamma21))*log(((gamma21 + 1)*x^2)/(2 + (gamma21 - 1)*x^2) );
    M21 = fsolve(y,0.006);
    
    T_star = T20/(0.5*(gamma20 + 1)/(1 + (gamma20 - 1)*0.5*M20^2));
    T21= T_star*(0.5*(gamma21+ 1)/(1 + (gamma21 - 1)*0.5*M21^2));
    
    P_star = P20/((1/M20)*sqrt(0.5*(gamma20 + 1)/(1 + (gamma20 - 1)*0.5*M20^2)));
    P21 = P_star*((1/M21)*sqrt(0.5*(gamma21 + 1)/(1 + (gamma21 - 1)*0.5*M21^2)));

    rho_star = rho20/((1/M20)*sqrt( 2*(1 + (gamma20 - 1)*0.5*M20^2)/(gamma20 + 1)));
    rho21 = rho_star*((1/M21)*sqrt( 2*(1 + (gamma21 - 1)*0.5*M21^2)/(gamma21 + 1)));
    
    c21 = sqrt(gamma21*R*T21);
    v21 = c21*M21;

    gamma21_new = gamma_N2(find(T==round(T21)),find(abs(P - round(P21,1)) < 0.001));

    err = abs(gamma21- gamma21_new);

    gamma21 = gamma21_new;

end

clear gamma21_new


%% After green fitting ( 1/2"-> 3/4 ") (point 22)

P22 = 1e-5*(P21*1e5 - 1*rho21*v21^2);
T22 = T21;


%% After green fitting and before prenumatic valve (point 22 and 23)

d22_23_ext = 19.05*1e-3;
t22_23 = 1.5*1e-3;
d22_23_int = d22_23_ext - 2*t22_23;
A22_23 = pi*(d22_23_int/2)^2;

eps22_23_rel = eps/d22_23_int;               % Relative roughness of stainless steel [-]


T = (floor(T22)-5):0.5:(ceil(T22));
P = (floor(P22)-3):0.1:(ceil(P22));
data = nistdata('N2',T,P);

rho_N2 = data.Rho*data.Mw;           % Density of Nitrogen [kg/m^3] 
cp_N2 = data.Cp/data.Mw;             % Specific heat at constant pressure of Nitrogen [J/kgK]
cv_N2 = data.Cv/data.Mw;             % Specific heat at constant volume of Nitrogen [J/kgK]
gamma_N2 = cp_N2./cv_N2;             % Ratio of specific heats [-]
mu_N2 = data.mu;                     % Viscosity of Nitrogen [Pa*s]

L22_23 = 30*1e-2;
rho22 = rho_N2(find(T==round(T22)),find(abs(P - round(P22,1)) < 0.001));       % Density downstream the pipe bending after the pressure regulator [kg/m^3]
gamma22 = gamma_N2(find(T==round(T22)),find(abs(P - round(P22,1)) < 0.001));   % Ratio of specific heats downstream the pipe bending after the pressure regulator [-]
gamma23 = gamma_N2(find(T==round(T22)),find(abs(P - round(P22,1)) < 0.001));
mu22 = mu_N2(find(T==round(T22)),find(abs(P - round(P22,1)) < 0.001));         % Viscosity downstream the pipe bending after the pressure regulator [Pa*s]
v22 = m_dot_N2/(A22_23*rho22);                     % Gas velocity downstream the pipe bending after pressure regulator [m/s]
c22 = (gamma22*R*T22)^0.5;                         % Sound speed downstream the pipe bending after pressure regulator [m/s]
M22 = v22/c22;           
                  
Re22 = (rho22*v22*d22_23_int)/mu22;                    % Reynolds number downstream the manual ball valve [-]

if Re22< 2300

        lambda = 64/Re22;

    else 

        z = @(x) 1/sqrt(x) + 2*log10(2.51/(Re22*sqrt(x)) + eps22_23_rel/3.71);   % Colebrook-White correlation
        lambda = fsolve(z,0.0004);

end

iter = 0;
err = 1;

while err > 1e-3

    iter = iter + 1;

    g_M22 = (1 - M22^2)/(gamma22*M22^2) + ((gamma22 + 1)/(2*gamma22))*log(((gamma22 + 1)*M22^2)/(2 + (gamma22 - 1)*M22^2) );
    g_M23 = g_M22- lambda*(L22_23/d22_23_int);

    y = @(x) g_M23 - (1 - x^2)/(gamma23*x^2) + ((gamma23 + 1)/(2*gamma23))*log(((gamma23 + 1)*x^2)/(2 + (gamma23 - 1)*x^2) );
    M23 = fsolve(y,0.006);
    
    T_star = T22/(0.5*(gamma22 + 1)/(1 + (gamma22 - 1)*0.5*M22^2));
    T23= T_star*(0.5*(gamma23+ 1)/(1 + (gamma23 - 1)*0.5*M23^2));
    
    P_star = P22/((1/M22)*sqrt(0.5*(gamma22 + 1)/(1 + (gamma22 - 1)*0.5*M22^2)));
    P23 = P_star*((1/M23)*sqrt(0.5*(gamma23 + 1)/(1 + (gamma23 - 1)*0.5*M23^2)));

    rho_star = rho22/((1/M22)*sqrt( 2*(1 + (gamma22 - 1)*0.5*M22^2)/(gamma22 + 1)));
    rho23 = rho_star*((1/M23)*sqrt( 2*(1 + (gamma23 - 1)*0.5*M23^2)/(gamma23 + 1)));
    
    c23 = sqrt(gamma23*R*T23);
    v23 = c23*M23;

    gamma23_new = gamma_N2(find(T==round(T23)),find(abs(P - round(P23,1)) < 0.001));

    err = abs(gamma23- gamma23_new);

    gamma23 = gamma23_new;

end

clear gamma23_new


%% Before servovalve (point 23) and after servovalve (point 24)
G_g = rho23/1000;                      % Nitrogen specific gravity [-]
q_N2 = (m_dot_N2/rho23)*1000;           % Nitrogen volumetric flow rate [L/s]
C_V = 3.8;                             % Flow coefficient ball valve

P24 = P23 - (G_g*(q_N2*60)^2)/(14.42*C_V)^2;     % Pressure downstream the servovalve (ball valve approx) [bar]
T24 = T23;                        % Temperature downstream the servovalve (ball valve approx) [K]

%% After servovalve and before check valve (point 24 and 25)

d24_25_ext = 19.05*1e-3;
t24_25 = 1.5*1e-3;
d24_25_int = d24_25_ext - 2*t24_25;
A24_25 = pi*(d24_25_int/2)^2;

eps24_25_rel = eps/d24_25_int;               % Relative roughness of stainless steel [-]


T = (floor(T24)-5):0.5:(ceil(T24));
P = (floor(P24)-3):0.1:(ceil(P24));
data = nistdata('N2',T,P);

rho_N2 = data.Rho*data.Mw;           % Density of Nitrogen [kg/m^3] 
cp_N2 = data.Cp/data.Mw;             % Specific heat at constant pressure of Nitrogen [J/kgK]
cv_N2 = data.Cv/data.Mw;             % Specific heat at constant volume of Nitrogen [J/kgK]
gamma_N2 = cp_N2./cv_N2;             % Ratio of specific heats [-]
mu_N2 = data.mu;                     % Viscosity of Nitrogen [Pa*s]

L24_25 = 30*1e-2;
rho24 = rho_N2(find(T==round(T24)),find(abs(P - round(P24,1)) < 0.001));       % Density downstream the pipe bending after the pressure regulator [kg/m^3]
gamma24 = gamma_N2(find(T==round(T24)),find(abs(P - round(P24,1)) < 0.001));   % Ratio of specific heats downstream the pipe bending after the pressure regulator [-]
gamma25 = gamma_N2(find(T==round(T24)),find(abs(P - round(P24,1)) < 0.001));
mu24 = mu_N2(find(T==round(T24)),find(abs(P - round(P24,1)) < 0.001));         % Viscosity downstream the pipe bending after the pressure regulator [Pa*s]
v24 = m_dot_N2/(A24_25*rho24);                     % Gas velocity downstream the pipe bending after pressure regulator [m/s]
c24 = (gamma24*R*T24)^0.5;                         % Sound speed downstream the pipe bending after pressure regulator [m/s]
M24 = v24/c24;           
                  
Re24 = (rho24*v24*d24_25_int)/mu24;                    % Reynolds number downstream the manual ball valve [-]

if Re24< 2300

        lambda = 64/Re24;

    else 

        z = @(x) 1/sqrt(x) + 2*log10(2.51/(Re24*sqrt(x)) + eps24_25_rel/3.71);   % Colebrook-White correlation
        lambda = fsolve(z,0.0004);

end

iter = 0;
err = 1;

while err > 1e-3

    iter = iter + 1;

    g_M24 = (1 - M24^2)/(gamma24*M24^2) + ((gamma24 + 1)/(2*gamma24))*log(((gamma24 + 1)*M24^2)/(2 + (gamma24 - 1)*M24^2) );
    g_M25 = g_M24- lambda*(L24_25/d24_25_int);

    y = @(x) g_M25 - (1 - x^2)/(gamma25*x^2) + ((gamma25 + 1)/(2*gamma25))*log(((gamma25+ 1)*x^2)/(2 + (gamma25 - 1)*x^2) );
    M25 = fsolve(y,0.006);
    
    T_star = T24/(0.5*(gamma24 + 1)/(1 + (gamma24 - 1)*0.5*M24^2));
    T25= T_star*(0.5*(gamma25+ 1)/(1 + (gamma25 - 1)*0.5*M25^2));
    
    P_star = P24/((1/M24)*sqrt(0.5*(gamma24 + 1)/(1 + (gamma24 - 1)*0.5*M24^2)));
    P25 = P_star*((1/M25)*sqrt(0.5*(gamma25 + 1)/(1 + (gamma25 - 1)*0.5*M25^2)));

    rho_star = rho24/((1/M24)*sqrt( 2*(1 + (gamma24 - 1)*0.5*M24^2)/(gamma24 + 1)));
    rho25 = rho_star*((1/M25)*sqrt( 2*(1 + (gamma25 - 1)*0.5*M25^2)/(gamma25 + 1)));
    
    c25 = sqrt(gamma25*R*T25);
    v25 = c25*M25;

    gamma25_new = gamma_N2(find(T==round(T25)),find(abs(P - round(P25,1)) < 0.001));

    err = abs(gamma25- gamma25_new);

    gamma25 = gamma25_new;

end

clear gamma25_new

%% Before check valve (point 25) and after check valve (point 26)
G_g = rho25/1000;                      % Nitrogen specific gravity [-]
q_N2 = (m_dot_N2/rho25)*1000;           % Nitrogen volumetric flow rate [L/s]

C_V = 1.68;                             % Flow coefficient check valve
P26 = P25 - (G_g*(q_N2*60)^2)/(14.42*C_V)^2;         % Pressure downstream the check valve [bar]
T26 = T25;


%% After checkvalve and before mixing chamber (point 26 and 27)

d26_27_ext = 19.05*1e-3;
t26_27 = 1.5*1e-3;
d26_27_int = d26_27_ext - 2*t24_25;
A26_27 = pi*(d26_27_int/2)^2;

eps26_27_rel = eps/d26_27_int;               % Relative roughness of stainless steel [-]


T = (floor(T26)-10):0.5:(ceil(T26));
P = (floor(P26)-2):0.1:(ceil(P26));
data = nistdata('N2',T,P);

rho_N2 = data.Rho*data.Mw;           % Density of Nitrogen [kg/m^3] 
cp_N2 = data.Cp/data.Mw;             % Specific heat at constant pressure of Nitrogen [J/kgK]
cv_N2 = data.Cv/data.Mw;             % Specific heat at constant volume of Nitrogen [J/kgK]
gamma_N2 = cp_N2./cv_N2;             % Ratio of specific heats [-]
mu_N2 = data.mu;                     % Viscosity of Nitrogen [Pa*s]

L26_27 = 30*1e-2;
rho26 = rho_N2(find(T==round(T26)),find(abs(P - round(P26,1)) < 0.001));       % Density downstream the pipe bending after the pressure regulator [kg/m^3]
gamma26 = gamma_N2(find(T==round(T26)),find(abs(P - round(P26,1)) < 0.001));   % Ratio of specific heats downstream the pipe bending after the pressure regulator [-]
gamma27 =gamma_N2(find(T==round(T26)),find(abs(P - round(P26,1)) < 0.001));
mu26 = mu_N2(find(T==round(T26)),find(abs(P - round(P26,1)) < 0.001));         % Viscosity downstream the pipe bending after the pressure regulator [Pa*s]
v26 = m_dot_N2/(A26_27*rho26);                     % Gas velocity downstream the pipe bending after pressure regulator [m/s]
c26 = (gamma26*R*T26)^0.5;                         % Sound speed downstream the pipe bending after pressure regulator [m/s]
M26 = v26/c26;           
                  
Re26 = (rho26*v26*d26_27_int)/mu26;                    % Reynolds number downstream the manual ball valve [-]

if Re26< 2300

        lambda = 64/Re26;

    else 

        z = @(x) 1/sqrt(x) + 2*log10(2.51/(Re26*sqrt(x)) + eps26_27_rel/3.71);   % Colebrook-White correlation
        lambda = fsolve(z,0.0004);

end

iter = 0;
err = 1;

while err > 1e-3

    iter = iter + 1;

    g_M26 = (1 - M26^2)/(gamma26*M26^2) + ((gamma26 + 1)/(2*gamma26))*log(((gamma26 + 1)*M26^2)/(2 + (gamma26 - 1)*M26^2) );
    g_M27 = g_M26- lambda*(L26_27/d26_27_int);

    y = @(x) g_M27 - (1 - x^2)/(gamma27*x^2) + ((gamma27 + 1)/(2*gamma27))*log(((gamma27+ 1)*x^2)/(2 + (gamma27 - 1)*x^2) );
    M27 = fsolve(y,0.006);
    
    T_star = T26/(0.5*(gamma26 + 1)/(1 + (gamma26 - 1)*0.5*M26^2));
    T27= T_star*(0.5*(gamma27+ 1)/(1 + (gamma27 - 1)*0.5*M27^2));
    
    P_star = P26/((1/M26)*sqrt(0.5*(gamma26 + 1)/(1 + (gamma26 - 1)*0.5*M26^2)));
    P27 = P_star*((1/M27)*sqrt(0.5*(gamma27+ 1)/(1 + (gamma27 - 1)*0.5*M27^2)));

    rho_star = rho26/((1/M26)*sqrt( 2*(1 + (gamma26 - 1)*0.5*M26^2)/(gamma26 + 1)));
    rho27 = rho_star*((1/M27)*sqrt( 2*(1 + (gamma27 - 1)*0.5*M27^2)/(gamma27 + 1)));
    
    c27 = sqrt(gamma27*R*T27);
    v27 = c27*M27;

    gamma27_new = gamma_N2(find(T==round(T27)),find(abs(P - round(P27,1)) < 0.001));

    err = abs(gamma27- gamma27_new);

    gamma27 = gamma27_new;

end

clear gamma27_new

%% Before mixing chamber (point 27) and after mixing chamber (point 28)
Texit = 680; 

rho_mc_entrance = rho_N2(find(T==round(T27)),find(abs(P - round(P27,1)) < 0.001)); 
v_mc_entrance = m_dot_N2/(A26_27*rho_mc_entrance);
Ptank = P27*1e5 - rho_mc_entrance*v_mc_entrance^2;   %[Pa]

T = [floor(Texit)-200:1:ceil(Texit)+5];
P = [floor(P27)-1:0.1:ceil(P27)];
data = nistdata('N2',T,P);
                                      %[bar]
rho_N2 = data.Rho*data.Mw; 
cp_N2 = data.Cp/data.Mw;             
cv_N2 = data.Cv/data.Mw;            
gamma_N2 = cp_N2./cv_N2; 

m_dot_SRP = 18*1e-3;

rho_mc_exit = rho_N2(find(T==round(Texit)),find(abs(P - round(Ptank*1e-5,1)) < 0.001)); 
m_dot_mix = m_dot_N2 + m_dot_SRP;
v_mc_exit = m_dot_mix/(A26_27*rho_mc_exit);
P_exit = Ptank - 0.5*rho_mc_exit*v_mc_exit^2;              %[Pa]
P28 = P_exit*1e-5;                                      %[bar]
T28 = Texit;
gamma28 = gamma_N2(find(T==round(T28)),find(abs(P - round(P28,1)) < 0.001));

v28 = v_mc_exit;
rho28 = rho_N2(find(T==round(T28)),find(abs(P - round(P28,1)) < 0.001)); 
c28 = (gamma28*R*T28)^0.5;
M28 = v28/c28; 


%% After cross fitting (point 29)

P29 = 1e-5*(P28*1e5 - 2*rho28*v28^2);
T29 = T28;

%% After cross fitting and before injector (point 29 and 30)

d29_30_ext = 19.05*1e-3;
t29_30 = 1.5*1e-3;
d29_30_int = d29_30_ext - 2*t29_30;
A29_30 = pi*(d29_30_int/2)^2;

eps29_30_rel = eps/d29_30_int;               % Relative roughness of stainless steel [-]


T = (floor(T29)-50):0.5:(ceil(T29));
P = (floor(P29)-2):0.1:(ceil(P29));
data = nistdata('N2',T,P);

rho_N2 = data.Rho*data.Mw;           % Density of Nitrogen [kg/m^3] 
cp_N2 = data.Cp/data.Mw;             % Specific heat at constant pressure of Nitrogen [J/kgK]
cv_N2 = data.Cv/data.Mw;             % Specific heat at constant volume of Nitrogen [J/kgK]
gamma_N2 = cp_N2./cv_N2;             % Ratio of specific heats [-]
mu_N2 = data.mu;                     % Viscosity of Nitrogen [Pa*s]

L29_30= 30*1e-2;
rho29 = rho_N2(find(T==round(T29)),find(abs(P - round(P29,1)) < 0.001));       % Density downstream the pipe bending after the pressure regulator [kg/m^3]
gamma29 = gamma_N2(find(T==round(T29)),find(abs(P - round(P29,1)) < 0.001));   % Ratio of specific heats downstream the pipe bending after the pressure regulator [-]
gamma30 =gamma_N2(find(T==round(T29)),find(abs(P - round(P29,1)) < 0.001));
mu29 = mu_N2(find(T==round(T29)),find(abs(P - round(P29,1)) < 0.001));         % Viscosity downstream the pipe bending after the pressure regulator [Pa*s]
v29 = m_dot_N2/(A29_30*rho29);                     % Gas velocity downstream the pipe bending after pressure regulator [m/s]
c29 = (gamma29*R*T29)^0.5;                         % Sound speed downstream the pipe bending after pressure regulator [m/s]
M29 = v29/c29;           
                  
Re29 = (rho29*v29*d29_30_int)/mu29;                    % Reynolds number downstream the manual ball valve [-]

if Re29< 2300

        lambda = 64/Re29;

    else 

        z = @(x) 1/sqrt(x) + 2*log10(2.51/(Re29*sqrt(x)) + eps29_30_rel/3.71);   % Colebrook-White correlation
        lambda = fsolve(z,0.0004);

end

iter = 0;
err = 1;

while err > 1e-3

    iter = iter + 1;

    g_M29 = (1 - M29^2)/(gamma29*M29^2) + ((gamma29 + 1)/(2*gamma29))*log(((gamma29 + 1)*M29^2)/(2 + (gamma29 - 1)*M29^2) );
    g_M30 = g_M29- lambda*(L29_30/d29_30_int);

    y = @(x) g_M30 - (1 - x^2)/(gamma30*x^2) + ((gamma30 + 1)/(2*gamma30))*log(((gamma30+ 1)*x^2)/(2 + (gamma30 - 1)*x^2) );
    M30 = fsolve(y,0.006);
    
    T_star = T29/(0.5*(gamma29 + 1)/(1 + (gamma29 - 1)*0.5*M29^2));
    T30= T_star*(0.5*(gamma30+ 1)/(1 + (gamma30 - 1)*0.5*M30^2));
    
    P_star = P29/((1/M29)*sqrt(0.5*(gamma29 + 1)/(1 + (gamma29 - 1)*0.5*M29^2)));
    P30 = P_star*((1/M30)*sqrt(0.5*(gamma30+ 1)/(1 + (gamma30 - 1)*0.5*M30^2)));

    rho_star = rho29/((1/M29)*sqrt( 2*(1 + (gamma29 - 1)*0.5*M29^2)/(gamma29 + 1)));
    rho30 = rho_star*((1/M30)*sqrt( 2*(1 + (gamma30 - 1)*0.5*M30^2)/(gamma30 + 1)));
    
    c30 = sqrt(gamma30*R*T30);
    v30 = c30*M30;

    gamma30_new = gamma_N2(find(T==round(T30)),find(abs(P - round(P30,1)) < 0.001));

    err = abs(gamma30- gamma30_new);

    gamma30 = gamma30_new;

end

clear gamma30_new

%% Injector

P_tot30 = P30*(1 + 0.5*(gamma30 - 1)*M30^2);
M_throat = 1;

z = @(x) x/A29_30 - (M30/M_throat)*sqrt( ((1 + 0.5*(gamma30 - 1)*M_throat^2)/(1 + 0.5*(gamma30 - 1)*M30^2))^((gamma30 + 1)/(gamma30 - 1)) );
A_throat_int = fsolve(z,0.08);
d_inj=sqrt(4*A_throat_int/pi)

P_chamber = P_tot30/(1 + 0.5*(gamma30 - 1)*M_throat^2)