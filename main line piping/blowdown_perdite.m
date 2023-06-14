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
P1 = 10.8;  

T = (floor(T1)-20):0.5:(ceil(T1));
P = (floor(P1)-8):0.1:(ceil(P1));

data = nistdata('N2',T,P);

rho_N2 = data.Rho*data.Mw;           % Density of Nitrogen [kg/m^3] 
cp_N2 = data.Cp/data.Mw;             % Specific heat at constant pressure of Nitrogen [J/kgK]
cv_N2 = data.Cv/data.Mw;             % Specific heat at constant volume of Nitrogen [J/kgK]
gamma_N2 = cp_N2./cv_N2;             % Ratio of specific heats [-]
mu_N2 = data.mu;                     % Viscosity of Nitrogen [Pa*s]
m_dot_N2 = 20*1e-3;                  % Nitrogen mass flow rate [kg/s]
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

T = (floor(T1)-3):0.5:(ceil(T1));
P = (floor(P2)-5):0.1:(ceil(P2));
data = nistdata('N2',T,P);

rho_N2 = data.Rho*data.Mw;           % Density of Nitrogen [kg/m^3] 
cp_N2 = data.Cp/data.Mw;             % Specific heat at constant pressure of Nitrogen [J/kgK]
cv_N2 = data.Cv/data.Mw;             % Specific heat at constant volume of Nitrogen [J/kgK]
gamma_N2 = cp_N2./cv_N2;             % Ratio of specific heats [-]
mu_N2 = data.mu; 

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


%% After pink elbow (point 4)

P4 = 1e-5*(P3*1e5 - 0.7*rho3*v3^2);

%% 12 mm tube before first T fitting (point 4 and 5)

d4_5_ext = 12*1e-3;
t4_5 = 1.5*1e-3;
d4_5_int = d4_5_ext - 2*t4_5;
A4_5 = pi*(d4_5_int/2)^2;

eps = 0.015*1e-3;                    % Absolute roughness of stainless steel [m]
eps4_5_rel = eps/d4_5_int;               % Relative roughness of stainless steel [-]

T4 = T3;
T = (floor(T4)-3):0.5:(ceil(T4));
P = (floor(P4)-4):0.1:(ceil(P4));
data = nistdata('N2',T,P);

rho_N2 = data.Rho*data.Mw;           % Density of Nitrogen [kg/m^3] 
cp_N2 = data.Cp/data.Mw;             % Specific heat at constant pressure of Nitrogen [J/kgK]
cv_N2 = data.Cv/data.Mw;             % Specific heat at constant volume of Nitrogen [J/kgK]
gamma_N2 = cp_N2./cv_N2;             % Ratio of specific heats [-]
mu_N2 = data.mu; 

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


%% After first T-fitting (point 6)

P6 = 1e-5*(P5*1e5 - 1.3*rho5*v5^2);


%% 12 mm tube after first T fitting (point 6 and 7)

d6_7_ext = 12*1e-3;
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


%% After manual ball valve (point 8)

G_g = rho7/1000;                      % Nitrogen specific gravity [-]
q_N2 = (m_dot_N2/rho7)*1000;           % Nitrogen volumetric flow rate [L/s]
C_V = 3.8;

P8 = P7 - (G_g*(q_N2*60)^2)/(14.42*C_V)^2;
T8 = T7;

%% Between manual ball valve and second T fitting (point 8 and 9)

d8_9_ext = 12*1e-3;
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



%% After second T-fitting (point 10)

P10 = 1e-5*(P9*1e5 - 1.3*rho9*v9^2);
T10 = T9;

%% After second T-fitting  (point 10) and before pink elbow (point 11) 

d10_11_ext = 12*1e-3;
t10_11 = 1.5*1e-3;
d10_11_int = d10_11_ext - 2*t10_11;
A10_11 = pi*(d10_11_int/2)^2;

eps10_11_rel = eps/d10_11_int;               % Relative roughness of stainless steel [-]


T = (floor(T10)-5):0.5:(ceil(T10));
P = (floor(P10)-3):0.1:(ceil(P10));
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

%% After pink elbow (point 12)

P12 = 1e-5*(P11*1e5 - 0.7*rho11*v11^2);
T12=T11;

%% Between pink elbow and pneumatic valve (point 12 and 13)

d12_13_ext = 12*1e-3;
t12_13 = 1.5*1e-3;
d12_13_int = d12_13_ext - 2*t12_13;
A12_13 = pi*(d12_13_int/2)^2;

eps12_13_rel = eps/d12_13_int;               % Relative roughness of stainless steel [-]


T = (floor(T12)-5):0.5:(ceil(T12));
P = (floor(P12)-3):0.1:(ceil(P12));
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

%% Before penumatic valve (point 13) and after pneumatic valve (point 14)
G_g = rho13/1000;                      % Nitrogen specific gravity [-]
q_N2 = (m_dot_N2/rho13)*1000;           % Nitrogen volumetric flow rate [L/s]

C_V = 1.68;                             % Flow coefficient check valve
P14 = P13 - (G_g*(q_N2*60)^2)/(14.42*C_V)^2;         % Pressure downstream the check valve [bar]
T14 = T13;


%% After penumatic valve and before check valve (point 14 and 15)

d14_15_ext = 12*1e-3;
t14_15 = 1.5*1e-3;
d14_15_int = d14_15_ext - 2*t14_15;
A14_15 = pi*(d14_15_int/2)^2;

eps14_15_rel = eps/d14_15_int;               % Relative roughness of stainless steel [-]


T = (floor(T14)-5):0.5:(ceil(T14));
P = (floor(P14)-3):0.1:(ceil(P14));
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

%% Before check valve (point 15) and after check valve (point 16)
G_g = rho15/1000;                      % Nitrogen specific gravity [-]
q_N2 = (m_dot_N2/rho15)*1000;           % Nitrogen volumetric flow rate [L/s]

C_V = 1.68;                             % Flow coefficient check valve
P16 = P15 - (G_g*(q_N2*60)^2)/(14.42*C_V)^2;         % Pressure downstream the check valve [bar]
T16 = T15;



%% After check valve and before blow down (point 16 and 17)

d16_17_ext = 12*1e-3;
t16_17 = 1.5*1e-3;
d16_17_int = d16_17_ext - 2*t16_17;
A16_17 = pi*(d16_17_int/2)^2;

eps16_17_rel = eps/d16_17_int;               % Relative roughness of stainless steel [-]


T = (floor(T16)-5):0.5:(ceil(T16));
P = (floor(P16)-3):0.1:(ceil(P16));
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
