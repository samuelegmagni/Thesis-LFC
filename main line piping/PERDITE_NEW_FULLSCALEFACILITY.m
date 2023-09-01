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
P1 = 29;  

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

m_dot_N2 = 70*1e-3;                  % Nitrogen mass flow rate [kg/s]
R_N2 = 8314/28;                         % Specific ideal gas constant [J/kgK]

d1_ext = 6.35*1e-3;                   % Pipe external diameter [m]
t1 = 1*1e-3;                        % Thickness of the tube  [m]
d1_int = d1_ext - 2*t1;             % Pipe internal diameter [m]
A1_int = pi*(d1_int/2)^2;            % Internal cross sectional area [m^2]

rho1 = rho_N2( find(abs(T - round(T1,1))==min(abs(T - round(T1,1)))) ,find( abs(P - round(P1,1))==min(abs(P - round(P1,1)))) ); % Density downstream the pressure regulator [kg/m^3]
v1 = m_dot_N2/(A1_int*rho1);                                             % Velocity downstream the pressure regulator [m/s]
P2 = 1e-5*(P1*1e5 - 1*rho1*v1^2);                                  % Pressure downstream the pipe bending after the pressure regulator [bar]
gamma1 = gamma_N2( find(abs(T - round(T1,1))==min(abs(T - round(T1,1)))) ,find( abs(P - round(P1,1))==min(abs(P - round(P1,1)))) );  % Ratio of specific heats downstream the pipe bending after the pressure regulator [-]
mu1 = mu_N2( find(abs(T - round(T1,1))==min(abs(T - round(T1,1)))) ,find( abs(P - round(P1,1))==min(abs(P - round(P1,1)))) );       % Viscosity downstream the pipe bending after the pressure regulator [Pa*s]
c1 = (gamma1*R_N2*T1)^0.5;                         % Sound speed downstream the pipe bending after pressure regulator [m/s]
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
rho2 = rho_N2( find(abs(T - round(T2,1))==min(abs(T - round(T2,1)))) ,find( abs(P - round(P2,1))==min(abs(P - round(P2,1)))) );      % Density downstream the pipe bending after the pressure regulator [kg/m^3]
gamma2 = gamma_N2( find(abs(T - round(T2,1))==min(abs(T - round(T2,1)))) ,find( abs(P - round(P2,1))==min(abs(P - round(P2,1)))) );  % Ratio of specific heats downstream the pipe bending after the pressure regulator [-]
mu2 = mu_N2( find(abs(T - round(T2,1))==min(abs(T - round(T2,1)))) ,find( abs(P - round(P2,1))==min(abs(P - round(P2,1)))) );       % Viscosity downstream the pipe bending after the pressure regulator [Pa*s]
v2 = m_dot_N2/(A2_3*rho2);                     % Gas velocity downstream the pipe bending after pressure regulator [m/s]
c2 = (gamma2*R_N2*T2)^0.5;                         % Sound speed downstream the pipe bending after pressure regulator [m/s]
M2 = v2/c2;                                     % Mach number downstream the pipe bending after pressure regulator [-]
Re2 = (rho2*v2*d2_3_int)/mu2;

L2_3 = 5*1e-2;
gamma3 = gamma_N2( find(abs(T - round(T2,1))==min(abs(T - round(T2,1)))) ,find( abs(P - round(P2,1))==min(abs(P - round(P2,1)))) );  % Ratio of specific heats downstream the pipe bending after the pressure regulator [-]

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
    
    c3 = sqrt(gamma3*R_N2*T3);
    v3 = M3*c3;

    gamma3_new = gamma_N2( find(abs(T - round(T3,1))==min(abs(T - round(T3,1)))) ,find( abs(P - round(P3,1))==min(abs(P - round(P3,1)))) );

    err = abs(gamma3 - gamma3_new);

    gamma3 = gamma3_new;

end

clear gamma3_new

%% After green fitting (12 mm -> 3/4") (point 4)

P4 = 1e-5*(P3*1e5 - 1*rho3*v3^2);
T4 = T3;

%% 3/4" mm tube before blue elbow (point 4 and 5)

d4_5_ext = 19.05*1e-3;
t4_5 = 1.5*1e-3;
d4_5_int = d4_5_ext - 2*t4_5;
A4_5 = pi*(d4_5_int/2)^2;

eps = 0.015*1e-3;                    % Absolute roughness of stainless steel [m]
eps4_5_rel = eps/d4_5_int;               % Relative roughness of stainless steel [-]

rho4 = rho_N2( find(abs(T - round(T4,1))==min(abs(T - round(T4,1)))) ,find( abs(P - round(P4,1))==min(abs(P - round(P4,1)))) );      % Density downstream the pipe bending after the pressure regulator [kg/m^3]
gamma4 = gamma_N2( find(abs(T - round(T4,1))==min(abs(T - round(T4,1)))) ,find( abs(P - round(P4,1))==min(abs(P - round(P4,1)))) );   % Ratio of specific heats downstream the pipe bending after the pressure regulator [-]
mu4 = mu_N2( find(abs(T - round(T4,1))==min(abs(T - round(T4,1)))) ,find( abs(P - round(P4,1))==min(abs(P - round(P4,1)))) );          % Viscosity downstream the pipe bending after the pressure regulator [Pa*s]
v4 = m_dot_N2/(A4_5*rho4);                     % Gas velocity downstream the pipe bending after pressure regulator [m/s]
c4 = (gamma4*R_N2*T4)^0.5;                         % Sound speed downstream the pipe bending after pressure regulator [m/s]
M4 = v4/c4;                                     % Mach number downstream the pipe bending after pressure regulator [-]
Re4 = (rho4*v4*d4_5_int)/mu4;

L4_5 = 5*1e-2;
gamma5 = gamma_N2( find(abs(T - round(T4,1))==min(abs(T - round(T4,1)))) ,find( abs(P - round(P4,1))==min(abs(P - round(P4,1)))) );   % Ratio of specific heats downstream the pipe bending after the pressure regulator [-]

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
    
    c5 = sqrt(gamma5*R_N2*T5);
    v5 = M5*c5;

    gamma5_new = gamma_N2( find(abs(T - round(T5,1))==min(abs(T - round(T5,1)))) ,find( abs(P - round(P5,1))==min(abs(P - round(P5,1)))) ); 

    err = abs(gamma5 - gamma5_new);

    gamma5 = gamma5_new;

end

clear gamma5_new

%% After blue elbow (point 6)

P6 = 1e-5*(P5*1e5 - 0.7*rho5*v5^2);
T6 = T5;

%% 3/4" mm tube after blue elbow (point 6 and 7)

d6_7_ext = 19.05*1e-3;
t6_7 = 1.5*1e-3;
d6_7_int = d6_7_ext - 2*t6_7;
A6_7 = pi*(d6_7_int/2)^2;

eps = 0.015*1e-3;                    % Absolute roughness of stainless steel [m]
eps6_7_rel = eps/d6_7_int;               % Relative roughness of stainless steel [-]

L6_7 = 5*1e-2;
gamma6 = gamma_N2( find(abs(T - round(T6,1))==min(abs(T - round(T6,1)))) ,find( abs(P - round(P6,1))==min(abs(P - round(P6,1)))) );   % Ratio of specific heats downstream the pipe bending after the pressure regulator [-]
gamma7 = gamma_N2( find(abs(T - round(T6,1))==min(abs(T - round(T6,1)))) ,find( abs(P - round(P6,1))==min(abs(P - round(P6,1)))) );    % Ratio of specific heats downstream the pipe bending after the pressure regulator [-]
mu6 = mu_N2( find(abs(T - round(T6,1))==min(abs(T - round(T6,1)))) ,find( abs(P - round(P6,1))==min(abs(P - round(P6,1)))) );    % Ratio of specific heats downstream the pipe bending after the pressure regulator [-]
rho6 = rho_N2( find(abs(T - round(T6,1))==min(abs(T - round(T6,1)))) ,find( abs(P - round(P6,1))==min(abs(P - round(P6,1)))) );    % Ratio of specific heats downstream the pipe bending after the pressure regulator [-]

v6 = m_dot_N2/(A6_7*rho6);                     % Gas velocity downstream the pipe bending after pressure regulator [m/s]
c6 = (gamma6*R_N2*T6)^0.5;                         % Sound speed downstream the pipe bending after pressure regulator [m/s]
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
    
    c7 = sqrt(gamma7*R_N2*T7);
    v7 = M7*c7;

    gamma7_new = gamma_N2( find(abs(T - round(T7,1))==min(abs(T - round(T7,1)))) ,find( abs(P - round(P7,1))==min(abs(P - round(P7,1)))) ); 

    err = abs(gamma7 - gamma7_new);

    gamma7 = gamma7_new;

end

clear gamma7_new

%% After first T-fitting (point 8)

P8 = 1e-5*(P7*1e5 - 1.3*rho7*v7^2);
T8 = T7;

%% Between T-fitting and manual ball valve (point 8 and 9)

d8_9_ext = 19.05*1e-3;
t8_9 = 1.5*1e-3;
d8_9_int = d8_9_ext - 2*t8_9;
A8_9 = pi*(d8_9_int/2)^2;

eps = 0.015*1e-3;                    % Absolute roughness of stainless steel [m]
eps8_9_rel = eps/d8_9_int;               % Relative roughness of stainless steel [-]

L8_9 = 5*1e-2;
gamma8 = gamma_N2( find(abs(T - round(T8,1))==min(abs(T - round(T8,1)))) ,find( abs(P - round(P8,1))==min(abs(P - round(P8,1)))) );   % Ratio of specific heats downstream the pipe bending after the pressure regulator [-]
gamma9 = gamma_N2( find(abs(T - round(T8,1))==min(abs(T - round(T8,1)))) ,find( abs(P - round(P8,1))==min(abs(P - round(P8,1)))) );  % Ratio of specific heats downstream the pipe bending after the pressure regulator [-]
mu8 = mu_N2( find(abs(T - round(T8,1))==min(abs(T - round(T8,1)))) ,find( abs(P - round(P8,1))==min(abs(P - round(P8,1)))) );  % Ratio of specific heats downstream the pipe bending after the pressure regulator [-]
rho8 = rho_N2( find(abs(T - round(T8,1))==min(abs(T - round(T8,1)))) ,find( abs(P - round(P8,1))==min(abs(P - round(P8,1)))) );   % Ratio of specific heats downstream the pipe bending after the pressure regulator [-]

v8 = m_dot_N2/(A8_9*rho8);                     % Gas velocity downstream the pipe bending after pressure regulator [m/s]
c8 = (gamma8*R_N2*T8)^0.5;                         % Sound speed downstream the pipe bending after pressure regulator [m/s]
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
    
    c9 = sqrt(gamma9*R_N2*T9);
    v9 = M9*c9;

    gamma9_new = gamma_N2( find(abs(T - round(T9,1))==min(abs(T - round(T9,1)))) ,find( abs(P - round(P9,1))==min(abs(P - round(P9,1)))) );

    err = abs(gamma9 - gamma9_new);

    gamma9 = gamma9_new;

end

clear gamma9_new

%% After manual ball valve (point 10)

G_g = rho9/1000;                      % Nitrogen specific gravity [-]
q_N2 = (m_dot_N2/rho9)*1000;           % Nitrogen volumetric flow rate [L/s]
% q_N2_std = (P4*q_N2*T4*60)/(1*273.15); % Nitrogen volumetric flow rate at std conditions [std L/min]
C_V = 7.1;

P10 = P9 - (G_g*(q_N2*60)^2)/(14.42*C_V)^2;
T10 = T9;

%% After first manual ball valve (point 10) and before second T-fitting (point 11) 

d10_11_ext = 19.05*1e-3;
t10_11 = 1.5*1e-3;
d10_11_int = d10_11_ext - 2*t10_11;
A10_11 = pi*(d10_11_int/2)^2;

eps10_11_rel = eps/d10_11_int;               % Relative roughness of stainless steel [-]

L10_11 = 5*1e-2;
rho10 = rho_N2( find(abs(T - round(T10,1))==min(abs(T - round(T10,1)))) ,find( abs(P - round(P10,1))==min(abs(P - round(P10,1)))) );     % Density downstream the manual ball valve [kg/m^3]
gamma10 = gamma_N2( find(abs(T - round(T10,1))==min(abs(T - round(T10,1)))) ,find( abs(P - round(P10,1))==min(abs(P - round(P10,1)))) ); % Ratio of specific heats downstream the manual ball valve  [-]
gamma11 = gamma_N2( find(abs(T - round(T10,1))==min(abs(T - round(T10,1)))) ,find( abs(P - round(P10,1))==min(abs(P - round(P10,1)))) );
mu10 = mu_N2( find(abs(T - round(T10,1))==min(abs(T - round(T10,1)))) ,find( abs(P - round(P10,1))==min(abs(P - round(P10,1)))) );       % Viscosity downstream the manual ball valve [Pa*s]
v10 = m_dot_N2/(A10_11*rho10);                     % Gas velocity downstream the manual ball valve [m/s]
c10 = (gamma10*R_N2*T10)^0.5;                         % Sound speed downstream the manual ball valve [m/s]
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
    
    c11 = sqrt(gamma11*R_N2*T11);
    v11 = c11*M11;

    gamma11_new = gamma_N2( find(abs(T - round(T11,1))==min(abs(T - round(T11,1)))) ,find( abs(P - round(P11,1))==min(abs(P - round(P11,1)))) );

    err = abs(gamma11 - gamma11_new);

    gamma11 = gamma11_new;

end

clear gamma11_new

%% After second T-fitting (point 12)

P12 = 1e-5*(P11*1e5 - 1.3*rho11*v11^2);
T12 = T11;

%% Between second and elbow (point 12 and 13)

d12_13_ext = 19.05*1e-3;
t12_13 = 1.5*1e-3;
d12_13_int = d12_13_ext - 2*t12_13;
A12_13 = pi*(d12_13_int/2)^2;

eps12_13_rel = eps/d12_13_int;               % Relative roughness of stainless steel [-]

L12_13 = 5*1e-2;
rho12 = rho_N2( find(abs(T - round(T12,1))==min(abs(T - round(T12,1)))) ,find( abs(P - round(P12,1))==min(abs(P - round(P12,1)))) );    % Density downstream the manual ball valve [kg/m^3]
gamma12 = gamma_N2( find(abs(T - round(T12,1))==min(abs(T - round(T12,1)))) ,find( abs(P - round(P12,1))==min(abs(P - round(P12,1)))) ); % Ratio of specific heats downstream the manual ball valve  [-]
gamma13 = gamma_N2( find(abs(T - round(T12,1))==min(abs(T - round(T12,1)))) ,find( abs(P - round(P12,1))==min(abs(P - round(P12,1)))) );
mu12 = mu_N2( find(abs(T - round(T12,1))==min(abs(T - round(T12,1)))) ,find( abs(P - round(P12,1))==min(abs(P - round(P12,1)))) );       % Viscosity downstream the manual ball valve [Pa*s]
v12 = m_dot_N2/(A12_13*rho12);                     % Gas velocity downstream the manual ball valve [m/s]
c12 = (gamma12*R_N2*T12)^0.5;                         % Sound speed downstream the manual ball valve [m/s]
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
    
    c13 = sqrt(gamma13*R_N2*T13);
    v13 = c13*M13;

    gamma13_new = gamma_N2( find(abs(T - round(T13,1))==min(abs(T - round(T13,1)))) ,find( abs(P - round(P13,1))==min(abs(P - round(P13,1)))) );

    err = abs(gamma13 - gamma13_new);

    gamma13 = gamma13_new;

end

clear gamma13_new

%% After elbow (point 14)

P14 = 1e-5*(P13*1e5 - 0.7*rho13*v13^2);
T14 = T13;

%% After check valve (point 15)

rho14 = rho_N2( find(abs(T - round(T14,1))==min(abs(T - round(T14,1)))) ,find( abs(P - round(P14,1))==min(abs(P - round(P14,1)))) );    % Density downstream the manual ball valve [kg/m^3]

G_g = rho14/1000;                      % Nitrogen specific gravity [-]
q_N2 = (m_dot_N2/rho14)*1000;           % Nitrogen volumetric flow rate [L/s]

C_V = 2.2;                             % Flow coefficient check valve
P15 = P14 - (G_g*(q_N2*60)^2)/(14.42*C_V)^2;         % Pressure downstream the check valve [bar]
T15 = T14;

%% After check valve and before T-fitting (point 15 and 16)

d15_16_ext = 19.05*1e-3;
t15_16 = 1.5*1e-3;
d15_16_int = d15_16_ext - 2*t15_16;
A15_16 = pi*(d15_16_int/2)^2;

eps15_16_rel = eps/d15_16_int;               % Relative roughness of stainless steel [-]

L15_16 = 5*1e-2;
rho15 = rho_N2( find(abs(T - round(T15,1))==min(abs(T - round(T15,1)))) ,find( abs(P - round(P15,1))==min(abs(P - round(P15,1)))) );    % Density downstream the manual ball valve [kg/m^3]
gamma15 = gamma_N2( find(abs(T - round(T15,1))==min(abs(T - round(T15,1)))) ,find( abs(P - round(P15,1))==min(abs(P - round(P15,1)))) ); % Ratio of specific heats downstream the manual ball valve  [-]
gamma16 = gamma_N2( find(abs(T - round(T15,1))==min(abs(T - round(T15,1)))) ,find( abs(P - round(P15,1))==min(abs(P - round(P15,1)))) );
mu15 = mu_N2( find(abs(T - round(T15,1))==min(abs(T - round(T15,1)))) ,find( abs(P - round(P15,1))==min(abs(P - round(P15,1)))) );       % Viscosity downstream the manual ball valve [Pa*s]
v15 = m_dot_N2/(A15_16*rho15);                     % Gas velocity downstream the manual ball valve [m/s]
c15 = (gamma15*R_N2*T15)^0.5;                         % Sound speed downstream the manual ball valve [m/s]
M15 = v15/c15;                                     % Mach number downstream the manual ball valve [-]
Re15 = (rho15*v15*d15_16_int)/mu15;                    % Reynolds number downstream the manual ball valve [-]

if Re15 < 2300

        lambda = 64/Re15;

    else 

        z = @(x) 1/sqrt(x) + 2*log10(2.51/(Re15*sqrt(x)) + eps15_16_rel/3.71);   % Colebrook-White correlation
        lambda = fsolve(z,0.0004);

end

iter = 0;
err = 1;

while err > 1e-3

    iter = iter + 1;

    g_M15 = (1 - M15^2)/(gamma15*M15^2) + ((gamma15 + 1)/(2*gamma15))*log(((gamma15 + 1)*M15^2)/(2 + (gamma15 - 1)*M15^2) );
    g_M16 = g_M15 - lambda*(L15_16/d15_16_int);

    y = @(x) g_M16 - (1 - x^2)/(gamma16*x^2) + ((gamma16 + 1)/(2*gamma16))*log(((gamma16 + 1)*x^2)/(2 + (gamma16 - 1)*x^2) );
    M16 = fsolve(y,0.006);
    
    T_star = T15/(0.5*(gamma15 + 1)/(1 + (gamma15 - 1)*0.5*M15^2));
    T16 = T_star*(0.5*(gamma16 + 1)/(1 + (gamma16 - 1)*0.5*M16^2));
    
    P_star = P15/((1/M15)*sqrt(0.5*(gamma15 + 1)/(1 + (gamma15 - 1)*0.5*M15^2)));
    P16 = P_star*((1/M16)*sqrt(0.5*(gamma16 + 1)/(1 + (gamma16 - 1)*0.5*M16^2)));

    rho_star = rho15/((1/M15)*sqrt( 2*(1 + (gamma15 - 1)*0.5*M15^2)/(gamma15 + 1)));
    rho16 = rho_star*((1/M16)*sqrt( 2*(1 + (gamma16 - 1)*0.5*M16^2)/(gamma16 + 1)));
    
    c16 = sqrt(gamma16*R_N2*T16);
    v16 = c16*M16;

    gamma16_new = gamma_N2( find(abs(T - round(T16,1))==min(abs(T - round(T16,1)))) ,find( abs(P - round(P16,1))==min(abs(P - round(P16,1)))) );

    err = abs(gamma16 - gamma16_new);

    gamma16 = gamma16_new;

end

clear gamma16_new

%% After T-fitting that merges the two lines (point 17 and 18)

P17 = 1e-5*(P16*1e5 - 2.5*rho16*v16^2);
T17 = T16;

m_dot_new = 2*m_dot_N2;

d17_18_ext = 19.05*1e-3;
t17_18 = 1.5*1e-3;
d17_18_int = d17_18_ext - 2*t17_18;
A17_18 = pi*(d17_18_int/2)^2;

eps17_18_rel = eps/d17_18_int;               % Relative roughness of stainless steel [-]

L17_18 = 15*1e-2;
rho17 = rho_N2( find(abs(T - round(T17,1))==min(abs(T - round(T17,1)))) ,find( abs(P - round(P17,1))==min(abs(P - round(P17,1)))) );    % Density downstream the manual ball valve [kg/m^3]
gamma17 = gamma_N2( find(abs(T - round(T17,1))==min(abs(T - round(T17,1)))) ,find( abs(P - round(P17,1))==min(abs(P - round(P17,1)))) );  % Ratio of specific heats downstream the manual ball valve  [-]
gamma18 = gamma_N2( find(abs(T - round(T17,1))==min(abs(T - round(T17,1)))) ,find( abs(P - round(P17,1))==min(abs(P - round(P17,1)))) ); 
mu17 = mu_N2( find(abs(T - round(T17,1))==min(abs(T - round(T17,1)))) ,find( abs(P - round(P17,1))==min(abs(P - round(P17,1)))) );        % Viscosity downstream the manual ball valve [Pa*s]
v17 = m_dot_new/(A17_18*rho17);                     % Gas velocity downstream the manual ball valve [m/s]
c17 = (gamma17*R_N2*T17)^0.5;                         % Sound speed downstream the manual ball valve [m/s]
M17 = v17/c17;                                     % Mach number downstream the manual ball valve [-]
Re17 = (rho17*v17*d17_18_int)/mu17;                    % Reynolds number downstream the manual ball valve [-]

if Re17 < 2300

        lambda = 64/Re17;

    else 

        z = @(x) 1/sqrt(x) + 2*log10(2.51/(Re17*sqrt(x)) + eps17_18_rel/3.71);   % Colebrook-White correlation
        lambda = fsolve(z,0.0004);

end

iter = 0;
err = 1;

while err > 1e-3

    iter = iter + 1;

    g_M17 = (1 - M17^2)/(gamma17*M17^2) + ((gamma17 + 1)/(2*gamma17))*log(((gamma17 + 1)*M17^2)/(2 + (gamma17 - 1)*M17^2) );
    g_M18 = g_M17 - lambda*(L17_18/d17_18_int);

    y = @(x) g_M18 - (1 - x^2)/(gamma18*x^2) + ((gamma18 + 1)/(2*gamma18))*log(((gamma18 + 1)*x^2)/(2 + (gamma18 - 1)*x^2) );
    M18 = fsolve(y,0.006);
    
    T_star = T17/(0.5*(gamma17 + 1)/(1 + (gamma17 - 1)*0.5*M17^2));
    T18 = T_star*(0.5*(gamma18 + 1)/(1 + (gamma18 - 1)*0.5*M18^2));
    
    P_star = P17/((1/M17)*sqrt(0.5*(gamma17 + 1)/(1 + (gamma17 - 1)*0.5*M17^2)));
    P18 = P_star*((1/M18)*sqrt(0.5*(gamma18 + 1)/(1 + (gamma18 - 1)*0.5*M18^2)));

    rho_star = rho17/((1/M17)*sqrt( 2*(1 + (gamma17 - 1)*0.5*M17^2)/(gamma17 + 1)));
    rho18 = rho_star*((1/M18)*sqrt( 2*(1 + (gamma18 - 1)*0.5*M18^2)/(gamma18 + 1)));
    
    c18 = sqrt(gamma18*R_N2*T18);
    v18 = c18*M18;

    gamma18_new = gamma_N2( find(abs(T - round(T18,1))==min(abs(T - round(T18,1)))) ,find( abs(P - round(P18,1))==min(abs(P - round(P18,1)))) );

    err = abs(gamma18 - gamma18_new);

    gamma18 = gamma18_new;

end

clear gamma18_new

%% After T-fitting and before Venturi (point 19 and 20)

P19 = 1e-5*(P18*1e5 - 1.3*rho18*v18^2);
T19 = T18;

d19_20_ext = 19.05*1e-3;
t19_20 = 1.5*1e-3;
d19_20_int = d19_20_ext - 2*t19_20;
A19_20 = pi*(d19_20_int/2)^2;

eps19_20_rel = eps/d19_20_int;               % Relative roughness of stainless steel [-]

L19_20 = 5*1e-2;
rho19 = rho_N2( find(abs(T - round(T19,1))==min(abs(T - round(T19,1)))) ,find( abs(P - round(P19,1))==min(abs(P - round(P19,1)))) );    % Density downstream the manual ball valve [kg/m^3]
gamma19 = gamma_N2( find(abs(T - round(T19,1))==min(abs(T - round(T19,1)))) ,find( abs(P - round(P19,1))==min(abs(P - round(P19,1)))) );  % Ratio of specific heats downstream the manual ball valve  [-]
gamma20 = gamma_N2( find(abs(T - round(T19,1))==min(abs(T - round(T19,1)))) ,find( abs(P - round(P19,1))==min(abs(P - round(P19,1)))) ); 
mu19 = mu_N2( find(abs(T - round(T19,1))==min(abs(T - round(T19,1)))) ,find( abs(P - round(P19,1))==min(abs(P - round(P19,1)))) );       % Viscosity downstream the manual ball valve [Pa*s]
v19 = m_dot_new/(A19_20*rho19);                     % Gas velocity downstream the manual ball valve [m/s]
c19 = (gamma19*R_N2*T19)^0.5;                         % Sound speed downstream the manual ball valve [m/s]
M19 = v19/c19;                                     % Mach number downstream the manual ball valve [-]
Re19 = (rho19*v19*d19_20_int)/mu19;                    % Reynolds number downstream the manual ball valve [-]

if Re19 < 2300

        lambda = 64/Re19;

    else 

        z = @(x) 1/sqrt(x) + 2*log10(2.51/(Re19*sqrt(x)) + eps19_20_rel/3.71);   % Colebrook-White correlation
        lambda = fsolve(z,0.0004);

end

iter = 0;
err = 1;

while err > 1e-3

    iter = iter + 1;

    g_M19 = (1 - M19^2)/(gamma19*M19^2) + ((gamma19 + 1)/(2*gamma19))*log(((gamma19 + 1)*M19^2)/(2 + (gamma19 - 1)*M19^2) );
    g_M20 = g_M19 - lambda*(L19_20/d19_20_int);

    y = @(x) g_M20 - (1 - x^2)/(gamma20*x^2) + ((gamma20 + 1)/(2*gamma20))*log(((gamma20 + 1)*x^2)/(2 + (gamma20 - 1)*x^2) );
    M20 = fsolve(y,0.006);
    
    T_star = T19/(0.5*(gamma19 + 1)/(1 + (gamma19 - 1)*0.5*M19^2));
    T20 = T_star*(0.5*(gamma20 + 1)/(1 + (gamma20 - 1)*0.5*M20^2));
    
    P_star = P19/((1/M19)*sqrt(0.5*(gamma19 + 1)/(1 + (gamma19 - 1)*0.5*M19^2)));
    P20 = P_star*((1/M20)*sqrt(0.5*(gamma20 + 1)/(1 + (gamma20 - 1)*0.5*M20^2)));

    rho_star = rho19/((1/M19)*sqrt( 2*(1 + (gamma19 - 1)*0.5*M19^2)/(gamma19 + 1)));
    rho20 = rho_star*((1/M20)*sqrt( 2*(1 + (gamma20 - 1)*0.5*M20^2)/(gamma20 + 1)));
    
    c20 = sqrt(gamma20*R_N2*T20);
    v20 = c20*M20;

    gamma20_new = gamma_N2( find(abs(T - round(T20,1))==min(abs(T - round(T20,1)))) ,find( abs(P - round(P20,1))==min(abs(P - round(P20,1)))) );

    err = abs(gamma20 - gamma20_new);

    gamma20 = gamma20_new;

end

clear gamma20_new


%% Across venturi channel (point 20 and 21)


d_throat_int = 7.5*1e-3;
A_throat_int = 0.25*pi*d_throat_int^2;

d21_22_ext = 19.05*1e-3;
t21_22 = 1.5*1e-3;
d21_22_int = d21_22_ext - 2*t21_22;
A21_22 = pi*(d21_22_int/2)^2;

T_tot = T20*(1 + ((gamma20 - 1)/2)*M20^2);
P_tot20 = P20*(1 + ((gamma20 - 1)/2)*M20^2)^(gamma20/(gamma20 - 1));

z = @(x) A_throat_int/A21_22 - (M20/x)*sqrt( ((1 + 0.5*(gamma20 - 1)*x^2)/(1 + 0.5*(gamma20 - 1)*M20^2))^((gamma20 + 1)/(gamma20 - 1)) );
M20_1 = fsolve(z,0.8);
T20_1 = T_tot/(1 + ((gamma20 - 1)/2)*M20_1^2);
P20_1 = P_tot20/(1 + ((gamma20 - 1)/2)*M20_1^2)^(gamma20/(gamma20 - 1));
c20_1 = sqrt(gamma20*R_N2*T20_1);
v20_1 = c20_1*M20_1;
rho20_1 = (rho20*v20*A21_22)/(A_throat_int*v20_1);

mu20_1 = mu_N2( find(abs(T - round(T20_1,1))==min(abs(T - round(T20_1,1)))) ,find( abs(P - round(P20_1,1))==min(abs(P - round(P20_1,1)))) );
gamma20_1 = gamma_N2( find(abs(T - round(T20_1,1))==min(abs(T - round(T20_1,1)))) ,find( abs(P - round(P20_1,1))==min(abs(P - round(P20_1,1)))) );
gamma20_2 = gamma_N2( find(abs(T - round(T20_1,1))==min(abs(T - round(T20_1,1)))) ,find( abs(P - round(P20_1,1))==min(abs(P - round(P20_1,1)))) );
Re20_1 = (rho20_1*v20_1*d_throat_int)/mu20_1;

eps20_1rel = eps/d_throat_int; 

if Re20_1 < 2300

        lambda = 64/Re20_1;

    else 

        z = @(x) 1/sqrt(x) + 2*log10(2.51/(Re20_1*sqrt(x)) + eps20_1rel/3.71);   % Colebrook-White correlation
        lambda = fsolve(z,0.0004);

end

M20_2 = 1;
iter = 0;
err = 1;

while err > 1e-3

    iter = iter + 1;

    g_M20_1 = (1 - M20_1^2)/(gamma20_1*M20_1^2) + ((gamma20_1 + 1)/(2*gamma20_1))*log(((gamma20_1 + 1)*M20_1^2)/(2 + (gamma20_1 - 1)*M20_1^2) );
    L_throat = 1000*(g_M20_1*d_throat_int)/lambda;       % [mm]

    T_star = T20_1/(0.5*(gamma20_1 + 1)/(1 + (gamma20_1 - 1)*0.5*M20_1^2));
    T20_2 = T_star*(0.5*(gamma20_2 + 1)/(1 + (gamma20_2 - 1)*0.5*M20_2^2));

    P_star = P20_1/((1/M20_1)*sqrt(0.5*(gamma20_1 + 1)/(1 + (gamma20_1 - 1)*0.5*M20_1^2)));
    P20_2 = P_star*((1/M20_2)*sqrt(0.5*(gamma20_2 + 1)/(1 + (gamma20_2 - 1)*0.5*M20_2^2)));

    P_tot_star = (M20_1*P_tot20)/( ((1 + 0.5*(gamma20_1 - 1)*M20_1^2)/(0.5*(gamma20_1 + 1)))^( (0.5*(gamma20_1 + 1))/(gamma20_1 - 1) ));
    P_tot20_2 = (P_tot_star/M20_2)*( ((1 + 0.5*(gamma20_2 - 1)*M20_2^2)/(0.5*(gamma20_2 + 1)))^( (0.5*(gamma20_2 + 1))/(gamma20_2 - 1) ));

    rho_star = rho20_1/((1/M20_1)*sqrt( 2*(1 + (gamma20_1 - 1)*0.5*M20_1^2)/(gamma20_1 + 1)));
    rho20_2 = rho_star*((1/M20_2)*sqrt( 2*(1 + (gamma20_2 - 1)*0.5*M20_2^2)/(gamma20_2 + 1)));

    c20_2 = sqrt(gamma20_2*R_N2*T20_2);
    v20_2 = M20_2*c20_2;

    gamma20_2_new = gamma_N2( find(abs(T - round(T20_2,1))==min(abs(T - round(T20_2,1)))) ,find( abs(P - round(P20_2,1))==min(abs(P - round(P20_2,1)))) );

    err = abs(gamma20_2 - gamma20_2_new);

    gamma20_2 = gamma20_2_new;

end

clear gamma20_2_new

z = @(x) A21_22/A_throat_int - (M20_2/x)*sqrt( ((1 + 0.5*(gamma20_2 - 1)*x^2)/(1 + 0.5*(gamma20_2 - 1)*M20_2^2))^((gamma20_2 + 1)/(gamma20_2 - 1)) );
M21 = fsolve(z,0.8);

P21 = P_tot20_2/(1 + ((gamma20_2 - 1)/2)*M21^2)^(gamma20_2/(gamma20_2 - 1));
T21 = T_tot/(1 + ((gamma20_2 - 1)/2)*M21^2);

%% After Venturi channel, before servovalve (point 21 and 22)

eps21_22_rel = eps/d21_22_int;               % Relative roughness of stainless steel [-]

L21_22 = 30*1e-2;
rho21 = rho_N2(find(abs(T - round(T21,1))==min(abs(T - round(T21,1)))) ,find( abs(P - round(P21,1))==min(abs(P - round(P21,1)))) );       % Density downstream the pipe bending after the pressure regulator [kg/m^3]
gamma21 = gamma_N2(find(abs(T - round(T21,1))==min(abs(T - round(T21,1)))) ,find( abs(P - round(P21,1))==min(abs(P - round(P21,1)))) );   % Ratio of specific heats downstream the pipe bending after the pressure regulator [-]
gamma22 = gamma_N2(find(abs(T - round(T21,1))==min(abs(T - round(T21,1)))) ,find( abs(P - round(P21,1))==min(abs(P - round(P21,1)))) );
mu21 = mu_N2(find(abs(T - round(T21,1))==min(abs(T - round(T21,1)))) ,find( abs(P - round(P21,1))==min(abs(P - round(P21,1)))) );       % Viscosity downstream the pipe bending after the pressure regulator [Pa*s]
v21 = m_dot_new/(A21_22*rho21);                     % Gas velocity downstream the pipe bending after pressure regulator [m/s]
c21 = (gamma21*R_N2*T21)^0.5;           % Sound speed downstream the pipe bending after pressure regulator [m/s]       

x1 = 0;
x2 = 0.0348*1e3;
x3 = x2 + L_throat;
x4 = x3 + 0.0489*1000;
d_vect = [x1 x2 x3 x4];
v_vect = [v20 v20_1 v20_2 v21];
P_vect = [P20 P20_1 P20_2 P21];
T_vect = [T20 T20_1 T20_2 T21];
rho_vect = [rho20 rho20_1 rho20_2 rho21];
M_vect = [M20 M20_1 M20_2 M21];


figure()
plot(d_vect(1),v_vect(1),'ro','linewidth',1.5)
grid on
hold on
plot(d_vect(2),v_vect(2),'bo','linewidth',1.5)
plot(d_vect(3),v_vect(3),'go','linewidth',1.5)
plot(d_vect(4),v_vect(4),'ko','linewidth',1.5)
legend('$v_{in,conv}$','$v_{in,throat}$','$v_{out,throat}$','$v_{out,div}$')
xlabel('Position, $x_i$ $[mm]$')
ylabel('Velocity, $v_i$ $[m/s]$')

figure()
plot(d_vect(1),M_vect(1),'ro','linewidth',1.5)
grid on
hold on
plot(d_vect(2),M_vect(2),'bo','linewidth',1.5)
plot(d_vect(3),M_vect(3),'go','linewidth',1.5)
plot(d_vect(4),M_vect(4),'ko','linewidth',1.5)
legend('$M_{in,conv}$','$M_{in,throat}$','$M_{out,throat}$','$M_{out,div}$')
xlabel('Position, $x_i$ $[mm]$')
ylabel('Mach Number, $M_i$ $[-]$')

figure()
plot(d_vect(1),P_vect(1),'ro','linewidth',1.5)
grid on
hold on
plot(d_vect(2),P_vect(2),'bo','linewidth',1.5)
plot(d_vect(3),P_vect(3),'go','linewidth',1.5)
plot(d_vect(4),P_vect(4),'ko','linewidth',1.5)
legend('$P_{in,conv}$','$P_{in,throat}$','$P_{out,throat}$','$P_{out,div}$')
xlabel('Position, $x_i$ $[mm]$')
ylabel('Pressure, $P_i$ $[bar]$')


figure()
plot(d_vect(1),T_vect(1),'ro','linewidth',1.5)
grid on
hold on
plot(d_vect(2),T_vect(2),'bo','linewidth',1.5)
plot(d_vect(3),T_vect(3),'go','linewidth',1.5)
plot(d_vect(4),T_vect(4),'ko','linewidth',1.5)
legend('$T_{in,conv}$','$T_{in,throat}$','$T_{out,throat}$','$T_{out,div}$')
xlabel('Position, $x_i$ $[mm]$')
ylabel('Temperature, $T_i$ $[K]$')

figure()
plot(d_vect(1),rho_vect(1),'ro','linewidth',1.5)
grid on
hold on
plot(d_vect(2),rho_vect(2),'bo','linewidth',1.5)
plot(d_vect(3),rho_vect(3),'go','linewidth',1.5)
plot(d_vect(4),rho_vect(4),'ko','linewidth',1.5)
legend('$\rho_{in,conv}$','$\rho_{in,throat}$','$\rho_{out,throat}$','$\rho_{out,div}$')
xlabel('Position, $x_i$ $[mm]$')
ylabel('Density, $\rho_i$ $[kg/m^3]$')



Re21 = (rho21*v21*d21_22_int)/mu21;                    % Reynolds number downstream the manual ball valve [-]

if Re21 < 2300

        lambda = 64/Re21;

    else 

        z = @(x) 1/sqrt(x) + 2*log10(2.51/(Re21*sqrt(x)) + eps21_22_rel/3.71);   % Colebrook-White correlation
        lambda = fsolve(z,0.0004);

end

iter = 0;
err = 1;

while err > 1e-3

    iter = iter + 1;

    g_M21 = (1 - M21^2)/(gamma21*M21^2) + ((gamma21 + 1)/(2*gamma21))*log(((gamma21 + 1)*M21^2)/(2 + (gamma21 - 1)*M21^2) );
    g_M22 = g_M21 - lambda*(L21_22/d21_22_int);

    y = @(x) g_M22 - (1 - x^2)/(gamma22*x^2) + ((gamma22 + 1)/(2*gamma22))*log(((gamma22 + 1)*x^2)/(2 + (gamma22 - 1)*x^2) );
    M22 = fsolve(y,0.006);
    
    T_star = T21/(0.5*(gamma21 + 1)/(1 + (gamma21 - 1)*0.5*M21^2));
    T22 = T_star*(0.5*(gamma22 + 1)/(1 + (gamma22 - 1)*0.5*M22^2));
    
    P_star = P21/((1/M21)*sqrt(0.5*(gamma21 + 1)/(1 + (gamma21 - 1)*0.5*M21^2)));
    P22 = P_star*((1/M22)*sqrt(0.5*(gamma22 + 1)/(1 + (gamma22 - 1)*0.5*M22^2)));

    rho_star = rho21/((1/M21)*sqrt( 2*(1 + (gamma21 - 1)*0.5*M21^2)/(gamma21 + 1)));
    rho22 = rho_star*((1/M22)*sqrt( 2*(1 + (gamma22 - 1)*0.5*M22^2)/(gamma22 + 1)));
    
    c22 = sqrt(gamma22*R_N2*T22);
    v22 = c22*M22;

    gamma22_new = gamma_N2(find(abs(T - round(T22,1))==min(abs(T - round(T22,1)))) ,find( abs(P - round(P22,1))==min(abs(P - round(P22,1)))) );

    err = abs(gamma22 - gamma22_new);

    gamma22 = gamma22_new;

end

clear gamma22_new

%% Before servovalve (point 22) and after servovalve (point 23)

G_g = rho22/1000;                      % Nitrogen specific gravity [-]
q_N2 = (m_dot_new/rho22)*1000;           % Nitrogen volumetric flow rate [L/s]
C_V = 7.1;                             % Flow coefficient ball valve

P23 = P22 - (G_g*(q_N2*60)^2)/(14.42*C_V)^2;     % Pressure downstream the servovalve (ball valve approx) [bar]
T23 = T22;                        % Temperature downstream the servovalve (ball valve approx) [K]

%% After servovalve and before check valve (point 23 and 24)

d23_24_ext = 19.05*1e-3;
t23_24 = 1.5*1e-3;
d23_24_int = d23_24_ext - 2*t23_24;
A23_24 = pi*(d23_24_int/2)^2;

eps23_24_rel = eps/d23_24_int;               % Relative roughness of stainless steel [-]

L23_24 = 30*1e-2;
rho23 = rho_N2(find(abs(T - round(T23,1))==min(abs(T - round(T23,1)))) ,find( abs(P - round(P23,1))==min(abs(P - round(P23,1)))) );       % Density downstream the pipe bending after the pressure regulator [kg/m^3]
gamma23 = gamma_N2(find(abs(T - round(T23,1))==min(abs(T - round(T23,1)))) ,find( abs(P - round(P23,1))==min(abs(P - round(P23,1)))) );    % Ratio of specific heats downstream the pipe bending after the pressure regulator [-]
gamma24 = gamma_N2(find(abs(T - round(T23,1))==min(abs(T - round(T23,1)))) ,find( abs(P - round(P23,1))==min(abs(P - round(P23,1)))) ); 
mu23 = mu_N2(find(abs(T - round(T23,1))==min(abs(T - round(T23,1)))) ,find( abs(P - round(P23,1))==min(abs(P - round(P23,1)))) );        % Viscosity downstream the pipe bending after the pressure regulator [Pa*s]
v23 = m_dot_new/(A23_24*rho23);                     % Gas velocity downstream the pipe bending after pressure regulator [m/s]
c23 = (gamma23*R_N2*T23)^0.5;                         % Sound speed downstream the pipe bending after pressure regulator [m/s]
M23 = v23/c23;           
                  
Re23 = (rho23*v23*d23_24_int)/mu23;                    % Reynolds number downstream the manual ball valve [-]

if Re23 < 2300

        lambda = 64/Re23;

    else 

        z = @(x) 1/sqrt(x) + 2*log10(2.51/(Re23*sqrt(x)) + eps23_24_rel/3.71);   % Colebrook-White correlation
        lambda = fsolve(z,0.0004);

end

iter = 0;
err = 1;

while err > 1e-3

    iter = iter + 1;

    g_M23 = (1 - M23^2)/(gamma23*M23^2) + ((gamma23 + 1)/(2*gamma23))*log(((gamma23 + 1)*M23^2)/(2 + (gamma23 - 1)*M23^2) );
    g_M24 = g_M23 - lambda*(L23_24/d23_24_int);

    y = @(x) g_M24 - (1 - x^2)/(gamma24*x^2) + ((gamma24 + 1)/(2*gamma24))*log(((gamma24 + 1)*x^2)/(2 + (gamma24 - 1)*x^2) );
    M24 = fsolve(y,0.006);
    
    T_star = T23/(0.5*(gamma23 + 1)/(1 + (gamma23 - 1)*0.5*M23^2));
    T24 = T_star*(0.5*(gamma24 + 1)/(1 + (gamma24 - 1)*0.5*M24^2));
    
    P_star = P23/((1/M23)*sqrt(0.5*(gamma23 + 1)/(1 + (gamma23 - 1)*0.5*M23^2)));
    P24 = P_star*((1/M24)*sqrt(0.5*(gamma24 + 1)/(1 + (gamma24 - 1)*0.5*M24^2)));

    rho_star = rho23/((1/M23)*sqrt( 2*(1 + (gamma23 - 1)*0.5*M23^2)/(gamma23 + 1)));
    rho24 = rho_star*((1/M24)*sqrt( 2*(1 + (gamma24 - 1)*0.5*M24^2)/(gamma24 + 1)));
    
    c24 = sqrt(gamma24*R_N2*T24);
    v24 = c24*M24;

    gamma24_new = gamma_N2(find(abs(T - round(T24,1))==min(abs(T - round(T24,1)))) ,find( abs(P - round(P24,1))==min(abs(P - round(P24,1)))) );

    err = abs(gamma24 - gamma24_new);

    gamma24 = gamma24_new;

end

clear gamma24_new

%% After check valve (point 25)

G_g = rho24/1000;                      % Nitrogen specific gravity [-]
q_N2 = (m_dot_new/rho24)*1000;           % Nitrogen volumetric flow rate [L/s]

C_V = 2.2;                             % Flow coefficient check valve
P25 = P24 - (G_g*(q_N2*60)^2)/(14.42*C_V)^2;         % Pressure downstream the check valve [bar]
T25 = T24;

%% After check valve and before cross fitting (point 25 and 26)

d25_26_ext = 19.05*1e-3;
t25_26 = 1.5*1e-3;
d25_26_int = d25_26_ext - 2*t25_26;
A25_26 = pi*(d25_26_int/2)^2;

eps25_26_rel = eps/d25_26_int;               % Relative roughness of stainless steel [-]

L25_26 = 5*1e-2;
rho25 = rho_N2(find(abs(T - round(T25,1))==min(abs(T - round(T25,1)))) ,find( abs(P - round(P25,1))==min(abs(P - round(P25,1)))) );       % Density downstream the pipe bending after the pressure regulator [kg/m^3]
gamma25 = gamma_N2(find(abs(T - round(T25,1))==min(abs(T - round(T25,1)))) ,find( abs(P - round(P25,1))==min(abs(P - round(P25,1)))) );    % Ratio of specific heats downstream the pipe bending after the pressure regulator [-]
gamma26 = gamma_N2(find(abs(T - round(T25,1))==min(abs(T - round(T25,1)))) ,find( abs(P - round(P25,1))==min(abs(P - round(P25,1)))) ); 
mu25 = mu_N2(find(abs(T - round(T25,1))==min(abs(T - round(T25,1)))) ,find( abs(P - round(P25,1))==min(abs(P - round(P25,1)))) );        % Viscosity downstream the pipe bending after the pressure regulator [Pa*s]
v25 = m_dot_new/(A25_26*rho23);                     % Gas velocity downstream the pipe bending after pressure regulator [m/s]
c25 = (gamma25*R_N2*T25)^0.5;                         % Sound speed downstream the pipe bending after pressure regulator [m/s]
M25 = v25/c25;           
                  
Re25 = (rho25*v25*d25_26_int)/mu25;                    % Reynolds number downstream the manual ball valve [-]

if Re25 < 2300

        lambda = 64/Re25;

    else 

        z = @(x) 1/sqrt(x) + 2*log10(2.51/(Re25*sqrt(x)) + eps25_26_rel/3.71);   % Colebrook-White correlation
        lambda = fsolve(z,0.0004);

end

iter = 0;
err = 1;

while err > 1e-3

    iter = iter + 1;

    g_M25 = (1 - M25^2)/(gamma25*M25^2) + ((gamma25 + 1)/(2*gamma25))*log(((gamma25 + 1)*M25^2)/(2 + (gamma25 - 1)*M25^2) );
    g_M26 = g_M25 - lambda*(L25_26/d25_26_int);

    y = @(x) g_M26 - (1 - x^2)/(gamma26*x^2) + ((gamma26 + 1)/(2*gamma26))*log(((gamma26 + 1)*x^2)/(2 + (gamma26 - 1)*x^2) );
    M26 = fsolve(y,0.006);
    
    T_star = T25/(0.5*(gamma25 + 1)/(1 + (gamma25 - 1)*0.5*M25^2));
    T26 = T_star*(0.5*(gamma26 + 1)/(1 + (gamma26 - 1)*0.5*M26^2));
    
    P_star = P25/((1/M25)*sqrt(0.5*(gamma25 + 1)/(1 + (gamma25 - 1)*0.5*M25^2)));
    P26 = P_star*((1/M26)*sqrt(0.5*(gamma26 + 1)/(1 + (gamma26 - 1)*0.5*M26^2)));

    rho_star = rho25/((1/M25)*sqrt( 2*(1 + (gamma25 - 1)*0.5*M25^2)/(gamma25 + 1)));
    rho26 = rho_star*((1/M26)*sqrt( 2*(1 + (gamma26 - 1)*0.5*M26^2)/(gamma26 + 1)));
    
    c26 = sqrt(gamma26*R_N2*T26);
    v26 = c26*M26;

    gamma26_new = gamma_N2(find(abs(T - round(T26,1))==min(abs(T - round(T26,1)))) ,find( abs(P - round(P26,1))==min(abs(P - round(P26,1)))) );

    err = abs(gamma26 - gamma26_new);

    gamma26 = gamma26_new;

end

clear gamma26_new

%% After cross fitting (point 27)

P27 = 1e-5*(P26*1e5 - 2*rho26*v26^2);
T27 = T26;

%% After cross fitting and before the injector (point 27 and 28)

d27_28_ext = 19.05*1e-3;
t27_28 = 1.5*1e-3;
d27_28_int = d27_28_ext - 2*t27_28;
A27_28 = pi*(d27_28_int/2)^2;

eps27_28_rel = eps/d27_28_int;               % Relative roughness of stainless steel [-]

L27_28 = 5*1e-2;
rho27 = rho_N2(find(abs(T - round(T27,1))==min(abs(T - round(T27,1)))) ,find( abs(P - round(P27,1))==min(abs(P - round(P27,1)))) );       % Density downstream the pipe bending after the pressure regulator [kg/m^3]
gamma27 = gamma_N2(find(abs(T - round(T27,1))==min(abs(T - round(T27,1)))) ,find( abs(P - round(P27,1))==min(abs(P - round(P27,1)))) );    % Ratio of specific heats downstream the pipe bending after the pressure regulator [-]
gamma28 = gamma_N2(find(abs(T - round(T27,1))==min(abs(T - round(T27,1)))) ,find( abs(P - round(P27,1))==min(abs(P - round(P27,1)))) ); 
mu27 = mu_N2(find(abs(T - round(T27,1))==min(abs(T - round(T27,1)))) ,find( abs(P - round(P27,1))==min(abs(P - round(P27,1)))) );        % Viscosity downstream the pipe bending after the pressure regulator [Pa*s]
v27 = m_dot_new/(A27_28*rho23);                     % Gas velocity downstream the pipe bending after pressure regulator [m/s]
c27 = (gamma27*R_N2*T27)^0.5;                         % Sound speed downstream the pipe bending after pressure regulator [m/s]
M27 = v27/c27;           
                  
Re27 = (rho27*v27*d27_28_int)/mu27;                    % Reynolds number downstream the manual ball valve [-]

if Re27 < 2300

        lambda = 64/Re27;

    else 

        z = @(x) 1/sqrt(x) + 2*log10(2.51/(Re27*sqrt(x)) + eps27_28_rel/3.71);   % Colebrook-White correlation
        lambda = fsolve(z,0.0004);

end

iter = 0;
err = 1;

while err > 1e-3

    iter = iter + 1;

    g_M27 = (1 - M27^2)/(gamma27*M27^2) + ((gamma27 + 1)/(2*gamma27))*log(((gamma27 + 1)*M27^2)/(2 + (gamma27 - 1)*M27^2) );
    g_M28 = g_M27 - lambda*(L27_28/d27_28_int);

    y = @(x) g_M28 - (1 - x^2)/(gamma28*x^2) + ((gamma28 + 1)/(2*gamma28))*log(((gamma28 + 1)*x^2)/(2 + (gamma28 - 1)*x^2) );
    M28 = fsolve(y,0.006);
    
    T_star = T27/(0.5*(gamma27 + 1)/(1 + (gamma27 - 1)*0.5*M27^2));
    T28 = T_star*(0.5*(gamma28 + 1)/(1 + (gamma28 - 1)*0.5*M28^2));
    
    P_star = P27/((1/M27)*sqrt(0.5*(gamma27 + 1)/(1 + (gamma27 - 1)*0.5*M27^2)));
    P28 = P_star*((1/M28)*sqrt(0.5*(gamma28 + 1)/(1 + (gamma28 - 1)*0.5*M28^2)));

    rho_star = rho27/((1/M27)*sqrt( 2*(1 + (gamma27 - 1)*0.5*M27^2)/(gamma27 + 1)));
    rho28 = rho_star*((1/M28)*sqrt( 2*(1 + (gamma28 - 1)*0.5*M28^2)/(gamma28 + 1)));
    
    c28 = sqrt(gamma28*R_N2*T28);
    v28 = c28*M28;

    gamma28_new = gamma_N2(find(abs(T - round(T28,1))==min(abs(T - round(T28,1)))) ,find( abs(P - round(P28,1))==min(abs(P - round(P28,1)))) );

    err = abs(gamma28 - gamma28_new);

    gamma28 = gamma28_new;

end

clear gamma28_new