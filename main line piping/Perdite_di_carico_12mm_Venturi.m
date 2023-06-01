%% Pressure drops in N2 line with 12 mm diameter tubes

clc
clear


set(0,'DefaultTextFontSize',12);              % Settings for the plot
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultLegendInterpreter','Latex');
set(groot,'DefaultAxesTickLabelInterpreter','Latex');
set(0,'DefaultTextInterpreter','Latex');
set(0,'DefaultLegendFontSize',12);

d_p_ext = 12*1e-3;                   % Pipe external diameter [m]
t = 1.5*1e-3;                        % Thickness of the tube  [m]
d_p_int = d_p_ext - 2*t;             % Pipe internal diameter [m]
A_int = pi*(d_p_int/2)^2;            % Internal cross sectional area [m^2]

eps = 0.015*1e-3;                    % Absolute roughness of stainless steel [m]
eps_rel = eps/d_p_int;               % Relative roughness of stainless steel [-]

T1 = 298;                                       % Temperature downstream the pressure regulator [K]
P_reg = 42;                                     % Pressure downstream the pressure regulator [bar]

T = (floor(T1)-3):0.5:(ceil(T1));
P = (floor(P_reg)-10):0.1:(ceil(P_reg));
data = nistdata('N2',T,P);

rho_N2 = data.Rho*data.Mw;           % Density of Nitrogen [kg/m^3] 
cp_N2 = data.Cp/data.Mw;             % Specific heat at constant pressure of Nitrogen [J/kgK]
cv_N2 = data.Cv/data.Mw;             % Specific heat at constant volume of Nitrogen [J/kgK]
gamma_N2 = cp_N2./cv_N2;             % Ratio of specific heats [-]
mu_N2 = data.mu;                     % Viscosity of Nitrogen [Pa*s]
m_dot_N2 = 140*1e-3;                  % Nitrogen mass flow rate [kg/s]
R = 8314/28;                         % Specific ideal gas constant [J/kgK]


L1 = 0.4;                             % Length of the tube [m]

%% After pressure regulator (point 1) and before first manual ball valve (point 2)

rho_reg = rho_N2(find(T==round(T1)),find(abs(P - round(P_reg,1)) < 0.001)); % Density downstream the pressure regulator [kg/m^3]
v_reg = m_dot_N2/(A_int*rho_reg);                                             % Velocity downstream the pressure regulator [m/s]
P1 = 1e-5*(P_reg*1e5 - 1.3*rho_reg*v_reg^2);                                  % Pressure downstream the pipe bending after the pressure regulator [bar]
rho1 = rho_N2(find(T==round(T1)),find(abs(P - round(P1,1)) < 0.001));       % Density downstream the pipe bending after the pressure regulator [kg/m^3]
gamma1 = gamma_N2(find(T==round(T1)),find(abs(P - round(P1,1)) < 0.001));   % Ratio of specific heats downstream the pipe bending after the pressure regulator [-]
gamma2 = gamma_N2(find(T==round(T1)),find(abs(P - round(P1,1)) < 0.001));
mu1 = mu_N2(find(T==round(T1)),find(abs(P - round(P1,1)) < 0.001));         % Viscosity downstream the pipe bending after the pressure regulator [Pa*s]
v1 = m_dot_N2/(A_int*rho1);                     % Gas velocity downstream the pipe bending after pressure regulator [m/s]
c1 = (gamma1*R*T1)^0.5;                         % Sound speed downstream the pipe bending after pressure regulator [m/s]
M1 = v1/c1;                                     % Mach number downstream the pipe bending after pressure regulator [-]
Re1 = (rho1*v1*d_p_int)/mu1;                    % Reynolds number downstream the pipe bending after pressure regulator [-]

if Re1 < 2300

        lambda = 64/Re1;

    else 

        z = @(x) 1/sqrt(x) + 2*log10(2.51/(Re1*sqrt(x)) + eps_rel/3.71);   % Colebrook-White correlation
        lambda = fsolve(z,0.0004);

end

iter = 0;
err = 1;

while err > 1e-3

    iter = iter + 1;

    g_M1 = (1 - M1^2)/(gamma1*M1^2) + ((gamma1 + 1)/(2*gamma1))*log(((gamma1 + 1)*M1^2)/(2 + (gamma1 - 1)*M1^2) );
    g_M2 = g_M1 - lambda*(L1/d_p_int);
    
    y = @(x) g_M2 - (1 - x^2)/(gamma2*x^2) + ((gamma2 + 1)/(2*gamma2))*log(((gamma2 + 1)*x^2)/(2 + (gamma2 - 1)*x^2) );
    M2 = fsolve(y,0.006);
    
    T_star = T1/(0.5*(gamma1 + 1)/(1 + (gamma1 - 1)*0.5*M1^2));
    T2 = T_star*(0.5*(gamma2 + 1)/(1 + (gamma2 - 1)*0.5*M2^2));
    
    P_star = P1/((1/M1)*sqrt(0.5*(gamma1 + 1)/(1 + (gamma1 - 1)*0.5*M1^2)));
    P2 = P_star*((1/M2)*sqrt(0.5*(gamma2 + 1)/(1 + (gamma2 - 1)*0.5*M2^2)));
    
    rho_star = rho1/((1/M1)*sqrt( 2*(1 + (gamma1 - 1)*0.5*M1^2)/(gamma1 + 1)));
    rho2 = rho_star*((1/M2)*sqrt( 2*(1 + (gamma2 - 1)*0.5*M2^2)/(gamma2 + 1)));
    
    c2 = sqrt(gamma2*R*T2);
    v2 = M2*c2;

    gamma2_new = gamma_N2(find(T==round(T2)),find(abs(P - round(P2,1)) < 0.001));

    err = abs(gamma2 - gamma2_new);

    gamma2 = gamma2_new;

end

clear gamma2_new

%% After first manual ball valve (point 3) and before second manual ball valve (point 4)

P3 = 1e-5*(P2*1e5 - 0.1*rho2*v2^2);                % Pressure drop related to the T fitting (dump line)  [bar]
T3 = T2;

T = (floor(T3)-10):0.5:(ceil(T3));
P = (floor(P3)-9):0.1:(ceil(P3));
data = nistdata('N2',T,P);

rho_N2 = data.Rho*data.Mw;           % Density of Nitrogen [kg/m^3] 
cp_N2 = data.Cp/data.Mw;             % Specific heat at constant pressure of Nitrogen [J/kgK]
cv_N2 = data.Cv/data.Mw;             % Specific heat at constant volume of Nitrogen [J/kgK]
gamma_N2 = cp_N2./cv_N2;             % Ratio of specific heats [-]
mu_N2 = data.mu;                     % Viscosity of Nitrogen [Pa*s]

rho3 = rho_N2(find(T==round(T3)),find(abs(P - round(P3,1)) < 0.001));     % Density downstream the pressure regulator [kg/m^3]
gamma3 = gamma_N2(find(T==round(T3)),find(abs(P - round(P3,1)) < 0.001)); % Ratio of specific heats downstream the pressure regulator [-]
gamma4 = gamma_N2(find(T==round(T3)),find(abs(P - round(P3,1)) < 0.001));
mu3 = mu_N2(find(T==round(T3)),find(abs(P - round(P3,1)) < 0.001));       % Viscosity downstream the pressure regulator [Pa*s]
v3 = m_dot_N2/(A_int*rho3);                     % Gas velocity downstream the pressure regulator [m/s]
c3 = (gamma3*R*T3)^0.5;                         % Sound speed downstream the pressure regulator [m/s]
M3 = v3/c3;                                     % Mach number downstream the pressure regulator [-]
Re3 = (rho3*v3*d_p_int)/mu3;                    % Reynolds number downstream the pressure regulator [-]


if Re3 < 2300

        lambda = 64/Re3;

    else 

        z = @(x) 1/sqrt(x) + 2*log10(2.51/(Re3*sqrt(x)) + eps_rel/3.71);   % Colebrook-White correlation
        lambda = fsolve(z,0.0004);

end

iter = 0;
err = 1;

while err > 1e-3

    iter = iter + 1;

    g_M3 = (1 - M3^2)/(gamma3*M3^2) + ((gamma3 + 1)/(2*gamma3))*log(((gamma3 + 1)*M3^2)/(2 + (gamma3 - 1)*M3^2) );
    g_M4 = g_M3 - lambda*(L1/d_p_int);
    
    y = @(x) g_M4 - (1 - x^2)/(gamma4*x^2) + ((gamma4 + 1)/(2*gamma4))*log(((gamma4 + 1)*x^2)/(2 + (gamma4 - 1)*x^2) );
    M4 = fsolve(y,0.006);
    
    T_star = T3/(0.5*(gamma3 + 1)/(1 + (gamma3 - 1)*0.5*M3^2));
    T4 = T_star*(0.5*(gamma4 + 1)/(1 + (gamma4 - 1)*0.5*M4^2));
    
    P_star = P3/((1/M3)*sqrt(0.5*(gamma3 + 1)/(1 + (gamma3 - 1)*0.5*M3^2)));
    P4 = P_star*((1/M4)*sqrt(0.5*(gamma4 + 1)/(1 + (gamma4 - 1)*0.5*M4^2)));
    
    rho_star = rho3/((1/M3)*sqrt( 2*(1 + (gamma3 - 1)*0.5*M3^2)/(gamma3 + 1)));
    rho4 = rho_star*((1/M4)*sqrt( 2*(1 + (gamma4 - 1)*0.5*M4^2)/(gamma4 + 1)));
    
    c4 = sqrt(gamma4*R*T4);
    v4 = M4*c4;

    gamma4_new = gamma_N2(find(T==round(T4)),find(abs(P - round(P4,1)) < 0.001));

    err = abs(gamma4 - gamma4_new);

    gamma4 = gamma4_new;

end

clear gamma4_new

%% Before second manual ball valve (point 4) and after second manual ball valve (point 5)
G_g = rho4/1000;                      % Nitrogen specific gravity [-]
q_N2 = (m_dot_N2/rho4)*1000;           % Nitrogen volumetric flow rate [L/s]
% q_N2_std = (P4*q_N2*T4*60)/(1*273.15); % Nitrogen volumetric flow rate at std conditions [std L/min]
C_V = 3.8;

P5 = P4 - (G_g*(q_N2*60)^2)/(14.42*C_V)^2;
T5 = T4;

%% After second manual ball valve (point 5) and before Venturi channel (point 6) 
T = (floor(T5)-10):0.5:(ceil(T5));
P = (floor(P5)-9):0.1:(ceil(P5));
data = nistdata('N2',T,P);

rho_N2 = data.Rho*data.Mw;           % Density of Nitrogen [kg/m^3] 
cp_N2 = data.Cp/data.Mw;             % Specific heat at constant pressure of Nitrogen [J/kgK]
cv_N2 = data.Cv/data.Mw;             % Specific heat at constant volume of Nitrogen [J/kgK]
gamma_N2 = cp_N2./cv_N2;             % Ratio of specific heats [-]
mu_N2 = data.mu;                     % Viscosity of Nitrogen [Pa*s]


rho5 = rho_N2(find(T==round(T5)),find(abs(P - round(P5,1)) < 0.001));     % Density downstream the manual ball valve [kg/m^3]
gamma5 = gamma_N2(find(T==round(T5)),find(abs(P - round(P5,1)) < 0.001)); % Ratio of specific heats downstream the manual ball valve  [-]
gamma6 = gamma_N2(find(T==round(T5)),find(abs(P - round(P5,1)) < 0.001));
mu5 = mu_N2(find(T==round(T5)),find(abs(P - round(P5,1)) < 0.001));       % Viscosity downstream the manual ball valve [Pa*s]
v5 = m_dot_N2/(A_int*rho5);                     % Gas velocity downstream the manual ball valve [m/s]
c5 = (gamma5*R*T5)^0.5;                         % Sound speed downstream the manual ball valve [m/s]
M5 = v5/c5;                                     % Mach number downstream the manual ball valve [-]
Re5 = (rho5*v5*d_p_int)/mu5;                    % Reynolds number downstream the manual ball valve [-]

if Re5 < 2300

        lambda = 64/Re5;

    else 

        z = @(x) 1/sqrt(x) + 2*log10(2.51/(Re5*sqrt(x)) + eps_rel/3.71);   % Colebrook-White correlation
        lambda = fsolve(z,0.0004);

end

iter = 0;
err = 1;

while err > 1e-3

    iter = iter + 1;

    g_M5 = (1 - M5^2)/(gamma5*M5^2) + ((gamma5 + 1)/(2*gamma5))*log(((gamma5 + 1)*M5^2)/(2 + (gamma5 - 1)*M5^2) );
    g_M6 = g_M5 - lambda*(L1/d_p_int);

    y = @(x) g_M6 - (1 - x^2)/(gamma6*x^2) + ((gamma6 + 1)/(2*gamma6))*log(((gamma6 + 1)*x^2)/(2 + (gamma6 - 1)*x^2) );
    M6 = fsolve(y,0.006);
    
    T_star = T5/(0.5*(gamma5 + 1)/(1 + (gamma5 - 1)*0.5*M5^2));
    T6 = T_star*(0.5*(gamma6 + 1)/(1 + (gamma6 - 1)*0.5*M6^2));
    
    P_star = P5/((1/M5)*sqrt(0.5*(gamma5 + 1)/(1 + (gamma5 - 1)*0.5*M5^2)));
    P6 = P_star*((1/M6)*sqrt(0.5*(gamma6 + 1)/(1 + (gamma6 - 1)*0.5*M6^2)));

    rho_star = rho5/((1/M5)*sqrt( 2*(1 + (gamma5 - 1)*0.5*M5^2)/(gamma5 + 1)));
    rho6 = rho_star*((1/M6)*sqrt( 2*(1 + (gamma6 - 1)*0.5*M6^2)/(gamma6 + 1)));
    
    c6 = sqrt(gamma6*R*T6);
    v6 = c6*M6;

    gamma6_new = gamma_N2(find(T==round(T6)),find(abs(P - round(P6,1)) < 0.001));

    err = abs(gamma6 - gamma6_new);

    gamma6 = gamma6_new;

end

clear gamma6_new

%% Before Venturi channel (point 6) and after Venturi channel (point 7)

d_throat_int = 5.5*1e-3;
A_throat_int = 0.25*pi*d_throat_int^2;

T_tot = T6*(1 + ((gamma6 - 1)/2)*M6^2);
P_tot6 = P6*(1 + ((gamma6 - 1)/2)*M6^2)^(gamma6/(gamma6 - 1));

z = @(x) A_throat_int/A_int - (M6/x)*sqrt( ((1 + 0.5*(gamma6 - 1)*x^2)/(1 + 0.5*(gamma6 - 1)*M6^2))^((gamma6 + 1)/(gamma6 - 1)) );
M6_1 = fsolve(z,0.8);
T6_1 = T_tot/(1 + ((gamma6 - 1)/2)*M6_1^2);
P6_1 = P_tot6/(1 + ((gamma6 - 1)/2)*M6_1^2)^(gamma6/(gamma6 - 1));
c6_1 = sqrt(gamma6*R*T6_1);
v6_1 = c6_1*M6_1;
rho6_1 = (rho6*v6*A_int)/(A_throat_int*v6_1);

T = (floor(T6_1)-35):0.5:(ceil(T6_1));
P = (floor(P6_1)-10):0.1:(ceil(P6_1));

data = nistdata('N2',T,P);

cp_N2 = data.Cp/data.Mw;             % Specific heat at constant pressure of Nitrogen [J/kgK]
cv_N2 = data.Cv/data.Mw;             % Specific heat at constant volume of Nitrogen [J/kgK]
gamma_N2 = cp_N2./cv_N2;             % Ratio of specific heats [-]
mu_N2 = data.mu;                     % Viscosity of Nitrogen [Pa*s]

mu6_1 = mu_N2(find(T==round(T6_1)),find(abs(P - round(P6_1,1)) < 0.001));
gamma6_1 = gamma_N2(find(T==round(T6_1)),find(abs(P - round(P6_1,1)) < 0.001));
gamma6_2 = gamma_N2(find(T==round(T6_1)),find(abs(P - round(P6_1,1)) < 0.001));
Re6_1 = (rho6_1*v6_1*d_throat_int)/mu6_1;

if Re6_1 < 2300

        lambda = 64/Re6_1;

    else 

        z = @(x) 1/sqrt(x) + 2*log10(2.51/(Re6_1*sqrt(x)) + eps_rel/3.71);   % Colebrook-White correlation
        lambda = fsolve(z,0.0004);

end

M6_2 = 1;
iter = 0;
err = 1;

while err > 1e-3

    iter = iter + 1;

    g_M6_1 = (1 - M6_1^2)/(gamma6_1*M6_1^2) + ((gamma6_1 + 1)/(2*gamma6_1))*log(((gamma6_1 + 1)*M6_1^2)/(2 + (gamma6_1 - 1)*M6_1^2) );
    L_throat = 1000*(g_M6_1*d_throat_int)/lambda;

    T_star = T6_1/(0.5*(gamma6_1 + 1)/(1 + (gamma6_1 - 1)*0.5*M6_1^2));
    T6_2 = T_star*(0.5*(gamma6_2 + 1)/(1 + (gamma6_2 - 1)*0.5*M6_2^2));

    P_star = P6_1/((1/M6_1)*sqrt(0.5*(gamma6_1 + 1)/(1 + (gamma6_1 - 1)*0.5*M6_1^2)));
    P6_2 = P_star*((1/M6_2)*sqrt(0.5*(gamma6_2 + 1)/(1 + (gamma6_2 - 1)*0.5*M6_2^2)));

    P_tot_star = (M6_1*P_tot6)/( ((1 + 0.5*(gamma6_1 - 1)*M6_1^2)/(0.5*(gamma6_1 + 1)))^( (0.5*(gamma6_1 + 1))/(gamma6_1 - 1) ));
    P_tot6_2 = (P_tot_star/M6_2)*( ((1 + 0.5*(gamma6_2 - 1)*M6_2^2)/(0.5*(gamma6_2 + 1)))^( (0.5*(gamma6_2 + 1))/(gamma6_2 - 1) ));

    rho_star = rho6_1/((1/M6_1)*sqrt( 2*(1 + (gamma6_1 - 1)*0.5*M6_1^2)/(gamma6_1 + 1)));
    rho6_2 = rho_star*((1/M6_2)*sqrt( 2*(1 + (gamma6_2 - 1)*0.5*M6_2^2)/(gamma6_2 + 1)));

    c6_2 = sqrt(gamma6_2*R*T6_2);
    v6_2 = M6_2*c6_2;

    gamma6_2_new = gamma_N2(find(T==round(T6_2)),find(abs(P - round(P6_2,1)) < 0.001));

    err = abs(gamma6_2 - gamma6_2_new);

    gamma6_2 = gamma6_2_new;

end

clear gamma6_2_new

z = @(x) A_int/A_throat_int - (M6_2/x)*sqrt( ((1 + 0.5*(gamma6_2 - 1)*x^2)/(1 + 0.5*(gamma6_2 - 1)*M6_2^2))^((gamma6_2 + 1)/(gamma6_2 - 1)) );
M7 = fsolve(z,0.8);

P7 = P_tot6_2/(1 + ((gamma6_2 - 1)/2)*M7^2)^(gamma6_2/(gamma6_2 - 1));
T7 = T_tot/(1 + ((gamma6_2 - 1)/2)*M7^2);

%% After Venturi channel (point 7) and before servo valve (point 8)

T = (floor(T7)-12):0.5:(ceil(T7));
P = (floor(P7)-12):0.1:(ceil(P7));
data = nistdata('N2',T,P);

rho_N2 = data.Rho*data.Mw;           % Density of Nitrogen [kg/m^3] 
cp_N2 = data.Cp/data.Mw;             % Specific heat at constant pressure of Nitrogen [J/kgK]
cv_N2 = data.Cv/data.Mw;             % Specific heat at constant volume of Nitrogen [J/kgK]
gamma_N2 = cp_N2./cv_N2;             % Ratio of specific heats [-]
mu_N2 = data.mu;                     % Viscosity of Nitrogen [Pa*s]


rho7 = rho_N2(find(T==round(T7)),find(abs(P - round(P7,1)) < 0.001));     % Density downstream the mass flow meter (needle valve approx) [kg/m^3]
gamma7 = gamma_N2(find(T==round(T7)),find(abs(P - round(P7,1)) < 0.001)); % Ratio of specific heats the mass flow meter (needle valve approx)  [-]
gamma8 = gamma_N2(find(T==round(T7)),find(abs(P - round(P7,1)) < 0.001));
mu7 = mu_N2(find(T==round(T7)),find(abs(P - round(P7,1)) < 0.001));      % Viscosity downstream the mass flow meter (needle valve approx) [Pa*s]
v7 = m_dot_N2/(A_int*rho7);                       % Gas velocity downstream the mass flow meter (needle valve approx) [m/s]
c7 = (gamma7*R*T7)^0.5;                           % Sound speed downstream the mass flow meter (needle valve approx) [m/s]
M7 = v7/c7;                                       % Mach number downstream the mass flow meter (needle valve approx) [-]
Re7 = (rho7*v7*d_p_int)/mu7;                      % Reynolds number downstream the mass flow meter (needle valve approx) [-]

%% Venturi channel plots
x1 = 0;
x2 = 2.4;
x3 = x2 + 8;
x4 = x3 + 3.4;
d_vect = [x1 x2 x3 x4];
v_vect = [v6 v6_1 v6_2 v7];
P_vect = [P6 P6_1 P6_2 P7];
T_vect = [T6 T6_1 T6_2 T7];
rho_vect = [rho6 rho6_1 rho6_2 rho7];
M_vect = [M6 M6_1 M6_2 M7];

figure()
plot(d_vect,v_vect,'ro','linewidth',1.5)
grid on
xlabel('Position, $x_i$ $[cm]$')
ylabel('Velocity, $v_i$ $[m/s]$')
title('Velocity evolution - Venturi channel 12 mm diameter line')

figure()
plot(d_vect,M_vect,'ro','linewidth',1.5)
grid on
xlabel('Position, $x_i$ $[cm]$')
ylabel('Mach Number, $M_i$ $[m/s]$')
title('Mach evolution - Venturi channel 12 mm diameter line')

figure()
plot(d_vect,P_vect,'ro','linewidth',1.5)
grid on
xlabel('Position, $x_i$ $[cm]$')
ylabel('Pressure, $P_i$ $[bar]$')
title('Pressure evolution - Venturi channel 12 mm diameter line')


figure()
plot(d_vect,T_vect,'ro','linewidth',1.5)
grid on
xlabel('Position, $x_i$ $[cm]$')
ylabel('Temperature, $T_i$ $[K]$')
title('Temperature evolution - Venturi channel 12 mm diameter line')

figure()
plot(d_vect,rho_vect,'ro','linewidth',1.5)
grid on
xlabel('Position, $x_i$ $[cm]$')
ylabel('Density, $\rho_i$ $[kg/m^3]$')
title('Density evolution - Venturi channel 12 mm diameter line')

%%

L2 = 0.5;

if Re7 < 2300

        lambda = 64/Re7;

    else 

        z = @(x) 1/sqrt(x) + 2*log10(2.51/(Re7*sqrt(x)) + eps_rel/3.71);   % Colebrook-White correlation
        lambda = fsolve(z,0.0004);

end

iter = 0;
err = 1;

while err > 1e-3

    iter = iter + 1;
    
    g_M7 = (1 - M7^2)/(gamma7*M7^2) + ((gamma7 + 1)/(2*gamma7))*log(((gamma7 + 1)*M7^2)/(2 + (gamma7 - 1)*M7^2) );
    g_M8 = g_M7 - lambda*(L2/d_p_int);
    y = @(x) g_M8 - (1 - x^2)/(gamma8*x^2) + ((gamma8 + 1)/(2*gamma8))*log(((gamma8 + 1)*x^2)/(2 + (gamma8 - 1)*x^2) );
    M8 = fsolve(y,0.006);

    T_star = T7/(0.5*(gamma7 + 1)/(1 + (gamma7 - 1)*0.5*M7^2));
    T8 = T_star*(0.5*(gamma8 + 1)/(1 + (gamma8 - 1)*0.5*M8^2));
    
    P_star = P7/((1/M7)*sqrt(0.5*(gamma7 + 1)/(1 + (gamma7 - 1)*0.5*M7^2)));
    P8 = P_star*((1/M8)*sqrt(0.5*(gamma8 + 1)/(1 + (gamma8 - 1)*0.5*M8^2)));
    
    rho_star = rho7/((1/M7)*sqrt( 2*(1 + (gamma7 - 1)*0.5*M7^2)/(gamma7 + 1)));
    rho8 = rho_star*((1/M8)*sqrt( 2*(1 + (gamma8 - 1)*0.5*M8^2)/(gamma8 + 1)));
    
    c8 = sqrt(gamma8*R*T8);
    v8 = c8*M8;

    gamma8_new = gamma_N2(find(T==round(T8)),find(abs(P - round(P8,1)) < 0.001));

    err = abs(gamma8 - gamma8_new);

    gamma8 = gamma8_new;

end

clear gamma8_new
%% Before servovalve (point 8) and after servovalve (point 9)
G_g = rho8/1000;                      % Nitrogen specific gravity [-]
q_N2 = (m_dot_N2/rho8)*1000;           % Nitrogen volumetric flow rate [L/s]
% q_N2_std = (P8*q_N2*T8*60)/(1*273.15); % Nitrogen volumetric flow rate at std conditions [std L/min]
C_V = 3.8;                             % Flow coefficient ball valve

P9 = P8 - (G_g*(q_N2*60)^2)/(14.42*C_V)^2;                                                                   % Pressure downstream the servovalve (ball valve approx) [bar]
T9 = T8;                                                                              % Temperature downstream the servovalve (ball valve approx) [K]

%% Before check valve (point 9) and after check valve (point 10)

T = (floor(T9)-3):0.5:(ceil(T9));
P = (floor(P9)-2):0.1:(ceil(P9));
data = nistdata('N2',T,P);

rho_N2 = data.Rho*data.Mw;           % Density of Nitrogen [kg/m^3] 
cp_N2 = data.Cp/data.Mw;             % Specific heat at constant pressure of Nitrogen [J/kgK]
cv_N2 = data.Cv/data.Mw;             % Specific heat at constant volume of Nitrogen [J/kgK]
gamma_N2 = cp_N2./cv_N2;             % Ratio of specific heats [-]
mu_N2 = data.mu;  

rho9 = rho_N2(find(T==round(T9)),find(abs(P - round(P9,1)) < 0.001));     % Density downstream the servovalve (ball valve approx) [kg/m^3]
gamma9 = gamma_N2(find(T==round(T9)),find(abs(P - round(P9,1)) < 0.001)); % Ratio of specific heats the servovalve (ball valve approx)  [-]
mu9 = mu_N2(find(T==round(T9)),find(abs(P - round(P9,1)) < 0.001));      % Viscosity downstream the servovalve (ball valve approx) [Pa*s]
v9 = m_dot_N2/(A_int*rho9);                       % Gas velocity downstream the servovalve (ball valve approx) [m/s]
c9 = (gamma9*R*T9)^0.5;                           % Sound speed downstream the servovalve (ball valve approx) [m/s]
M9 = v9/c9;                                       % Mach number downstream the servovalve (ball valve approx) [-]
Re9 = (rho9*v9*d_p_int)/mu9;                      % Reynolds number downstream the servovalve (ball valve approx) [-]

G_g = rho9/1000;                      % Nitrogen specific gravity [-]
q_N2 = (m_dot_N2/rho9)*1000;           % Nitrogen volumetric flow rate [L/s]
%q_N2_std = (P9*q_N2*T9*60)/(1*273.15); % Nitrogen volumetric flow rate at std conditions [std L/min]
C_V = 1.68;                             % Flow coefficient check valve

% z = @(x) q_N2_std - 6950*C_V*P9*(1 - (2/3)*(P9 - x)/P9)*sqrt((P9 - x)/(P9*G_g*T9));
P10 = P9 - (G_g*(q_N2*60)^2)/(14.42*C_V)^2;                                                                % Pressure downstream the check valve [bar]
T10 = T9;

%% After check valve (point 10) and before injector (point 11)

T = (floor(T10)-100):0.5:(ceil(T10));
P = (floor(P10)-10):0.1:(ceil(P10));
data = nistdata('N2',T,P);

rho_N2 = data.Rho*data.Mw;           % Density of Nitrogen [kg/m^3] 
cp_N2 = data.Cp/data.Mw;             % Specific heat at constant pressure of Nitrogen [J/kgK]
cv_N2 = data.Cv/data.Mw;             % Specific heat at constant volume of Nitrogen [J/kgK]
gamma_N2 = cp_N2./cv_N2;             % Ratio of specific heats [-]
mu_N2 = data.mu;  

rho10 = rho_N2(find(T==round(T10)),find(abs(P - round(P10,1)) < 0.001));     % Density downstream the check valve [kg/m^3]
gamma10 = gamma_N2(find(T==round(T10)),find(abs(P - round(P10,1)) < 0.001)); % Ratio of specific heats the check valve  [-]
gamma11 = gamma_N2(find(T==round(T10)),find(abs(P - round(P10,1)) < 0.001));
mu10 = mu_N2(find(T==round(T10)),find(abs(P - round(P10,1)) < 0.001));      % Viscosity downstream the check valve [Pa*s]
v10 = m_dot_N2/(A_int*rho10);                      % Gas velocity downstream the check valve [m/s]
c10 = (gamma10*R*T10)^0.5;                         % Sound speed downstream the check valve [m/s]
M10 = v10/c10;                                     % Mach number downstream the check valve [-]
Re10 = (rho10*v10*d_p_int)/mu10;                   % Reynolds number downstream the check valve [-]

if Re10 < 2300

        lambda = 64/Re10;

    else 

        z = @(x) 1/sqrt(x) + 2*log10(2.51/(Re10*sqrt(x)) + eps_rel/3.71);   % Colebrook-White correlation
        lambda = fsolve(z,0.0004);

end

iter = 0;
err = 1;
L_fin = 0.4;

while err > 1e-3

    iter = iter + 1;
    g_M10 = (1 - M10^2)/(gamma10*M10^2) + ((gamma10 + 1)/(2*gamma10))*log(((gamma10 + 1)*M10^2)/(2 + (gamma10 - 1)*M10^2) );
    g_M11 = g_M10 - lambda*(L_fin/d_p_int);

    if g_M11 < 0

        M11 = 1;
    else

        y = @(x) g_M11 - (1 - x^2)/(gamma11*x^2) + ((gamma11 + 1)/(2*gamma11))*log(((gamma11 + 1)*x^2)/(2 + (gamma11 - 1)*x^2) );
        M11 = fsolve(y,0.006);

    end

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

%% Injector pressure loss
P12 = 1;
delta_P_inj = (P11 - P12)*1e5;

N_inj = 1;
C_d = 0.61;                              % Sharp-edged orifice with diameter bigger than 2.5 mm
A_needed = m_dot_N2/(C_d*sqrt(2*delta_P_inj*rho11));
A_inj = A_needed./N_inj;
d_inj = sqrt((4*A_inj)/pi);
v_inj=C_d*sqrt(2*delta_P_inj/rho11);
A_slab = 30*30*10e-6;                    % Area of slab test facility [m^2]

%% Venturi tube's size

d_conv = d_p_int;
d_div = d_p_int;
alpha_conv = 7;
alpha_div = 5;
L_conv = (d_p_int - d_throat_int)/(2*tand(alpha_conv));
L_div = (d_p_int - d_throat_int)/(2*tand(alpha_div));

%% Figures


P_vect = [P_reg P1 P2 P3 P4 P5 P6 P7 P8 P9 P10 P11 P12];
x_reg = 0;
x1 = 0;
x2 = x1 + L1;
x3 = x1 + L1;
x4 = x3 + L1;
x5 = x3 + L1;
x6 = x5 + L1;
% x6_1 = x6 + L_conv;
% x6_2 = x6_1 + L_throat/1000;
x7 = x6 + L_conv + L_throat/1000 + L_div;
x8 = x7 + L2;
x9 = x7 + L2;
x10 = x7 + L2;
x11 = x10 + L_fin;
x = [x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11];

figure()
plot(x,P_vect,'ro','linewidth',1.5)
grid on
xlabel('Position on the line, $L_i$ [m]')
ylabel('Pressure, $P_i$ [bar]')
title('Pressure evolution: 12 mm diameter pipes')

figure()
plot(N_inj,d_inj*1e3,'ro','linewidth',1.5)
grid on
xlabel('Number of injectors, $N_{inj}$ [-]')
ylabel('Injector diameter, $d_{inj}$ [mm]')
title('Number of injectors vs diameter of injectors')