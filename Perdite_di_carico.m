clc
clear


set(0,'DefaultTextFontSize',12);              % Settings for the plot
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultLegendInterpreter','Latex');
set(groot,'DefaultAxesTickLabelInterpreter','Latex');
set(0,'DefaultTextInterpreter','Latex');
set(0,'DefaultLegendFontSize',12);

d_p_ext = 12*1e-3;                   % Pipe external diameter [m]
t = 0.9*1e-3;                        % Thickness of the tube  [m]
d_p_int = d_p_ext - 2*t;             % Pipe internal diameter [m]
A_int = pi*(d_p_int/2)^2;            % Internal cross sectional area [m^2]

eps = 0.015*1e-3;                    % Absolute roughness of stainless steel [m]
eps_rel = eps/d_p_ext;               % Relative roughness of stainless steel [-]

%%
T = 298.15:1:300.15;
P = 33:0.05:35;
data = nistdata('N2',T,P);

rho_N2 = data.Rho*data.Mw;           % Density of Nitrogen [kg/m^3] 
cp_N2 = data.Cp/data.Mw;             % Specific heat at constant pressure of Nitrogen [J/kgK]
cv_N2 = data.Cv/data.Mw;             % Specific heat at constant volume of Nitrogen [J/kgK]
gamma_N2 = cp_N2./cv_N2;             % Ratio of specific heats [-]
mu_N2 = data.mu;                     % Viscosity of Nitrogen [Pa*s]
% L_test = 30*1e-3;                    % Dimension of the square test chamber [m]
% rho = 1.832;                         % Density in full scale facility [kg/m^3]
% u = 22.12;                           % Velocity in full scale facility [m/s]
% mom_flux = rho*u^2;                  % Momentum flux in full scale facility [kg/ms^2]

m_dot_N2 = 60*1e-3;                    % Nitrogen mass flow rate [kg/s]
R = 8314/28;                           % Specific ideal gas constant [J/kgK]

% P_c = 7 
% delta_P_inj = 0.4*100*sqrt(10*P_c*1e5); % Pressure drop across the injection plate [Pa] 

L = 1;                                   % Length of the tube [m]

% d_inj = 0.5*1e-3;                    % Injector diameter (conical entrance) [m]
% C_d = 0.7;                           % Discharge coefficient
% A_inj = pi*0.25*d_inj^2;             % Injector cross sectional area

%for i = 1:length(lambda)
%f= @(x) ((x^2-(P_c*10^5+delta_P_inj)^2)/(2*x))-lambda(i)*(x*28/(8314*T))*(L/d_p_int)*(m_dot_N2/(A_int*(x*28/(8314*T))))^2*0.5;
%P_distr(i) = fzero(f,0.1);
%end


T1 = T(1);                                      % Temperature downstream the pressure regulator [K]
P1 = 35;                                        % Pressure downstream the pressure regulator [bar]
rho1 = rho_N2(find(T==298.15),find(P==35));     % Density downstream the pressure regulator [kg/m^3]
gamma1 = gamma_N2(find(T==298.15),find(P==35)); % Ratio of specific heats downstream the pressure regulator [-]
mu1 = mu_N2(find(T==298.15),find(P==35));       % Viscosity downstream the pressure regulator [Pa*s]
v1 = m_dot_N2/(A_int*rho1);                     % Gas velocity downstream the pressure regulator [m/s]
c1 = (gamma1*R*T1)^0.5;                         % Sound speed downstream the pressure regulator [m/s]
M1 = v1/c1;                                     % Mach number downstream the pressure regulator [-]
Re1 = (rho1*v1*d_p_int)/mu1;                    % Reynolds number downstream the pressure regulator [-]

if Re1 < 2300

        lambda = 64/Re1;

    else 

        z = @(x) 1/sqrt(x) + 2*log10(2.51/(Re1*sqrt(x)) + eps_rel/3.71);   % Colebrook-White correlation
        lambda = fsolve(z,0.0004);

end

g_M1 = (1 - M1^2)/(gamma1*M1^2) + ((gamma1 + 1)/(2*gamma1))*log(((gamma1 + 1)*M1^2)/(2 + (gamma1 - 1)*M1^2) );
g_M2 = g_M1 - lambda*(L/d_p_int);

y = @(x) g_M2 - (1 - x^2)/(gamma1*x^2) + ((gamma1 + 1)/(2*gamma1))*log(((gamma1 + 1)*x^2)/(2 + (gamma1 - 1)*x^2) );
M2 = fsolve(y,0.006);

T_star = T1/(0.5*(gamma1 + 1)/(1 + (gamma1 - 1)*0.5*M1^2));
T2 = T_star*(0.5*(gamma1 + 1)/(1 + (gamma1 - 1)*0.5*M2^2));

P_star = P1/((1/M1)*sqrt(0.5*(gamma1 + 1)/(1 + (gamma1 - 1)*0.5*M1^2)));
P2 = P_star*((1/M2)*sqrt(0.5*(gamma1 + 1)/(1 + (gamma1 - 1)*0.5*M2^2)));

rho_star = rho1/((1/M1)*sqrt( 2*(1 + (gamma1 - 1)*0.5*M1^2)/(gamma1 + 1)));
rho2 = rho_star*((1/M2)*sqrt( 2*(1 + (gamma1 - 1)*0.5*M2^2)/(gamma1 + 1)));

% gamma in un tratto di tubo è al momento considerata costante, anche se in
% realtà varia siccome variano pressione e temperatura, siccome una delle
% ipotesi del flusso di Fanno è che gamma sia costante. Si potrebbe fare un
% ciclo for usando come guess iniziale gamma1 per il punto 2, per poi 
% aggiornarla iterativamente finchè non si giunge a convergenza. Al
% momento, considerando P1 = 35 bar, sembra non ne valga la pena: la
% pressione nel punto 2 cambia veramente poco.

%%
c2 = sqrt(gamma1*R*T2);
v2 = M2*c2;

P3 = 1e-5*(P2*1e5 - 0.1*rho2*2^2);                % Pressure drop related to the T fitting (dump line)  [bar]
T3 = T2;

rho3 = rho_N2(find(T==298.15),find(P==34.2));     % Density downstream the pressure regulator [kg/m^3]
gamma3 = gamma_N2(find(T==298.15),find(P==34.2)); % Ratio of specific heats downstream the pressure regulator [-]
gamma4 = gamma_N2(find(T==298.15),find(P==33.35));
mu3 = mu_N2(find(T==298.15),find(P==34.2));       % Viscosity downstream the pressure regulator [Pa*s]
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

g_M3 = (1 - M3^2)/(gamma3*M3^2) + ((gamma3 + 1)/(2*gamma3))*log(((gamma3 + 1)*M3^2)/(2 + (gamma3 - 1)*M3^2) );
g_M4 = g_M3 - lambda*(L/d_p_int);

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
%%
G_g = rho4/1.205;                      % Nitrogen specific gravity [-]
q_N2 = (m_dot_N2/rho4)*1000;           % Nitrogen volumetric flow rate [L/s]
q_N2_std = (P4*q_N2*T4*60)/(1*273.15); % Nitrogen volumetric flow rate at std conditions [std L/min]
C_V = 5.2;

z = @(x) q_N2_std - 6950*C_V*P4*(1 - (2/3)*(P4 - x)/P4)*sqrt((P4 - x)/(P4*G_g*T4));
P5 = fsolve(z,32);
T5 = T4;

%% 
T = 297.45:0.5:298.45;
P = 29.5:0.05:30.5;
data = nistdata('N2',T,P);

rho_N2 = data.Rho*data.Mw;           % Density of Nitrogen [kg/m^3] 
cp_N2 = data.Cp/data.Mw;             % Specific heat at constant pressure of Nitrogen [J/kgK]
cv_N2 = data.Cv/data.Mw;             % Specific heat at constant volume of Nitrogen [J/kgK]
gamma_N2 = cp_N2./cv_N2;             % Ratio of specific heats [-]
mu_N2 = data.mu;  

rho5 = rho_N2(find(T==297.95),find(P==30.45));     % Density downstream the manual ball valve [kg/m^3]
gamma5 = gamma_N2(find(T==297.95),find(P==30.45)); % Ratio of specific heats downstream the manual ball valve  [-]
gamma6 = gamma_N2(find(T==297.95),find(P==29.55));
mu5 = mu_N2(find(T==297.95),find(P==30.45));       % Viscosity downstream the manual ball valve [Pa*s]
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

g_M5 = (1 - M5^2)/(gamma5*M5^2) + ((gamma5 + 1)/(2*gamma5))*log(((gamma5 + 1)*M5^2)/(2 + (gamma5 - 1)*M5^2) );
g_M6 = g_M5 - lambda*(L/d_p_int);

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

%%
G_g = rho6/1.205;                      % Nitrogen specific gravity [-]
q_N2 = (m_dot_N2/rho6)*1000;           % Nitrogen volumetric flow rate [L/s]
q_N2_std = (P6*q_N2*T6*60)/(1*273.15); % Nitrogen volumetric flow rate at std conditions [std L/min]
C_V = 1.8;                             % Flow coefficient ball valve

z = @(x) q_N2_std - 6950*C_V*P6*(1 - (2/3)*(P6 - x)/P6)*sqrt((P6 - x)/(P6*G_g*T6));
P7 = fsolve(z,20);                                                                    % Pressure downstream the mass flow meter (needle valve approx) [bar]
T7 = T6;                                                                              % Temperature downstream the mass flow meter (needle valve approx) [K]

%%
T = 297.25:0.5:298.75;
P = 13.3:0.1:14.8;
data = nistdata('N2',T,P);

rho_N2 = data.Rho*data.Mw;           % Density of Nitrogen [kg/m^3] 
cp_N2 = data.Cp/data.Mw;             % Specific heat at constant pressure of Nitrogen [J/kgK]
cv_N2 = data.Cv/data.Mw;             % Specific heat at constant volume of Nitrogen [J/kgK]
gamma_N2 = cp_N2./cv_N2;             % Ratio of specific heats [-]
mu_N2 = data.mu;  

rho7 = rho_N2(find(T==297.75),find(P==14.8));     % Density downstream the mass flow meter (needle valve approx) [kg/m^3]
gamma7 = gamma_N2(find(T==297.75),find(P==14.8)); % Ratio of specific heats the mass flow meter (needle valve approx)  [-]
gamma8 = gamma_N2(find(T==297.25),find(P==13.3));
mu7 = mu_N2(find(T==297.75),find(P==14.8));       % Viscosity downstream the mass flow meter (needle valve approx) [Pa*s]
v7 = m_dot_N2/(A_int*rho7);                       % Gas velocity downstream the mass flow meter (needle valve approx) [m/s]
c7 = (gamma7*R*T7)^0.5;                           % Sound speed downstream the mass flow meter (needle valve approx) [m/s]
M7 = v7/c7;                                       % Mach number downstream the mass flow meter (needle valve approx) [-]
Re7 = (rho7*v7*d_p_int)/mu7;                      % Reynolds number downstream the mass flow meter (needle valve approx) [-]

if Re7 < 2300

        lambda = 64/Re7;

    else 

        z = @(x) 1/sqrt(x) + 2*log10(2.51/(Re7*sqrt(x)) + eps_rel/3.71);   % Colebrook-White correlation
        lambda = fsolve(z,0.0004);

end

g_M7 = (1 - M7^2)/(gamma7*M7^2) + ((gamma7 + 1)/(2*gamma7))*log(((gamma7 + 1)*M7^2)/(2 + (gamma7 - 1)*M7^2) );
L = 1;
g_M8 = g_M7 - lambda*(L/d_p_int);

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

%%
G_g = rho8/1.205;                      % Nitrogen specific gravity [-]
q_N2 = (m_dot_N2/rho8)*1000;           % Nitrogen volumetric flow rate [L/s]
q_N2_std = (P8*q_N2*T8*60)/(1*273.15); % Nitrogen volumetric flow rate at std conditions [std L/min]
C_V = 5.2;                             % Flow coefficient ball valve

z = @(x) q_N2_std - 6950*C_V*P8*(1 - (2/3)*(P8 - x)/P8)*sqrt((P8 - x)/(P8*G_g*T8));
P9 = fsolve(z,13);                                                                    % Pressure downstream the servovalve (ball valve approx) [bar]
T9 = T8;                                                                              % Temperature downstream the servovalve (ball valve approx) [K]

%%
T = 297.25:0.5:298.75;
P = 9.35:0.1:9.45;
data = nistdata('N2',T,P);

rho_N2 = data.Rho*data.Mw;           % Density of Nitrogen [kg/m^3] 
cp_N2 = data.Cp/data.Mw;             % Specific heat at constant pressure of Nitrogen [J/kgK]
cv_N2 = data.Cv/data.Mw;             % Specific heat at constant volume of Nitrogen [J/kgK]
gamma_N2 = cp_N2./cv_N2;             % Ratio of specific heats [-]
mu_N2 = data.mu;  

rho9 = rho_N2(find(T==297.25),find(P==9.35));     % Density downstream the servovalve (ball valve approx) [kg/m^3]
gamma9 = gamma_N2(find(T==297.25),find(P==9.35)); % Ratio of specific heats the servovalve (ball valve approx)  [-]
mu9 = mu_N2(find(T==297.25),find(P==9.35));       % Viscosity downstream the servovalve (ball valve approx) [Pa*s]
v9 = m_dot_N2/(A_int*rho9);                       % Gas velocity downstream the servovalve (ball valve approx) [m/s]
c9 = (gamma9*R*T9)^0.5;                           % Sound speed downstream the servovalve (ball valve approx) [m/s]
M9 = v9/c9;                                       % Mach number downstream the servovalve (ball valve approx) [-]
Re9 = (rho9*v9*d_p_int)/mu9;                      % Reynolds number downstream the servovalve (ball valve approx) [-]

G_g = rho9/1.205;                      % Nitrogen specific gravity [-]
q_N2 = (m_dot_N2/rho9)*1000;           % Nitrogen volumetric flow rate [L/s]
q_N2_std = (P9*q_N2*T9*60)/(1*273.15); % Nitrogen volumetric flow rate at std conditions [std L/min]
C_V = 1.8;                             % Flow coefficient ball valve

z = @(x) q_N2_std - 6950*C_V*P9*(1 - (2/3)*(P9 - x)/P9)*sqrt((P9 - x)/(P9*G_g*T9));
P10 = fsolve(z,9);                                                                    % Pressure downstream the check valve [bar]
T10 = T9;

%%
T = 297.25:0.5:298.75;
P = 3.7:0.1:4.8;
data = nistdata('N2',T,P);

rho_N2 = data.Rho*data.Mw;           % Density of Nitrogen [kg/m^3] 
cp_N2 = data.Cp/data.Mw;             % Specific heat at constant pressure of Nitrogen [J/kgK]
cv_N2 = data.Cv/data.Mw;             % Specific heat at constant volume of Nitrogen [J/kgK]
gamma_N2 = cp_N2./cv_N2;             % Ratio of specific heats [-]
mu_N2 = data.mu;  

rho10 = rho_N2(find(T==297.25),find(P==4.7));     % Density downstream the check valve [kg/m^3]
gamma10 = gamma_N2(find(T==297.25),find(P==4.7)); % Ratio of specific heats the check valve  [-]
gamma11 = gamma_N2(find(T==297.25),find(P==3.7));
mu10 = mu_N2(find(T==297.25),find(P==4.7));       % Viscosity downstream the check valve [Pa*s]
v10 = m_dot_N2/(A_int*rho9);                       % Gas velocity downstream the check valve [m/s]
c10 = (gamma9*R*T9)^0.5;                           % Sound speed downstream the check valve [m/s]
M10 = v10/c10;                                     % Mach number downstream the check valve [-]
Re10 = (rho10*v10*d_p_int)/mu10;                   % Reynolds number downstream the check valve [-]

if Re10 < 2300

        lambda = 64/Re10;

    else 

        z = @(x) 1/sqrt(x) + 2*log10(2.51/(Re10*sqrt(x)) + eps_rel/3.71);   % Colebrook-White correlation
        lambda = fsolve(z,0.0004);

end

g_M10 = (1 - M10^2)/(gamma10*M10^2) + ((gamma10 + 1)/(2*gamma10))*log(((gamma10 + 1)*M10^2)/(2 + (gamma10 - 1)*M10^2) );
L = 1;
g_M11 = g_M10 - lambda*(L/d_p_int);

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