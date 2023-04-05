clc
clear


set(0,'DefaultTextFontSize',12);              % Settings for the plot
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultLegendInterpreter','Latex');
set(groot,'DefaultAxesTickLabelInterpreter','Latex');
set(0,'DefaultTextInterpreter','Latex');
set(0,'DefaultLegendFontSize',12);
T = 298.15:1:350.15;
P = 7:1:50;
data = nistdata('N2',T,P);
%% 
P_c = 7;                             % Pressure in the test chamber [bar]
d_p_ext = 12*1e-3;                   % Pipe external diameter [m]
t = 0.9*1e-3;                        % Thickness of the tube  [m]
d_p_int = d_p_ext - 2*t;             % Pipe internal diameter [m]
A_int = pi*(d_p_int/2)^2;            % Internal cross sectional area [m^2]

eps = 0.015*1e-3;                    % Absolute roughness of stainless steel [m]
eps_rel = eps/d_p_ext;               % Relative roughness of stainless steel [-]

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

delta_P_inj = 0.4*100*sqrt(10*P_c*1e5); % Pressure drop across the injection plate [Pa] 

L=1;                                   % Length of the tube [m]

% d_inj = 0.5*1e-3;                    % Injector diameter (conical entrance) [m]
% C_d = 0.7;                           % Discharge coefficient
% A_inj = pi*0.25*d_inj^2;             % Injector cross sectional area

% c = (gamma*(8314/28)*T).^0.5;        % Speed of sound [m/s]
% M = v_N2./c;                         % Mach number [-]
% 
% Re = (rho_N2.*v_N2*d_p_int)./mu_N2;  % Reynolds number [-]
% 
% 
% lambda = zeros(length(Re),1);
% 
% for i = 1:length(Re)
% 
%     if Re(i) < 2300
% 
%         lambda(i) = 64/Re(i);
% 
%     else 
% 
%         z = @(x) 1/sqrt(x) + 2*log10(2.51/(Re(i)*sqrt(x)) + eps_rel/3.71);   % Colebrook-White correlation
%         lambda(i) = fsolve(z,0.0004);
% 
%     end
% 
% end

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

% CADUTE DI PRESSIONE MOLTO BASSE SE TUBO LUNGO 1 METRO: PENSARE A
% DISCRETIZZARE CON UNO STEP PIU' BASSO LE PRESSIONI PER I DATI NEL NIST.
% OPPURE FORSE CONVIENE RICALCOLARE LE PROPRIETA' DEL NIST DI VOLTA IN
% VOLTA IN BASE AI RISULTATI CHE CI ESCONO, PIUTTOSTO CHE UTILIZZARE SEMPRE
% I DATI CALCOLATI ALL'INIZIO. IN QUESTO MODO POTREMMO DISCRETIZZARE DI
% VOLTA IN VOLTA UTILIZZANDO UNO STEP DIVERSO IN BASE ALLE ESIGENZE.