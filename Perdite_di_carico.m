clc
clear


set(0,'DefaultTextFontSize',12);              % Settings for the plot
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultLegendInterpreter','Latex');
set(groot,'DefaultAxesTickLabelInterpreter','Latex');
set(0,'DefaultTextInterpreter','Latex');
set(0,'DefaultLegendFontSize',12);
T = 298.15;
P = 7:1:20;
data = nistdata('N2',T,P);
%% 
P_c = 7;                             % Pressure in the test chamber [bar]
d_p_ext = 12*1e-3;                   % Pipe external diameter [m]
t = 0.9*1e-3;                        % Thickness of the tube  [m]
d_p_int = d_p_ext - 2*t;             % Pipe internal diameter [m]
A_int = pi*(d_p_int/2)^2;            % Internal cross sectional area [m^2]
T_amb = 298.15;                      % Ambient temperature [K]
rho_N2 = data.Rho*data.Mw;           % Density of Nitrogen [kg/m^3] 
cp_N2 = data.Cp/data.Mw;             % Specific heat at constant pressure of Nitrogen [J/kgK]
cv_N2 = data.Cv/data.Mw;             % Specific heat at constant volume of Nitrogen [J/kgK]
gamma = cp_N2./cv_N2;                % Ratio of specific heats [-]
mu_N2 = data.mu;                     % Viscosity of Nitrogen [Pa*s]
% L_test = 30*1e-3;                    % Dimension of the square test chamber [m]
% rho = 1.832;                         % Density in full scale facility [kg/m^3]
% u = 22.12;                           % Velocity in full scale facility [m/s]
% mom_flux = rho*u^2;                  % Momentum flux in full scale facility [kg/ms^2]
m_dot_N2 = 60*1e-3;                  % Nitrogen mass flow rate [kg/s]
v_N2 = m_dot_N2./(A_int*rho_N2);     % Velocity of the fluid inside the pipes [m/s]
delta_P_inj = 1e-5*40*sqrt(10*P_c*1e5);  

d_inj = 0.5*1e-3;                    % Injector diameter (conical entrance) [m]
C_d = 0.7;                           % Discharge coefficient
A_inj = pi*0.25*d_inj^2;             % Injector cross sectional area

c = (gamma*(8314/28)*T).^0.5;        % Speed of sound [m/s]
M = v_N2./c;                         % Mach number [-]

Re = (rho_N2.*v_N2*d_p_int)./mu_N2;  % Reynolds number [-]

eps = 0.015*1e-3;                     % Absolute roughness of stainless steel [m]
eps_rel = eps/d_p_ext;                % Relative roughness of stainless steel [-]

lambda = zeros(length(Re),1);

for i = 1:length(Re)
    
    if Re(i) < 2300
        
        lambda(i) = 64/Re(i);
    
    else 
        
        z = @(x) 1/sqrt(x) + 2*log10(2.51/(Re(i)*sqrt(x)) + eps_rel/3.71);   % Colebrook-White correlation
        lambda(i) = fsolve(z,0.0004);

    end

end