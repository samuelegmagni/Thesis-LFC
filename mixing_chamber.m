%% Mixing chamber sizing
clear
clc


% Define mass flow rates
mdot_N2 = 60*1e-3;
mdot_SRP = 17.8*1e-3;  

% Define pressures and temperatures
T = [280:2:600];
P = [1:0.2:13];
data = nistdata('N2',T,P);
Tmix = 280;                             
Ptank = 12;                                   % pressure in tank [bar]
Pburst = 4*Ptank*1e5;                         % [Pa]

% Define densities
rho_N2 = data.Rho*data.Mw; 
rho_N2_mix= rho_N2(find(T==round(Tmix)),find(abs(P - round(Ptank,1)) < 0.001)); 
rho_SRP_gas = 3.5928;                      % density from NASA CEA [kg/m^3]

% Define volumetric mass flow rate
qvol_N2 = mdot_N2/rho_N2_mix;
qvol_SRP = mdot_SRP/rho_SRP_gas;

% Define volume of mixing chamber in 1 s of flow
t_res=[ 0.02 0.025 0.03 0.035 0.04]
V_N2 = qvol_N2.*t_res;
V_SRP = qvol_SRP.*t_res;
Vtank = 1.05*(V_N2+V_SRP);

% Size mixing chamber
rho_steel = 8000;
sigma_ult = 505*1e6;
sigma_snerv = 215*1e6;
htank = 5*1e-2;
rtank = sqrt(Vtank./(pi*htank));
t_des = Pburst*rtank./sigma_snerv;
t_ASME = Pburst*rtank./(0.8*sigma_ult-0.6*Pburst);
Vtank_int = pi*(rtank-t_des).^2*htank;
mtank = (Vtank-Vtank_int).*rho_steel;

% Von Mises stresses
sigma_l = Pburst.*rtank/(2*t_des);
sigma_h = Pburst.*rtank/t_des;
sigma_vm = (1/sqrt(2))*sqrt((sigma_h-sigma_l)^2+sigma_l^2+sigma_h^2);

%Plot 
figure()
plot(t_res,rtank*1e2,'ro','linewidth',1.5)
grid on
xlabel('Residence time [s]')
ylabel('Tank radius [cm]')
title('Tank radius vs residence time in chamber')