clear
clc

format long

set(0,'DefaultTextFontSize',12);              % Settings for the plot
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultLegendInterpreter','Latex');
set(0,'DefaultAxesTickLabelInterpreter','Latex');
set(0,'DefaultTextInterpreter','Latex');
set(0,'DefaultLegendFontSize',12);

%%
T_vect = [298:5:850];
P_vect = [1:5:210];

data = nistdata('N2',T_vect,P_vect);

%% GE-C parameters

R_CC = 79.2 * 1e-3;
R_t = 14.42 * 1e-3;
eps_c = 30.16;
eps_e = 125;
L_CC = 22*1e-2;
m_dot = 0.798;
rho = 1.832;
T = 2787;
u = 22.12;
M_m = 21.23; %g/mol
c_p = 2519;  %J/kgK
mu = 9.76*1e-5;
lambda = 0.3372;

A_tot = pi*2*R_CC*L_CC + 2*pi*R_CC^2;

%% Full scale test section: Re similarity

D_test = 2*R_CC;

Re = rho*u*2*R_CC/mu;     %Re_paper = 65797

u_exp = Re*(data.mu./( data.Rho*data.Mw*D_test ));

m_dot_exp = data.Rho*data.Mw.*u_exp*0.25*pi*D_test^2;

adjust_vect = ones(length(T_vect),length(P_vect));

T_amb = 298;

q_dot_exp = m_dot_exp.*( (data.Cp/data.Mw).*T_vect' - ( (data.Cp(1,1)/data.Mw)*T_amb*adjust_vect ) );

figure()
contourf(data.P*1e-5,data.T,m_dot_exp); 
shading flat;
title('Mass flow rate $Re_D$ similarity')
xlabel('Pressure $[bar]$')
ylabel('Temperature $[K]$')
c = colorbar;
c.Label.String = 'm_{N_2} [kg/s]';

figure()
contourf(data.P*1e-5,data.T,1e-3*q_dot_exp/1);     % area di riferimento 1 m^2
shading flat;
title('Energy flux $Re_D$ similarity')
xlabel('Pressure $[bar]$')
ylabel('Temperature $[K]$')
c = colorbar;
c.Label.String = 'q_{N_2} [kW/m^2]';

%% Full scale test section: mass flux similarity

G = rho*u;      % G_paper = 40.54;

u_exp = G./(data.Rho*data.Mw);

m_dot_exp = data.Rho*data.Mw.*u_exp*0.25*pi*D_test^2;

q_dot_exp = m_dot_exp.*( (data.Cp/data.Mw).*T_vect' - ( (data.Cp(1,1)/data.Mw)*T_amb*adjust_vect ) );

figure()
contourf(data.P*1e-5,data.T,m_dot_exp); 
title('Mass flow rate: mass flux similarity')
xlabel('Pressure $[bar]$')
ylabel('Temperature $[K]$')
c = colorbar;
c.Label.String = 'm_{N_2} [kg/s]';

figure()
contourf(data.P*1e-5,data.T,1e-3*q_dot_exp/1);     % area di riferimento 1 m^2
title('Energy flux: mass flux similarity')
xlabel('Pressure $[bar]$')
ylabel('Temperature $[K]$')
c = colorbar;
c.Label.String = 'q_{N_2} [kW/m^2]';

%% Full scale test section: momentum flux similarity

mom_flux = rho*u^2;      % mom_flux_paper = 896.89;

u_exp = sqrt(mom_flux./(data.Rho*data.Mw));

m_dot_exp = data.Rho*data.Mw.*u_exp*0.25*pi*D_test^2;

q_dot_exp = m_dot_exp.*( (data.Cp/data.Mw).*T_vect' - ( (data.Cp(1,1)/data.Mw)*T_amb*adjust_vect ) );

figure()
contourf(data.P*1e-5,data.T,m_dot_exp); 
title('Mass flow rate: momentum flux similarity')
xlabel('Pressure $[bar]$')
ylabel('Temperature $[K]$')
c = colorbar;
c.Label.String = 'm_{N_2} [kg/s]';

figure()
contourf(data.P*1e-5,data.T,1e-3*q_dot_exp/1);     % area di riferimento 1 m^2
title('Energy flux: momentum flux similarity')
xlabel('Pressure $[bar]$')
ylabel('Temperature $[K]$')
c = colorbar;
c.Label.String = 'q_{N_2} [kW/m^2]';

%% Full scale test section: energy flux similarity

energy_flux = rho*u*c_p*T;      % energy_flux_paper = 284.6*1e6;

u_exp = energy_flux./( ((data.Rho*data.Mw).*(data.Cp/data.Mw)) .* (T_vect'*ones(1,length(P_vect))) );

m_dot_exp = data.Rho*data.Mw.*u_exp*0.25*pi*D_test^2;

q_dot_exp = m_dot_exp.*( (data.Cp/data.Mw).*T_vect' - ( (data.Cp(1,1)/data.Mw)*T_amb*adjust_vect ) );

figure()
contourf(data.P*1e-5,data.T,m_dot_exp); 
title('Mass flow rate: energy flux similarity')
xlabel('Pressure $[bar]$')
ylabel('Temperature $[K]$')
c = colorbar;
c.Label.String = 'm_{N_2} [kg/s]';
clim([6 16])

figure()
contourf(data.P*1e-5,data.T,1e-3*q_dot_exp/1);     % area di riferimento 1 m^2
title('Energy flux: energy flux similarity')
xlabel('Pressure $[bar]$')
ylabel('Temperature $[K]$')
c = colorbar;
c.Label.String = 'q_{N_2} [kW/m^2]';

%% Small scale test facility

T_small = [298:5:600];
P_small = [1:0.2:8];

data = nistdata('N2',T_small,P_small);

%% 30x30 momentum flux analogy
L_test = 30*1e-3;

mom_flux = rho*u^2;      % mom_flux_paper = 896.89;

u_exp = sqrt(mom_flux./(data.Rho*data.Mw));

m_dot_exp = data.Rho*data.Mw.*u_exp*L_test^2;

adjust_vect = ones(length(T_small),length(P_small));

q_dot_exp = m_dot_exp.*( (data.Cp/data.Mw).*T_small' - ( (data.Cp(1,1)/data.Mw)*T_small(1)*adjust_vect ) );

figure()
contourf(data.P*1e-5,data.T,m_dot_exp); 
title('Slab 30x30: mass flow rate momentum flux analogy')
xlabel('Pressure $[bar]$')
ylabel('Temperature $[K]$')
c = colorbar;
c.Label.String = 'm_{N_2} [kg/s]';

figure()
contourf(data.P*1e-5,data.T,1e-3*q_dot_exp/1);     % area di riferimento 1 m^2
title('Slab 30x30: energy flux momentum flux analogy')
xlabel('Pressure $[bar]$')
ylabel('Temperature $[K]$')
c = colorbar;
c.Label.String = 'q_{N_2} [kW/m^2]';

%% 49x30 momentum flux analogy
A_test = 49*30*1e-6;

mom_flux = rho*u^2;      % mom_flux_paper = 896.89;

u_exp = sqrt(mom_flux./(data.Rho*data.Mw));

m_dot_exp = data.Rho*data.Mw.*u_exp*A_test;

q_dot_exp = m_dot_exp.*( (data.Cp/data.Mw).*T_small' - ( (data.Cp(1,1)/data.Mw)*T_small(1)*adjust_vect ) );

figure()
contourf(data.P*1e-5,data.T,m_dot_exp); 
title('Slab 49x30: mass flow rate momentum flux analogy')
xlabel('Pressure $[bar]$')
ylabel('Temperature $[K]$')
c = colorbar;
c.Label.String = 'm_{N_2} [kg/s]';

figure()
contourf(data.P*1e-5,data.T,1e-3*q_dot_exp/1);     % area di riferimento 1 m^2
title('Slab 49x30: energy flux momentum flux analogy')
xlabel('Pressure $[bar]$')
ylabel('Temperature $[K]$')
c = colorbar;
c.Label.String = 'q_{N_2} [kW/m^2]';

%% 69x30 momentum flux analogy
A_test = 69*30*1e-6;

mom_flux = rho*u^2;      % mom_flux_paper = 896.89;

u_exp = sqrt(mom_flux./(data.Rho*data.Mw));

m_dot_exp = data.Rho*data.Mw.*u_exp*A_test;

q_dot_exp = m_dot_exp.*( (data.Cp/data.Mw).*T_small' - ( (data.Cp(1,1)/data.Mw)*T_small(1)*adjust_vect ) );

figure()
contourf(data.P*1e-5,data.T,m_dot_exp); 
title('Slab 69x30: mass flow rate momentum flux analogy')
xlabel('Pressure $[bar]$')
ylabel('Temperature $[K]$')
c = colorbar;
c.Label.String = 'm_{N_2} [kg/s]';

figure()
contourf(data.P*1e-5,data.T,1e-3*q_dot_exp/1);     % area di riferimento 1 m^2
title('Slab 69x30: energy flux momentum flux analogy')
xlabel('Pressure $[bar]$')
ylabel('Temperature $[K]$')
c = colorbar;
c.Label.String = 'q_{N_2} [kW/m^2]';

%% 69x69 momentum flux analogy
A_test = 69*69*1e-6;

mom_flux = rho*u^2;      % mom_flux_paper = 896.89;

u_exp = sqrt(mom_flux./(data.Rho*data.Mw));

m_dot_exp = data.Rho*data.Mw.*u_exp*A_test;

q_dot_exp = m_dot_exp.*( (data.Cp/data.Mw).*T_small' - ( (data.Cp(1,1)/data.Mw)*T_small(1)*adjust_vect ) );

figure()
contourf(data.P*1e-5,data.T,m_dot_exp); 
title('Slab 69x69: mass flow rate momentum flux analogy')
xlabel('Pressure $[bar]$')
ylabel('Temperature $[K]$')
c = colorbar;
c.Label.String = 'm_{N_2} [kg/s]';

figure()
contourf(data.P*1e-5,data.T,1e-3*q_dot_exp/1);     % area di riferimento 1 m^2
title('Slab 69x69: energy flux momentum flux analogy')
xlabel('Pressure $[bar]$')
ylabel('Temperature $[K]$')
c = colorbar;
c.Label.String = 'q_{N_2} [kW/m^2]';

%% 30x30 Reynolds analogy
L_test = 30*1e-3;

Re = rho*u*L_test/mu;     %Re_paper = 65797

u_exp = Re*(data.mu./( data.Rho*data.Mw*L_test ));

m_dot_exp = data.Rho*data.Mw.*u_exp*L_test^2;

adjust_vect = ones(length(T_small),length(P_small));

q_dot_exp = m_dot_exp.*( (data.Cp/data.Mw).*T_small' - ( (data.Cp(1,1)/data.Mw)*T_small(1)*adjust_vect ) );

figure()
contourf(data.P*1e-5,data.T,m_dot_exp); 
title('Slab 30x30: mass flow rate Reynolds analogy')
xlabel('Pressure $[bar]$')
ylabel('Temperature $[K]$')
c = colorbar;
c.Label.String = 'm_{N_2} [kg/s]';

figure()
contourf(data.P*1e-5,data.T,1e-3*q_dot_exp/1);     % area di riferimento 1 m^2
title('Slab 30x30: energy flux Reynolds analogy')
xlabel('Pressure $[bar]$')
ylabel('Temperature $[K]$')
c = colorbar;
c.Label.String = 'q_{N_2} [kW/m^2]';

%% 49x30 Reynolds analogy
A_test = 30*49*1e-6;

L_test= 4*30*49*1e-6/ ((49*2+60)*1e-3);

Re = rho*u*L_test/mu;     %Re_paper = 65797

u_exp = Re*(data.mu./( data.Rho*data.Mw*L_test ));

m_dot_exp = data.Rho*data.Mw.*u_exp*A_test;

adjust_vect = ones(length(T_small),length(P_small));

q_dot_exp = m_dot_exp.*( (data.Cp/data.Mw).*T_small' - ( (data.Cp(1,1)/data.Mw)*T_small(1)*adjust_vect ) );

figure()
contourf(data.P*1e-5,data.T,m_dot_exp); 
title('Slab 49x30: mass flow rate Reynolds analogy')
xlabel('Pressure $[bar]$')
ylabel('Temperature $[K]$')
c = colorbar;
c.Label.String = 'm_{N_2} [kg/s]';

figure()
contourf(data.P*1e-5,data.T,1e-3*q_dot_exp/1);     % area di riferimento 1 m^2
title('Slab 49x30: energy flux Reynolds analogy')
xlabel('Pressure $[bar]$')
ylabel('Temperature $[K]$')
c = colorbar;
c.Label.String = 'q_{N_2} [kW/m^2]';

%% 69x30 Reynolds analogy
A_test = 30*69*1e-6;

L_test= 4*30*69*1e-6/ ((69*2+60)*1e-3);

Re = rho*u*L_test/mu;     %Re_paper = 65797

u_exp = Re*(data.mu./( data.Rho*data.Mw*L_test ));

m_dot_exp = data.Rho*data.Mw.*u_exp*A_test;

adjust_vect = ones(length(T_small),length(P_small));

q_dot_exp = m_dot_exp.*( (data.Cp/data.Mw).*T_small' - ( (data.Cp(1,1)/data.Mw)*T_small(1)*adjust_vect ) );

figure()
contourf(data.P*1e-5,data.T,m_dot_exp); 
title('Slab 69x30: mass flow rate Reynolds analogy')
xlabel('Pressure $[bar]$')
ylabel('Temperature $[K]$')
c = colorbar;
c.Label.String = 'm_{N_2} [kg/s]';

figure()
contourf(data.P*1e-5,data.T,1e-3*q_dot_exp/1);     % area di riferimento 1 m^2
title('Slab 69x30: energy flux Reynolds analogy')
xlabel('Pressure $[bar]$')
ylabel('Temperature $[K]$')
c = colorbar;
c.Label.String = 'q_{N_2} [kW/m^2]';
%% 69x69 Reynolds analogy
L_test = 69*1e-3;

Re = rho*u*L_test/mu;     %Re_paper = 65797

u_exp = Re*(data.mu./( data.Rho*data.Mw*L_test ));

m_dot_exp = data.Rho*data.Mw.*u_exp*L_test^2;

adjust_vect = ones(length(T_small),length(P_small));

q_dot_exp = m_dot_exp.*( (data.Cp/data.Mw).*T_small' - ( (data.Cp(1,1)/data.Mw)*T_small(1)*adjust_vect ) );

figure()
contourf(data.P*1e-5,data.T,m_dot_exp); 
title('Slab 69x69: mass flow rate Reynolds analogy')
xlabel('Pressure $[bar]$')
ylabel('Temperature $[K]$')
c = colorbar;
c.Label.String = 'm_{N_2} [kg/s]';

figure()
contourf(data.P*1e-5,data.T,1e-3*q_dot_exp/1);     % area di riferimento 1 m^2
title('Slab 69x69: energy flux Reynolds analogy')
xlabel('Pressure $[bar]$')
ylabel('Temperature $[K]$')
c = colorbar;
c.Label.String = 'q_{N_2} [kW/m^2]';

%% 30x30 mass flux analogy
L_test = 30*1e-3;

G = rho*u;      % G_paper = 40.54;

u_exp = G./(data.Rho*data.Mw);

m_dot_exp = data.Rho*data.Mw.*u_exp*L_test^2;

adjust_vect = ones(length(T_small),length(P_small));

q_dot_exp = m_dot_exp.*( (data.Cp/data.Mw).*T_small' - ( (data.Cp(1,1)/data.Mw)*T_small(1)*adjust_vect ) );

figure()
contourf(data.P*1e-5,data.T,m_dot_exp); 
title('Slab 30x30: mass flow rate mass flux analogy')
xlabel('Pressure $[bar]$')
ylabel('Temperature $[K]$')
c = colorbar;
c.Label.String = 'm_{N_2} [kg/s]';

figure()
contourf(data.P*1e-5,data.T,1e-3*q_dot_exp/1);     % area di riferimento 1 m^2
title('Slab 30x30: energy flux mass flux analogy')
xlabel('Pressure $[bar]$')
ylabel('Temperature $[K]$')
c = colorbar;
c.Label.String = 'q_{N_2} [kW/m^2]';

%% 49x30 mass flux analogy
A_test = 30*49*1e-6;

G = rho*u;      % G_paper = 40.54;

u_exp = G./(data.Rho*data.Mw);

m_dot_exp = data.Rho*data.Mw.*u_exp*A_test;

adjust_vect = ones(length(T_small),length(P_small));

q_dot_exp = m_dot_exp.*( (data.Cp/data.Mw).*T_small' - ( (data.Cp(1,1)/data.Mw)*T_small(1)*adjust_vect ) );

figure()
contourf(data.P*1e-5,data.T,m_dot_exp); 
title('Slab 49x30: mass flow rate mass flux analogy')
xlabel('Pressure $[bar]$')
ylabel('Temperature $[K]$')
c = colorbar;
c.Label.String = 'm_{N_2} [kg/s]';

figure()
contourf(data.P*1e-5,data.T,1e-3*q_dot_exp/1);     % area di riferimento 1 m^2
title('Slab 49x30: energy flux mass flux analogy')
xlabel('Pressure $[bar]$')
ylabel('Temperature $[K]$')
c = colorbar;
c.Label.String = 'q_{N_2} [kW/m^2]';

%% 69x30 mass flux analogy
A_test = 30*69*1e-6;

G = rho*u;      % G_paper = 40.54;

u_exp = G./(data.Rho*data.Mw);

m_dot_exp = data.Rho*data.Mw.*u_exp*A_test;

adjust_vect = ones(length(T_small),length(P_small));

q_dot_exp = m_dot_exp.*( (data.Cp/data.Mw).*T_small' - ( (data.Cp(1,1)/data.Mw)*T_small(1)*adjust_vect ) );

figure()
contourf(data.P*1e-5,data.T,m_dot_exp); 
title('Slab 69x30: mass flow rate mass flux analogy')
xlabel('Pressure $[bar]$')
ylabel('Temperature $[K]$')
c = colorbar;
c.Label.String = 'm_{N_2} [kg/s]';

figure()
contourf(data.P*1e-5,data.T,1e-3*q_dot_exp/1);     % area di riferimento 1 m^2
title('Slab 69x30: energy flux mass flux analogy')
xlabel('Pressure $[bar]$')
ylabel('Temperature $[K]$')
c = colorbar;
c.Label.String = 'q_{N_2} [kW/m^2]';

%% 69x69 mass flux analogy
L_test = 69*1e-3;

G = rho*u;      % G_paper = 40.54;

u_exp = G./(data.Rho*data.Mw);

m_dot_exp = data.Rho*data.Mw.*u_exp*L_test^2;

adjust_vect = ones(length(T_small),length(P_small));

q_dot_exp = m_dot_exp.*( (data.Cp/data.Mw).*T_small' - ( (data.Cp(1,1)/data.Mw)*T_small(1)*adjust_vect ) );

figure()
contourf(data.P*1e-5,data.T,m_dot_exp); 
title('Slab 69x69: mass flow rate mass flux analogy')
xlabel('Pressure $[bar]$')
ylabel('Temperature $[K]$')
c = colorbar;
c.Label.String = 'm_{N_2} [kg/s]';

figure()
contourf(data.P*1e-5,data.T,1e-3*q_dot_exp/1);     % area di riferimento 1 m^2
title('Slab 69x69: energy flux mass flux analogy')
xlabel('Pressure $[bar]$')
ylabel('Temperature $[K]$')
c = colorbar;
c.Label.String = 'q_{N_2} [kW/m^2]';

%% 30x30 energy  flux analogy
L_test = 30*1e-3;

energy_flux = rho*u*c_p*T;      % energy_flux_paper = 284.6*1e6;

u_exp = energy_flux./( ((data.Rho*data.Mw).*(data.Cp/data.Mw)) .* (T_vect'*ones(1,length(P_vect))) );

m_dot_exp = data.Rho*data.Mw.*u_exp*L_test^2;

adjust_vect = ones(length(T_small),length(P_small));

q_dot_exp = m_dot_exp.*( (data.Cp/data.Mw).*T_small' - ( (data.Cp(1,1)/data.Mw)*T_small(1)*adjust_vect ) );

figure()
contourf(data.P*1e-5,data.T,m_dot_exp); 
title('Slab 30x30: mass flow rate energy flux analogy')
xlabel('Pressure $[bar]$')
ylabel('Temperature $[K]$')
c = colorbar;
c.Label.String = 'm_{N_2} [kg/s]';

figure()
contourf(data.P*1e-5,data.T,1e-3*q_dot_exp/1);     % area di riferimento 1 m^2
title('Slab 30x30: energy flux energy flux analogy')
xlabel('Pressure $[bar]$')
ylabel('Temperature $[K]$')
c = colorbar;
c.Label.String = 'q_{N_2} [kW/m^2]';

%% 49x30 energy  flux analogy
A_test = 30*49*1e-6;

energy_flux = rho*u*c_p*T;      % energy_flux_paper = 284.6*1e6;

u_exp = energy_flux./( ((data.Rho*data.Mw).*(data.Cp/data.Mw)) .* (T_vect'*ones(1,length(P_vect))) );

m_dot_exp = data.Rho*data.Mw.*u_exp*A_test;

adjust_vect = ones(length(T_small),length(P_small));

q_dot_exp = m_dot_exp.*( (data.Cp/data.Mw).*T_small' - ( (data.Cp(1,1)/data.Mw)*T_small(1)*adjust_vect ) );

figure()
contourf(data.P*1e-5,data.T,m_dot_exp); 
title('Slab 49x30: mass flow rate energy flux analogy')
xlabel('Pressure $[bar]$')
ylabel('Temperature $[K]$')
c = colorbar;
c.Label.String = 'm_{N_2} [kg/s]';

figure()
contourf(data.P*1e-5,data.T,1e-3*q_dot_exp/1);     % area di riferimento 1 m^2
title('Slab 49x30: energy flux energy flux analogy')
xlabel('Pressure $[bar]$')
ylabel('Temperature $[K]$')
c = colorbar;
c.Label.String = 'q_{N_2} [kW/m^2]';

%% 49x30 energy  flux analogy
A_test = 30*69*1e-6;

energy_flux = rho*u*c_p*T;      % energy_flux_paper = 284.6*1e6;

u_exp = energy_flux./( ((data.Rho*data.Mw).*(data.Cp/data.Mw)) .* (T_vect'*ones(1,length(P_vect))) );

m_dot_exp = data.Rho*data.Mw.*u_exp*A_test;

adjust_vect = ones(length(T_small),length(P_small));

q_dot_exp = m_dot_exp.*( (data.Cp/data.Mw).*T_small' - ( (data.Cp(1,1)/data.Mw)*T_small(1)*adjust_vect ) );

figure()
contourf(data.P*1e-5,data.T,m_dot_exp); 
title('Slab 69x30: mass flow rate energy flux analogy')
xlabel('Pressure $[bar]$')
ylabel('Temperature $[K]$')
c = colorbar;
c.Label.String = 'm_{N_2} [kg/s]';

figure()
contourf(data.P*1e-5,data.T,1e-3*q_dot_exp/1);     % area di riferimento 1 m^2
title('Slab 69x30: energy flux energy flux analogy')
xlabel('Pressure $[bar]$')
ylabel('Temperature $[K]$')
c = colorbar;
c.Label.String = 'q_{N_2} [kW/m^2]';

%% 69x69 energy  flux analogy
L_test = 69*1e-3;

energy_flux = rho*u*c_p*T;      % energy_flux_paper = 284.6*1e6;

u_exp = energy_flux./( ((data.Rho*data.Mw).*(data.Cp/data.Mw)) .* (T_vect'*ones(1,length(P_vect))) );

m_dot_exp = data.Rho*data.Mw.*u_exp*L_test^2;

adjust_vect = ones(length(T_small),length(P_small));

q_dot_exp = m_dot_exp.*( (data.Cp/data.Mw).*T_small' - ( (data.Cp(1,1)/data.Mw)*T_small(1)*adjust_vect ) );

figure()
contourf(data.P*1e-5,data.T,m_dot_exp); 
title('Slab 69x69: mass flow rate energy flux analogy')
xlabel('Pressure $[bar]$')
ylabel('Temperature $[K]$')
c = colorbar;
c.Label.String = 'm_{N_2} [kg/s]';

figure()
contourf(data.P*1e-5,data.T,1e-3*q_dot_exp/1);     % area di riferimento 1 m^2
title('Slab 69x69: energy flux energy flux analogy')
xlabel('Pressure $[bar]$')
ylabel('Temperature $[K]$')
c = colorbar;
c.Label.String = 'q_{N_2} [kW/m^2]';