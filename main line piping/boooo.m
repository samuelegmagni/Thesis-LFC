clc
clear all
m_dot_N2= 78*10e-3;
R_N2 = 8314/28;    
rhoN2= 115;
G_g = rhoN2/1000;                       % Nitrogen specific gravity [-]
q_N2 = (m_dot_N2/rhoN2)*1000;           % Nitrogen volumetric flow rate [L/s]
deltaP= 100-95.43;
CV= sqrt(G_g*(q_N2*60)^2/(14.42^2*deltaP));



% m_dot_N2= 78*10e-3;
% R_N2 = 8314/28;
% rhoN2= 113;
% G_g = rhoN2/1000;                       % Nitrogen specific gravity [-]
% q_N2 = (m_dot_N2/rhoN2)*1000;   
% Cv=1.41;
% deltaP = (G_g*(q_N2*60)^2)/(14.42*Cv)^2;  