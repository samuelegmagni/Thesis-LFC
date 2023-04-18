%% Mixing chamber sizing
clear
clc


Pgi=;
Pgf=12*1e5;   % pressure before injector 
Pburst=4*Pgi;
rho_steel=8000;
sigma_ult=505*1e6;
sigma_snerv=215*1e6;
htank=30*1e-2;
Vtank=;  %consider volume of N2 and SRP
rtank=sqrt(Vtank/(pi*htank));
t_des= Pburst*r/sigma_snerv;
t-ASME=Pburst*r/(sigma_ult-0.6*Pburst);