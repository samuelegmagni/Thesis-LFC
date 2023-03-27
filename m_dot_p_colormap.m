clear
clc
set(0,'DefaultTextFontSize',12);              % Settings for the plot
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultLegendInterpreter','Latex');
set(groot,'DefaultAxesTickLabelInterpreter','Latex');
set(0,'DefaultTextInterpreter','Latex');
set(0,'DefaultLegendFontSize',12);

T = [300:5:700];
P = [1:2:30];
T_amb = 298.15;
cp_g = 2363;                 % J/kgK
T_fl = 0.8*2305; 

data = nistdata('N2',T,P); 

L_test = 30*1e-3;
rho = 1.832;
u = 22.12;
mom_flux = rho*u^2;      % mom_flux_paper = 896.89;

u_exp = sqrt(mom_flux./(data.Rho*data.Mw));

m_dot_N2 = data.Rho*data.Mw.*u_exp*L_test^2;

cp_N2 = data.Cp/data.Mw;

m_dot_p = zeros(length(T),length(P));

for i = 1 : length(T)

    for m = 1 : length(P)

        f = @(x) m_dot_N2(i,m)*cp_N2(i,m)*(T(i) - T_amb) + x*cp_g*(T_fl) - (m_dot_N2(i,m) + x) * ( (cp_N2(i,m)*m_dot_N2(i,m))/(m_dot_N2(i,m) + x) + (cp_g*x)/(m_dot_N2(i,m) + x) )*(T(i));
        z = 1.2*fzero(f,0.5) ;
        m_dot_p(i,m) = z; 

    end

end

figure()
contourf(data.P*1e-5,data.T,m_dot_N2*1e3); 
title('Slab 30x30: mass flow rate momentum flux analogy')
xlabel('Pressure $[bar]$')
ylabel('Temperature $[K]$')
c = colorbar;
c.Label.String = 'm_{dot,N_2} [g/s]';

figure()
contourf(data.P*1e-5,data.T,m_dot_p*1e3); 
title('Slab 30x30: mass flow rate momentum flux analogy')
xlabel('Pressure $[bar]$')
ylabel('Temperature $[K]$')
c = colorbar;
c.Label.String = 'm_{dot,p} [g/s]';