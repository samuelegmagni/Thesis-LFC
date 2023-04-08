%% Injector sizing
m_dot=60*1e-3                        % Mass flow rate [kg/s]
%T = 260:0.5:293;
%P = 7:0.05:12;
%data = nistdata('N2',T,P);
%rho_N2 = data.Rho*data.Mw;           % Density of Nitrogen [kg/m^3] 
%rho_inj = rho_N2(find(T==round(T12,1)),find(abs(P12 - round(P1,1)) < 0.001));
deltaP_inj=1*10^(5);
rho_inj= 30;
Cd=0.8;
N_inj=20;
D_inj_hole= sqrt((2*sqrt(2)*m_dot)/(pi*N_inj*Cd*sqrt(rho_inj *deltaP_inj)));
v_inj= 4*m_dot/(pi*rho_inj*N_inj*D_inj_hole^2);
A_slab=30*30*10e-6;                  % Area of slab test facility [m^2]
A_inj_hole= pi* (D_inj_hole/2)^2;
A_inj_tot= N_inj*A_inj_hole;