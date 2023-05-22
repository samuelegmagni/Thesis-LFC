% Calculation of pressure losses through inlet and outlet tubes

%% Case 1 with data from flow simulation
clc
clear all
P1=12*1e5;
rho1=13.58;
v1=68;
P2=P1-rho1*v1^2;
rho2=5.5;
v2=269.98;
P3=P2-0.5*rho2*v2^2;


%% Case 2 with formula for incompressible flows
clc
clear all
mdot=60*1e-3;
Atube=pi*(4.5*1e-3)^2;
T = [290:5:600];
P = [9:0.1:12];
data = nistdata('N2',T,P);
Tmix = 295;                             
Ptank = 12;                        
rho_N2 = data.Rho*data.Mw; 
rho1= rho_N2(find(T==round(Tmix)),find(abs(P - round(Ptank,1)) < 0.001)); 
v1=mdot/(Atube*rho1);
P1=12*1e5;
P2=P1-rho1*v1^2;
Texit=550;
rho2= rho_N2(find(T==round(Texit)),find(abs(P - round(P2*1e-5,1)) < 0.001)); 
v2=mdot/(Atube*rho2);
P3=P2-0.5*rho2*v2^2;



%% Case 3 with formula for compressible flows (P3 not found)
clc
clear all
mdot=60*1e-3;
T = [290:5:600];
P = [9:0.1:12];
data = nistdata('N2',T,P);
Tmix = 295;                             
Ptank = 12;  
K1=1;
cp_N2 = data.Cp/data.Mw;            
cv_N2 = data.Cv/data.Mw; 
gamma_N2 = cp_N2./cv_N2;  
gamma1= gamma_N2(find(T==round(Tmix)),find(abs(P - round(Ptank,1)) < 0.001)); 
MW=28;
R = 8314/28;  
P1=12*1e5;
z1= @(x) mdot-K1*Atube*P1*sqrt((2*MW)/(R*Tmix)*(gamma1/(gamma1-1))*((x/P1)^(2/gamma1)-(x/P1)^((gamma1+1)/gamma1)));
P2=fsolve(z1,11*1e5);

K2=0.5;
Texit=550;
rho2= rho_N2(find(T==round(Texit)),find(abs(P - round(P2*1e-5,1)) < 0.001)); 
gamma2= gamma_N2(find(T==round(Texit)),find(abs(P - round(P2*1e-5,1)) < 0.001)); 
z2= @(x) mdot-K2*Atube*P2*sqrt((2*MW)/(R*Texit)*(gamma2/(gamma2-1))*((x/P2)^(2/gamma2)-(x/P2)^((gamma2+1)/gamma2)));
P3=fsolve(z2,9*1e5);