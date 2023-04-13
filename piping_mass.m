rho_steel=8027; %[kg/m^3]
Ltot=2.2;       %[m]


%% Case with 10 mm tube
d_tube1=10*1e-3;
A_cross1=(pi/4)*d_tube1^2;
V_tube1=A_cross1*Ltot;
m_tube1=V_tube1*rho_steel;    %[kg]

%% Case with 12 mm tube
d_tube2=12*1e-3;
A_cross2=(pi/4)*d_tube2^2;
V_tube2=A_cross2*Ltot;
m_tube2=V_tube2*rho_steel;     %[kg]

%% Case with 1/2 inch tube
d_tube3=12.7*1e-3;
A_cross3=(pi/4)*d_tube3^2;
V_tube3=A_cross3*Ltot;
m_tube3=V_tube3*rho_steel;      %[kg]

%% Case with 3/4 inch tube
d_tube4=19.05*1e-3;
A_cross4=(pi/4)*d_tube4^2;
V_tube4=A_cross4*Ltot;
m_tube4=V_tube4*rho_steel;      %[kg]