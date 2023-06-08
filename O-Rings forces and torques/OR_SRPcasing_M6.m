clear
clc

D_bolt = 6*1e-3;          
K = 0.2;                         % K is a factor that accounts for the material of the bolt, lubricants...
T_clamping = 4.55;               % Grade 8.8: il numero dipende dai trattamenti termici effettuati
F_bolt = T_clamping/(D_bolt*K);

n_bolt = 4;
F_tot_bolt = F_bolt*n_bolt

d_circ = 44.12;                                % O-Ring diameter [mm]
l_circ = pi*d_circ;
F_max = ((35/2.205)*9.81)/25.4;                % Force per unit length [N/mm]
F_OR = F_max*l_circ;

P_int = 12*1e5;
A_int = 0.25*pi*(d_circ*1e-3)^2;
F_press = P_int*A_int;

F_tot_needed = F_press + F_OR