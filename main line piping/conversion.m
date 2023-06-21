clc
clear all
T1 = 298;                                     
P1 = 30;  
R_N2 = 8314/28; 
rho=P1/R_N2 *T1;
m_dot_kg_s=60*10^-3;
m_dot_m3_h= m_dot_kg_s*3600/rho;