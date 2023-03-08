clear
clc

format long

set(0,'DefaultTextFontSize',12);              % Settings for the plot
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultLegendInterpreter','Latex');
set(0,'DefaultAxesTickLabelInterpreter','Latex');
set(0,'DefaultTextInterpreter','Latex');
set(0,'DefaultLegendFontSize',12);

%% GE-C parameters

R_CC = 79.2 * 1e-3;
R_t = 14.42 * 1e-3;
eps_c = 30.16;
eps_e = 125;
L_CC = 22*1e-3;
m_dot = 0.798;
rho = 1.832;
T = 2787;
u = 22.12;
M_m = 21.23; %g/mol
c_p = 2519;  %J/kgK
mu = 9.76*1e-5;
lambda = 0.3372;

%% Test section

L = 

Re = rho*u*2*R_CC/mu;     %Re_paper = 65797

u_exp_vect = zeros(1)
% ciaos
u_exp = sqrt(Re*)