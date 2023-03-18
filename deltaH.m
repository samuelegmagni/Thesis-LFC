deltaHreag=-2129.73; %KW/kg dato dal CEA
% per i prodotti tengo conto di CO,CO2,HCl,H20,N2 e H2 hanno 0
deltaHprod= 0.32848*(-110.53/(28.01*1e-3))+0.10562*(-393.52/(44.01*1e-3))+0.24448*(-92.31/(36.458*1e-3))+0.20013*(-241.83/(18.02*1e-3)); %KW/kg
deltaHreaz= deltaHprod-deltaHreag; %KW/kg
Qcomb=-deltaHreaz; %KW/kg