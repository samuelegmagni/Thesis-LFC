T = [295:1:600];
P = [10:1:12];

T_beforeinj=295.32;
T_fin=600;
P_beforeinj=12;
data = nistdata('N2',T,P); 
cp_N2 = data.Cp/data.Mw;
cp_N2_in = cp_N2(find(T==round(T_beforeinj)),find(abs(P - round(P_beforeinj,1)) < 0.001)); 
cp_N2_fin = cp_N2(find(T==round(T_fin)),find(abs(P - round(P_beforeinj,1)) < 0.001)); 