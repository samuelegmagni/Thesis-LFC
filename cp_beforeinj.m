T = [295:1:600];
P = [10:1:12];

T1=295.32;
T2=600;
P1=12;
P2=12;
data = nistdata('N2',T,P); 
cp_N2 = data.Cp/data.Mw;
cp_N2in = cp_N2(find(T==round(T1)),find(abs(P - round(P1,1)) < 0.001)); 
cp_N2fin = cp_N2(find(T==round(T2)),find(abs(P - round(P2,1)) < 0.001)); 