%1 l=1 dm^3
%conversione 60 g/s a l/min
% densita azoto a 650 K e 200 bar è 90 kg/m^3
% densita azoto in condizioni ambiente è 1.25 kg/m^3
mdotN21=30*78*1e-3/(45*1e-3)  % conversione da g/s a l/ min
mdotN22=3600*(78*1e-3/45)*(100/1.01325)*(273.15/650)  % conversione da g/s a nm^3/h

mdotN23=mdotN21*(100/1.01325)*(273.15/650) % conversione da g/s a Nl/ min