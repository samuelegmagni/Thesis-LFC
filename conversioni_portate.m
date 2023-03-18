%1 l=1 dm^3
%conversione 60 g/s a l/min
% densita azoto a 650 K e 200 bar è 90 kg/m^3
% densita azoto in condizioni ambiente è 1.25 kg/m^3
mdotN21=60*60*1e-3/(90*1e-3)  % conversione da g/s a l/ min
mdotN22=3600*(60*1e-3/90)*(200/1)*(273/650)  % conversione da g/s a nm^3/h

mdotN23=mdotN21*(200/1)*(273/650) % conversione da g/s a Nl/ min