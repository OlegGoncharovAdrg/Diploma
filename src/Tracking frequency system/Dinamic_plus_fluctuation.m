clear; close all; clc;
Frequency_traking_system_8_P_fluctuation;
Dinamic;

SKO = sqrt(RMS_Omega.^2+SKO_f.^2);

figure(3)
plot(Band, SKO);
xlabel('\Deltaf, Гц'); ylabel('СКО_{f}, Гц');
grid on;
[miny, minx] = min(SKO);
x = [0.5 minx*0.01];
y = [0.5 miny*0.14];
annotation('textarrow',x,y,'String',['Минимум СКО = ' num2str(miny) ' Гц  при полосе CC = ' num2str(0.1*minx) ' Гц']);