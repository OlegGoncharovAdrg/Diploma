clear; close all; clc;
Frequency_traking_system_8_P_fluctuation;
Dinamic;

SKO = sqrt(RMS_Omega.^2+SKO_f.^2);

figure(3)
plot(Band, SKO);
xlabel('\Deltaf, ��'); ylabel('���_{f}, ��');
grid on;
[miny, minx] = min(SKO);
x = [0.5 minx*0.01];
y = [0.5 miny*0.14];
annotation('textarrow',x,y,'String',['������� ��� = ' num2str(miny) ' ��  ��� ������ CC = ' num2str(0.1*minx) ' ��']);