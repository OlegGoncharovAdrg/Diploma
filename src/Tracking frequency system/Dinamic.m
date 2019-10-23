clear; close all; clc
 
T = 3*3.15e-3; % такт поступления отчетов корреляционных сумм
Tmax = 200;   % время моделирования
qcno_dB = 35;  % отношение сигнал/шум
q_cno = 10^(qcno_dB/10);

t = T:T:Tmax;  % массив временных отчетов
N = length(t);
 
F = [1 T;      % коэффициент шумов наблюдений
 0 1];
 
G = [0 0; % коэффициент формирующего шума
 0 T]; 
hold on
 for cnt = 1:2
   if cnt==1 % вывод флуктационной ошибки на график
       Frequency_traking_system_8_P_fluctuation; 
   else
       if cnt == 2 % вывод динамической ошибки на график
       sigma_a = 1.2; % СКЗ ускорения
       alpha = 0.1; % Коэффициент размерности
       Sfi = 2*(33*sigma_a)^2*alpha; % Спектральная плотность формирующего шума
       D_teta_fi = Sfi/T*1; % Дисперсия формирующего шума - динамическая ошибка
       Dteta_omega = 16*qcno_dB^3*T^3*(1+1/(2*qcno_dB*T))*1; % дисперсия шумов наблюдений - флуктационная 
       Sd = 4*q_cno^2*T^3; % крутизна дискриминационной характеристики
       ksi = sqrt(D_teta_fi) * randn(1, N); % Реализация формирующего шума
       eta = sqrt(Dteta_omega) * randn(1, N); % Реализация шумов наблюдений
 
     Band = 0.1:0.1:10;  % Полоса СС
     Band_for_plot = 10; % Полоса, при которой вывести графики
     RMS_Omega = nan(1,length(Band));
 
     for i = 1:length(Band) 
     K = nan(2, 1);
     K(1) = 8/3 * Band(i) * T; % Коэффициенты СС
     K(2) = 32/9 * Band(i)^2 * T; 
     Phi_0 = 5/180*pi;
 
     Xist = [0; 0]; % Начальные условия
     Xest = [10; 0];
 
     I = nan(1, N); Q = nan(1, N);
     Ud = nan(1, N-1); ErrOmega = nan(1, N);
     Phi = nan(1, N); Phi_op = nan(1, N);  ErrPhi = nan(1, N); 
     Phi(1)=Phi_0;
     Phi_op(1)=Phi(1);
     ErrOmega(1) = Xist(1) - Xest(1);
     ErrPhi(1)=Phi(1)-Phi_op(1);
     I(1) = 2 * q_cno * T * sinc(ErrOmega(1) *T/2)* cos(ErrPhi(1)+ErrOmega(1) *T/2); 
     Q(1) = -2 * q_cno * T * sinc(ErrOmega(1)*T/2)* sin(ErrPhi(1)+ErrOmega(1)*T/2);
 
     for k = 2:N
         Xist = F*Xist + G*[0;ksi(k)]; % Развитие оцениваемого процесса
         Phi(k)=Phi(k-1)+T*Xist(1); % формирование оценки и экстраполированной оценки фазы
         Phi_op(k)=Phi_op(k-1)+T*Xest(1);
         ErrPhi(k)=Phi(k)-Phi_op(k); % ошибка рассогласования по фазе
         Xextr = F*Xest; % Экстраполяция оцениваемого процесса    
         ErrOmega(k) = Xist(1) - Xextr(1); % Ошибка оценивания
         I(k) = 2 * q_cno * T* cos(ErrPhi(k)+ErrOmega(k)*T/2)* sinc(ErrOmega(k)*T/2); 
         Q(k) = -2 * q_cno * T* sin(ErrPhi(k)+ErrOmega(k)*T/2)* sinc(ErrOmega(k)*T/2);
         Ud(k) = (I(k)*Q(k-1) - I(k-1)*Q(k))/Sd + eta(k); % сигнал с выхода дискриминатора с шумом
         Xest = Xextr + K*Ud(k);
     end
     RMS_Omega(i) = sqrt(mean(ErrOmega.^2));
     end
     figure(2)
     plot(Band, RMS_Omega);
      xlabel('\Deltaf, Гц'); ylabel('СКО_{f}, Гц');
      grid on;
     end
   end
 end