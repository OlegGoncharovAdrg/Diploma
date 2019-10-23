clear; close all; clc

N0 = -206;
Fd = 1.3e6; % Частота дискретизации, МГц
Td = 1/Fd;
std_y = 10^(N0/10)/(2*Td); % СКО шума выборки
NN = 1;
T = 9.45e-3; % Интервал накопления, мс
Pf = 1e-3; % Вероятность ЛТ
M = 1; % Ячеек по частоте
N = 52; % Ячеек по задержке
L = T / Td; % Число суммирований
std_IQ = std_y * sqrt(L/2); % СКО шума корреляционных сумм
delta_tau = 1/2; % Шаг ячеек по задержке, символов
tau_tilda = (0:N-1)*delta_tau + delta_tau/2; % Опорные задержки
tau_max = N*delta_tau;
delta_f = 2/3 / T; % Шаг ячеек по частоте, Гц
omega_tilda = 2*pi*((0:M-1)*delta_f + delta_f/2); % Опорные частоты
omega_max = 2*pi*M*delta_f;
% Сетка, задающая центры ячеек
[tau_tilda_m, omega_tilda_m] = meshgrid(tau_tilda, omega_tilda);
% Число экспериментов для поиска порога 
J1 = 100000;
% Число экспериментов для расчета характеристик обнаружения
J2 = 10000;
X2max = nan(1, J1); % Инициализация памяти

signal = 'off';
for j = 1:J1
    experiment;
    if ~mod(100*j/J1, 10)
     fprintf('Task 1: Progress %.0f%%\n', 100*j/J1);
    end
end

X2max_incoher = nan(1, floor(J2/NN)); % Инициализация памяти
count_k = 1;
for j = 1:floor(J2/NN)
    X2max_incoher(j) = 0;
    for kk = count_k:count_k + NN - 1
    X2max_incoher(j) = X2max_incoher(j) + X2max(kk);
    end
    count_k = count_k + NN;
end

R = std_IQ^2; % Очень низкий порог
while sum(X2max_incoher > R) / floor(J1/NN) > Pf
R = R * 1.0005; % Увеличиваем на 0.002 дБ
end


% Определение inline-функции ro(dtau)
ro = inline('(1 - abs(dtau)) .* (abs(dtau)<1)', 'dtau');
qcno_dB = 25:0.5:45;
X2max = nan(1, J2); % Стирание прошлых результатов
X2max_incoher = nan(1, floor(J2/NN)); % Инициализация памяти
Pd = nan(1, length(qcno_dB));
signal = 'on';

for q = 1:length(qcno_dB)
qcno = 10^(qcno_dB(q)/10); % Перевод из дБ в разы
A = 2*std_y * sqrt(qcno*Td); % Расчет амплитуды для данного с/ш
% Истинная задержка для каждого эксперимента
tau = tau_max * rand(1, J2);
% Истинная частота для каждого эксперимента
omega = omega_max * rand(1, J2);
% Начальная фаза в каждом эксперименте случайна
dphi = 2*pi*rand(1, J2);
for j = 1:J2
	experiment;
    if ~mod(100*j/J2, 10)
      fprintf('Task 2: SNR=%.0f dBHz Progress %.0f%%\n', qcno_dB(q), 100*j/J2);
    end
end
count_k = 1;
for j = 1:floor(J2/NN)
    X2max_incoher(j) = 0;
    for kk = count_k:count_k + NN - 1
    X2max_incoher(j) = X2max_incoher(j) + X2max(kk);
    end
    count_k = count_k + NN;
end
Pd(q) = sum(X2max_incoher > R) / floor(J2/NN); % Среднее превышение порога 
end

figure(1); plot(qcno_dB, Pd);
xlabel('q_{c/n0}, dBHz'); ylabel('P_d'); grid on;
title('Характеристики обнаружения некогерентного накопления')