clear; close all; clc

N0 = -206;
Fd = 1.3e6; % ������� �������������, ���
Td = 1/Fd;
std_y = 10^(N0/10)/(2*Td); % ��� ���� �������
NN = 1;
T = 9.45e-3; % �������� ����������, ��
Pf = 1e-3; % ����������� ��
M = 1; % ����� �� �������
N = 52; % ����� �� ��������
L = T / Td; % ����� ������������
std_IQ = std_y * sqrt(L/2); % ��� ���� �������������� ����
delta_tau = 1/2; % ��� ����� �� ��������, ��������
tau_tilda = (0:N-1)*delta_tau + delta_tau/2; % ������� ��������
tau_max = N*delta_tau;
delta_f = 2/3 / T; % ��� ����� �� �������, ��
omega_tilda = 2*pi*((0:M-1)*delta_f + delta_f/2); % ������� �������
omega_max = 2*pi*M*delta_f;
% �����, �������� ������ �����
[tau_tilda_m, omega_tilda_m] = meshgrid(tau_tilda, omega_tilda);
% ����� ������������� ��� ������ ������ 
J1 = 100000;
% ����� ������������� ��� ������� ������������� �����������
J2 = 10000;
X2max = nan(1, J1); % ������������� ������

signal = 'off';
for j = 1:J1
    experiment;
    if ~mod(100*j/J1, 10)
     fprintf('Task 1: Progress %.0f%%\n', 100*j/J1);
    end
end

X2max_incoher = nan(1, floor(J2/NN)); % ������������� ������
count_k = 1;
for j = 1:floor(J2/NN)
    X2max_incoher(j) = 0;
    for kk = count_k:count_k + NN - 1
    X2max_incoher(j) = X2max_incoher(j) + X2max(kk);
    end
    count_k = count_k + NN;
end

R = std_IQ^2; % ����� ������ �����
while sum(X2max_incoher > R) / floor(J1/NN) > Pf
R = R * 1.0005; % ����������� �� 0.002 ��
end


% ����������� inline-������� ro(dtau)
ro = inline('(1 - abs(dtau)) .* (abs(dtau)<1)', 'dtau');
qcno_dB = 25:0.5:45;
X2max = nan(1, J2); % �������� ������� �����������
X2max_incoher = nan(1, floor(J2/NN)); % ������������� ������
Pd = nan(1, length(qcno_dB));
signal = 'on';

for q = 1:length(qcno_dB)
qcno = 10^(qcno_dB(q)/10); % ������� �� �� � ����
A = 2*std_y * sqrt(qcno*Td); % ������ ��������� ��� ������� �/�
% �������� �������� ��� ������� ������������
tau = tau_max * rand(1, J2);
% �������� ������� ��� ������� ������������
omega = omega_max * rand(1, J2);
% ��������� ���� � ������ ������������ ��������
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
Pd(q) = sum(X2max_incoher > R) / floor(J2/NN); % ������� ���������� ������ 
end

figure(1); plot(qcno_dB, Pd);
xlabel('q_{c/n0}, dBHz'); ylabel('P_d'); grid on;
title('�������������� ����������� �������������� ����������')