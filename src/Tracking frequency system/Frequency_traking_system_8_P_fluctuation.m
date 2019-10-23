clear; close all; clc

T=3*3.15e-3;
Tmax = 10;   % ����� �������������
qcno_dB = 35;  % ��������� ������/���
q_cno = 10^(qcno_dB/10);

t = T:T:Tmax;  % ������ ��������� �������
N = length(t);

F = [1 T;      % ����������� ����� ����������
 0 1];

G = [0 0; % ����������� ������������ ����
 0 T]; 

Sd = 4*q_cno^2*T^3; % �������� ����������������� ��������������

Band = 0.1:0.1:10;  % ������ ��
Band_for_plot = Band(length(Band)); % ������, ��� ������� ������� �������
hold on 
 for i = 1:length(Band) 
 K = nan(2, 1);
 K(2) = (2*Band(i))^2;
 K(1) = (2*K(2))^(1/2);
 K = K*T;
 CKO_I_n=sqrt(2*q_cno*T);
 Phi_0 = 5/180*pi;

 Xist = [0; 0]; % ��������� �������
 Xest = [0; 0];

 I = nan(1, N); Q = nan(1, N);
 Ud = nan(1, N-1); ErrOmega = nan(1, N);
 Phi = nan(1, N); Phi_op = nan(1, N);  ErrPhi = nan(1, N); 
 Phi(1)=Phi_0;
 Phi_op(1)=Phi(1);
 ErrOmega(1) = Xist(1) - Xest(1);
 ErrPhi(1)=Phi(1)-Phi_op(1);

 I(1) = 2 * q_cno * T * sinc(ErrOmega(1) *T/2)* cos(ErrPhi(1)+ErrOmega(1) *T/2)+CKO_I_n*randn; 
 Q(1) = -2 * q_cno * T * sinc(ErrOmega(1)*T/2)* sin(ErrPhi(1)+ErrOmega(1)*T/2)+CKO_I_n*randn;
 SKO_f(i) = (ErrOmega(1)/(2*pi))^2/N;
     for k = 2:N
         Phi(k)=Phi(k-1)+T*Xist(1); % ������������ ������ � ������������������ ������ ����
         Phi_op(k)=Phi_op(k-1)+T*Xest(1);
         ErrPhi(k)=Phi(k)-Phi_op(k); % ������ ��������������� �� ����
         Xextr = F*Xest; % ������������� ������������ ��������    
         ErrOmega(k) = Xist(1) - Xextr(1); % ������ ����������
         SKO_f(i) = SKO_f(i)+(ErrOmega(k)/(2*pi))^2/N;
         I(k) = 2 * q_cno * T* cos(ErrPhi(k)+ErrOmega(k)*T/2)* sinc(ErrOmega(k)*T/2)+CKO_I_n*randn; 
         Q(k) = -2 * q_cno * T* sin(ErrPhi(k)+ErrOmega(k)*T/2)* sinc(ErrOmega(k)*T/2)+CKO_I_n*randn;
         Ud(k) = (I(k)*Q(k-1) - I(k-1)*Q(k))/Sd ; % ������ � ������ �������������� � �����
         Xest = Xextr + K*Ud(k);
     end
 end
 figure(1);
 plot(Band, SKO_f);
  xlabel('\Deltaf, Hz'); ylabel('���_{f}, Hz');
 grid on;
 
 
 