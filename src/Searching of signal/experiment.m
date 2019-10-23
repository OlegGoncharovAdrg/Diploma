if strcmp(signal, 'off')
I = std_IQ*randn(M, N) ;
Q = std_IQ*randn(M, N) ;
X2 = I.^2 + Q.^2;
X2max(j) = max(max(X2));
elseif strcmp(signal, 'on')
dtau = tau(j) - tau_tilda_m;
domega = omega(j) - omega_tilda_m;
I = A*L/2 * ro(dtau) .* sinc(domega*T/2 /pi) .* cos(domega*T/2 ...
+ dphi(j)) + std_IQ*randn(M, N) ;
Q =- A*L/2 * ro(dtau) .* sinc(domega*T/2 /pi) .* sin(domega*T/2 ...
+ dphi(j)) + std_IQ*randn(M, N) ;
X2 = I.^2 + Q.^2;
X2max(j) = max(max(X2));
end



  
