clear all;
close all;

N  = 100;
dx = 1.0/N; 

alpha = 1;
beta  = [1 10 100 1000 10000];

eta = [0:dx:1]';
psi = zeros(N+1,1);

for i=1:5
  psi(:,i) = 4/alpha*atanh(tanh(alpha/4)*exp(-sqrt(alpha*beta(i))*eta));
end

figure(1); 
plot(eta, psi(:,1), 'b-.');

hold;
plot(eta, psi(:,2), 'm--');
plot(eta, psi(:,3), 'r:');
plot(eta, psi(:,4), 'c-.');
plot(eta, psi(:,5), 'k-');
hold off;

figure(2);
plot(eta(2:end-1), (psi(3:end,1)-psi(1:end-2,1))./(eta(3:end)-eta(1:end-2)))

