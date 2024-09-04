clc; clear; close all;

m = 0.00148; % 14.8dg = 0.00148kg


%F(t) = 12 * (sin(sqrt(t-4))*sin(sqrt(t-4))) * exp(-0.01 * (t-4)) * 0.001; %0.001N = 1 mN

t0 = 4; tf = 34;
N = 1501;
t = linspace(t0, tf, N);
dt = t(2) - t(1);

v = zeros(1,N);
x = zeros(1,N);
a=zeros(1,N);
F = zeros(1,N);

v(1) = 0;
x(1) = 0;
F(1) = 0;

for i = 1:N-1
    v(i + 1) = v(i) + (F(i)/m)*dt;
    x(i+1) = x(i) + v(i) * dt;
    F(i + 1) = 12 * (sin(sqrt(t(i+1)-4))*sin(sqrt(t(i+1)-4))) * exp(-0.01 * (t(i+1)-4)) * 0.001;
end
t = t(1:i); v=v(1:i); x=x(1:i);


figure(1)
plot(t,v,'-r');
xlabel('t(s)'); ylabel('v(m/s)'); grid;
title('Viteza in functie de timp');
axis([t0 tf min(v) max(v)]);
legend('v(m/s)');
vmax = max(v);
text(3/5*tf,3/5*vmax,['vf = ',num2str(round(v(i))),' m/s']);









