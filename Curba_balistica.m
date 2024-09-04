clc
clear all
close all
format long

Gv = 1;
Gc = 1;
Gt = 1;
Gd = 1;
if Gd==1, Gv=0; Gc=0; Gt=0; end;

g = 9.80665; % acceleratia gravitationala
ro = 7850; % densitatea otelului
r = 0.13; % raza proiectilului
m = 4*pi*ro*r^3 / 3; % masa proiectilului
G = m * g; % greutatea proiectilului

v0 = 1100; % viteza initiala
alpha0 = 43; % unghiul de lansare
eta = 1.81 * 1e-5; % coeficientul de vascozitate
b1 = 6*pi*eta*r; % coeficientul b1
c = 0.469; % coeficientul de forma
ro0 = 1.22; % densitatea aerului
b2 = 2*pi*c*r^2*ro0; % coeficientul b2

t0 = 0; tf = 2*v0/(g*sind(alpha0)); 
N = 1501; dt = tf/(N - 1);
t = linspace(t0, tf, N);

vx = zeros(1,N); vy = vx;
x = zeros(1,N); y = x;
vx(1) = v0 * cosd(alpha0);
vy(1) = v0 * sind(alpha0);
x(1) = 0; y(1) = 0;
for i = 1:N-1
    aux = 1 - dt*(b1 + b2*sqrt(vx(i)^2 + vy(i)^2))/m;
    vx(i + 1) = vx(i)*aux;
    vy(i + 1) = vy(i)*aux - g*dt;
    x(i + 1) = x(i) + vx(i)*dt;
    y(i + 1) = y(i) + vy(i)*dt;
    if y(i) < 0 break; end;
end

% componentele vitezei
if Gv==1
    figure(1);
    plot(t,vx,'-r',t,vy,'-b');
    xlabel('t(s)'); ylabel('v(m/s)'); grid on;
    title('Componentele vitezei in functie de timp');
    legend('vx','vy','location','northwest');
end

% legile de miscare
if Gc==1
    figure(2);
    plot(t, x/1e3, '-r', t, y/1e3, '-b');
    xlabel('t(s)'); ylabel('coord(km)'); grid on;
    title('Componentele vitezei in functie de timp');
    legend('x','y','location','northwest');
end

% traiectoria
if Gt==1
    figure(3)
    plot(x/1e3, y/1e3, '-k', 'LineWidth',2);
    xlabel('x(km)'); ylabel('y(km)'); grid on;
    title('Curba balistica');
    axis equal; axis tight;
end

tf = t(i); % timpul de zbor
b = x(i); % bataia
h = max(y); % inaltimea maxima
tu = t(y==h); % timpul de urcare
tc = tf - tu; % timpul de coborare
Q = 1/2*m*(v0^2 - vx(i)^2 - vy(i)^2); % caldura produsa prin frecare
afis = ['Timpul de zbor: ', num2str(tf), 's']; disp(afis);
afis = ['Bataia proiectilului: ', num2str(b/1e3), 'km']; disp(afis);
afis = ['Altitudinea maxima: ', num2str(h/1e3), 'km']; disp(afis);
afis = ['Timpul de urcare: ', num2str(tu), 's']; disp(afis);
afis = ['Timpul de coborare: ', num2str(tc), 's']; disp(afis);
afis = ['Caldura produsa: ', num2str(Q/1e6), 'MJ']; disp(afis);

%simularea dinamica
if Gd==1
    figure(4)
    set(4, 'Position', [50 50 850 600]);
    tic; simt = 0; % porneste cronometrul si retine timpul initial
    while simt < tf
        plot(x/1e3, y/1e3, '-b'); hold on;
        xlabel('x(km)'); ylabel('y(km)'); grid;
        title('Simularea miscarii');
        axis equal; axis tight;
        index = abs(t-simt)==min(abs(t-simt)); % cauta cel mai apropiat t din discretizare
        plot(x(index)/1e3, y(index)/1e3, '.b', 'MarkerSize', 10); hold off;
        text(b/2/1e3, h/3/1e3, ['vx=', num2str(round(vx(index))), 'm/s']);
        text((b/2-b/5)/1e3, h/3/1e3, ['vt=', num2str(round(t(index))), 's']);
        text((b/2+b/5)/1e3, h/3/1e3, ['vy=', num2str(round(vy(index))), 'm/s']);
        pause(1e-3);
        simt = toc;
    end
end