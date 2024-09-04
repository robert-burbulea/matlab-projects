% Simularea oscilatiilor de amplitudine oarecare ale unui

% pendul gravitational

clear; close all; clc;

g=9.80665; % m/s^2;

lungime=1.5; % m; lungimea firului

omega0=sqrt(g/lungime); % pulsatia oscilatiilor armonice (amplitudini mici)

T0=2*pi/omega0; % perioada proprie oscilatiilor armonice

ti=0; tf=10*T0; N=500; % acasa schimbati 500 cu 1000-2000

t=linspace(ti,tf,N); % sirul momentelor de timp

dt=t(2)-t(1); % pasul de timp

theta_i=pi/3; % unghiul de inceput facut de fir cu verticala y

theta_an=theta_i*sin(omega0*t+pi/2); % solutia armonica

figure(1); % legea unghiulara de miscare

set(1,'Position',[60 50 900 600]);

plot(t,theta_an*180/pi,'-r'); hold on;

xlabel('t/s'); ylabel('theta/grade'); grid;

title('Legea unghiulara de miscare');

theta_num=zeros(1,N); % prealocare solutie numerica

theta_num(1)=theta_i; theta_num(2)=theta_i; % conditii initiale

auxiliar=g/lungime*dt^2; % pentru optimizare ciclu

for i=2:1:N-1

    % relatia de recurenta de ordinul 2:

    theta_num(i+1)=2*theta_num(i)-theta_num(i-1)-auxiliar*sin(theta_num(i));

end

plot(t,theta_num*180/pi,'-','Color',[0 0.5 0]);

legend('Solutia analitica','Solutia numerica');

figure(2); % simularea dinamica a oscilatiei

xp=lungime*sin(theta_num);

yp=-lungime*cos(theta_num);

x_analitic = lungime*sin(theta_an);
y_analitic = -lungime*cos(theta_an);

for i=1:1:N
    
    plot([0 xp(i)],[0 yp(i)],'-k',xp(i),yp(i),'ob'); hold on;
    plot([0 x_analitic(i)],[0 y_analitic(i)],'-r',x_analitic(i),y_analitic(i),'ob');
    hold off;

    xlabel('x/m'); ylabel('y/m'); grid;

    title('Simularea oscilatiei pendulului');

    axis equal;

    axis([-lungime lungime -lungime 0]); % stabileste limitele graficului [xmin xmax ymin ymax] 

    Film(i)=getframe; % "pozeaza" si stocheaza continutul ferestrei grafice

    % sau getframe(); la unele versiuni;
    

end

figure(3); % in timp real

movie(Film,1,N/tf);