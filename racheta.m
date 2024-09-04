clc; clear; close all;


% Selectori pentru grafice (1-on; 0-off):
Gv = 1; % componentele vitezei ca functii de timp
Gc = 1; % legile de miscare
Gt = 1; % traiectoria (curba balistica)
Gd = 1; % reprezentarea dinamica
if Gd==1, Gv=0; Gc=0; Gt=0; end

%Parametrii fizici:
g = 9.81;
m0 = 194; %masa initiala totala a rachetei
D = 0.18; %diametrul rachetei
r = 0.09; %raza rachetei
mc = 0.72 * m0; %masa combustibil
mr = m0 - mc; %masa racheta
rho0 = 1.22; %densitatea aerului (kg/m^3)
eta = 1.81 * 1e-5; %coef. de vascozitate (Pa * s)

%Conditii initiale:
%Precizare: pt. o mai buna observare a fenomenului, am micsorat tau de la
%57 la o secunda
v0 = 16; %m/s
alpha0 = 53; %unghiul de lansare (grade)
tau = 1; %timp de ardere (s)
u0 = 3880; %viteza combustibilului (față de corpul rachetei) (m/s)
q = mc / tau; %debit masic

%definirea intervalului de timp de interes
t0 = 0; tf = 300;
N = 10000;
time = linspace(t0, tf, N); dt = time(2) - time(1);
t_cadere = 300;
pas_cadere = N;
pas_coborare = N;

%prealocare si valori de inceput:
position = zeros(2, N);
velocity = zeros(2, N);
acceleration = zeros(2, N);
mass = zeros(1, N);
Fr = zeros(2, N);
G = zeros(2, N);
Ft = zeros(2, N);
R = zeros(2, N);
dissipated_energy = zeros(2, N);

position(:, 1) = [0; 0];
velocity(:, 1) = [v0 * cosd(alpha0); v0 * sind(alpha0)];
alpha = alpha0;
mass(1) = mr + mc;


i = 1;
while i < N
    %calculam fortele
    v = norm(velocity(:, i));
    Fr(:, i) = -6.54 * eta * D * velocity(:, i) - 0.64 * rho0 * D^2 * v * velocity(:, i);
    G(:, i) = [0; -mass(i) * g];
    
    % daca nu mai exista combustibil, tractiunea inceteaza
    if mass(i) > mr
        Ft(:, i) = [q * u0 * cosd(alpha); q * u0 * sind(alpha)];
    else
        Ft(:, i) = [0; 0];
    end
    R(:, i) = Fr(:, i) + G(:, i) + Ft(:, i);
    
    %actualizam pozitia, viteza si unghiul
    acceleration(:, i) = R(:, i) / mass(i);
    position(:, i+1) = position(:, i) + velocity(:, i) * dt;
    velocity(:, i+1) = velocity(:, i) + acceleration(:, i) * dt;
    alpha = atan(velocity(2, i+1) / velocity(1, i+1)); %arctan(vy/vx);
    
    %vb = viteza teoretica, daca nu am avea frecare
    vb = velocity(:, i+1) - Fr(:, i) ./ mass(i) .* dt;
    va = velocity(:, i+1);
    dissipated_energy(:, i) = 1/2 * mass(i) * (vb.^2 - va.^2);
    
    mass(i+1) = mass(i) - q * dt;
    if mass(i + 1) < mr
        mass(i + 1) = mr;
    end
    
    if position(2, i) < 0
        if pas_cadere == N
            pas_cadere = i;
        end
    end
    if velocity(2, i) < 0
        if pas_coborare == N
            pas_coborare = i;
        end
    end
    
    
    i = i + 1;
end

x = position(1, 1:pas_cadere);
y = position(2, 1:pas_cadere);
vx = velocity(1, 1:pas_cadere);
vy = velocity(2, 1:pas_cadere);
time = time(1:pas_cadere);

if Gv==1 % legile de viteza
    figure(1);
    plot(time,vx,'-r',time,vy,'-b');
    xlabel('t(s)'); ylabel('v(m/s)'); grid;
    title('Componentele vitezei ca functii de timp');
    legend('vx','vy');
end
if Gc==1 % legile de miscare
    figure(2);
    plot(time,x/1e3,'-r',time,y/1e3,'-b');
    xlabel('t(s)'); ylabel('coord(km)'); grid;
    title('Coordonatele ca functii de timp');
    legend('x','y','Location','northwest');
end
if Gt==1 %traiectoria
   figure(3);
   plot(x/1e3, y/1e3);
 
   xlabel('x(km)'); ylabel('y(km)'); grid;
   title('Problema rachetei');
   %axis equal; axis tight;
end

%Afisarea unor valori de interes:
timp_urcare = pas_coborare * dt;
timp_coborare = (N - pas_coborare) * dt;
b = x(pas_cadere);
h = max(y);

%suma pe x si pe y a energiei, va da N coloane
total_dissipated_energy = sum(dissipated_energy);
%se aduna cele N coloane
total_dissipated_energy = sum(total_dissipated_energy);

afis=['Timpul total de miscare: ', num2str(pas_cadere * dt), ' s']; disp(afis);
afis=['Bataia rachetei: ', num2str(b), ' m']; disp(afis);
afis=['Altitudinea maxima: ', num2str(h), ' m']; disp(afis);
afis=['Timp urcare: ', num2str(timp_urcare), ' s']; disp(afis);
afis=['Timp coborare: ', num2str(timp_coborare), ' s']; disp(afis);
afis=['Energia disipata: ', num2str(total_dissipated_energy), ' J']; disp(afis); %gresit


    

if Gd == 1 %simularea dinamica
    figure(4);
    set(4, 'Position', [100 100 800 600]);
    
    tic; simt =0; %porneste cronometrul si retine timpul initial
    while simt < (pas_cadere * dt)
       plot(x/1e3, y, '-c'); hold on;
       xlim([0 10]); ylim([0 100]);
       xlabel('x(km)'); ylabel('y(m)'); grid;
       title('Simularea miscarii');
       %axis equal; axis tight;
       index=abs(time-simt)==min(abs(time-simt));
       plot(x(index)/1e3,y(index),'.b','MarkerSize',10); hold off
        text(b/2/1e3,h/3,['vx=',num2str(round(round(vx(index)))),' m/s']);
        text((b/2-b/5)/1e3,h/3,['t=',num2str(round(time(index))),' s']);
        text((b/2+b/5)/1e3,h/3,['vy=',num2str(round(round(vy(index)))),' m/s']);
       pause(1e-3);
       simt = toc;
    end
end



%{
Ce îmbunătățiri ați aduce modelului propus, pentru o reprezentare mai fidelă a realității? 
1. Motorul rachetei poate avea un debit masic variabil. Mai ales la aproape
de finalul arderii, cand ramane foarte putin combustibil, debit masic sa
fie foarte mic
2. Densitatea atmosferica si coef. de vascozitate sa scada odata cu
altitudinea. Astfel vor exista pierderi mai mici din cauza rezistentei
aerului.
3. Forma specifica a rachetei, diferita de un cilindru.
4. Variabilitatea acceleratiei gravitationale: la altitudini asa mici,
variabilitatea ar fi totusi neglijabila.
5. Dinamica de separare: racheta sa aiba 2 sau 3 recipiente de combustibil, pe
care sa le arunce doar la golirea acestora.
7. Temperatura motorului poate afecta forta de tractiune

%}


