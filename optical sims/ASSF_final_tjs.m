% This code solves the NLS equation with the split-step method
% idu/ dz - sgn(beta2) /2 d^2u/d(tau) ^2+N ^2*|u| ^2*|u| ^2u=0
% Modified from version written by Govind P. Agrawal in Marcg 2005 for the NLFO book

%  idu/ dz - sgn(beta2) /2 d^2u/d(tau) ^2+N ^2*|u|^2 * u=0
% SS
%     du/dz = i sgn(beta2)/2 d^2u/dtau^2 + i N^2 *|u|^2 u  
%     direct integration in z:  u(x0+dz) = u(x0) exp(O dz) 
%     1)  In freq domain: O = FT[i sgn(beta2)/2 d^2u/dtau^2] = i sgn(beta2)/2 omega^2
%         lindisp = exp(i sgn(beta2)/2 omega^2 dz)
%         u(f,x0+dz) = u(f,x0) * lindisp
%     2)  In time domain 0 = N ^2*|u| ^2 and u(t) = fff(u(f)
%         u(t,x0+dz) = u(t,x0) exp(i N ^2 *dz |u| ^2)
%
% FD:
%    du/dz = i sgn(beta2)/2 d^2u/dtau^2 + i N^2 *|u|^2 u  
% 1) disp
%         i sgn(beta2)/2 d^2u/dtau^2 = i sgn(beta2)/2 G u
%                    du/dz = i sgn(beta2)/2 G u
%                    du(n)/dz  = alpha * (u(n+1) - 2 u(n) + u(n-1))
%                    (u - up)/dz = i sgn(beta2)/2 G u
%                    (I - dz i sgn(beta2)/2 G) u = up
%                    u = (1 - dz i sgn(beta2)/2 G)\up
%   2) nonlinear
%       u(t,x0+dz) = u(t,x0) exp(i N ^2 *dz |u| ^2)

%FFFFFFIIIIXXXXXXXX
%Does not work for non zero sigma value if only one wave
%sigma set to 0 if only one wave. should recheck math to see if multiwave
%approach is valid.
clear;
close all;

set(0,'DefaultFigureWindowStyle','docked')

UserInput = 0;
% UserInput = input('User Input? 1 or 0');

% T = t - z/Vg
% tau = T/T0 -- T0 is pulse width
% tau = (t-z/Vg)/T0;
% Ld = T0^2/abs(Beta2);
% Lnl = 1/(gamma P0) -- P0 is pulse peak power
% 
% beta2 ~ 20 ps^2/km
% gamma ~ 1/(W Km)
% T0 > 100 ps
% P0 ~ 1 W
% L ~ 50Km
% alpha ~ 0.2dB/km
% beta1 = 1/Vg = 1/(c/ng) = ng/c  ~1.47/3.0e8 s/m
%

if (UserInput)

    distance= 0+ input('Enter fiber length (in units of L_D) = ');% based off of first wave; entire length, graphed at the end
    sigma=0;
    
    % specify input parameters
    beta1 = 1;
    beta2 =  0 + input('dispersion: 1 for normal, -1 for anomalous = '); % dispersion of the group velocity (responsible for pulse broadening)
    N = 0 + input('Nonlinear parameter N = '); %Solition order
    mshape = 0;%input('m = 0 for sech, m > 0 for super-Gaussian = '); % decided on input wave
    alpha = 0+ input('alpha, positive for loss, negative for gain = ');% loss or gain factor (same for each wave?)
    ss = 0+ input('self steepening factor = ');
    
    chirp0 = 0; % input pulse chirp (default value)
else
    beta2 = 20; % ps^2/Km
    T0 = 100;   % ps
    L = 50;     % km
    
    Ld = T0^2/abs(beta2);    % Km
    Ln1 = 1;                 % 1/(W Km)
    beta1 = 1.47/3.0e8*1000; % s/km
    beta1 = beta1/1e-12;     % ps/km
    
    Pl = T0/beta1;  % Km
    Plm = Pl*1000;
    
    % normalized Z = z/Ld  tau = T/T0 T = t - beta1 z or t = T + beta1 z
    %                                                      = T + beta1 Ld*Z        
    % 
    %                 t = tau*T0 + beta1*Ld*Z = tau*T0 + Toff;   
    
    distance=15;% 
    sigma=0;
    
    % specify input parameters (Normalized)
    beta2 = -0; % this captures actual beta2 value in Ld normalization
    N = 2.0;
    mshape = 1;%input('m = 0 for sech, m > 0 for super-Gaussian = '); % decided on input wave
    alpha = 0;
    ss = 0;

    chirp0 = 0; % input pulse chirp (default value)
end

% set simulation parameters
nt = 1024; Tmax = 32; % FFT points and window size for graphs
% step_num = round (20 * distance * N{1}^2); % No. of z steps
step_num = 500;
deltaz = distance/step_num; % step size in z
dtau = (2 * Tmax)/nt; % step size in tau (radians from z)

% tau and omega arrays
tau = (-nt/2:nt/2-1) * dtau; % temporal grid
omega = (pi/Tmax) * [(0:nt/2-1) (-nt/2:-1)]; % frequency grid

if mshape==0 % soliton sech shape creationg of input wave
    uu = sech(tau).* exp(-0.5i* chirp0* tau.^2);
else % super-Gaussian (only seems to be real numbers)
    uu = exp(-0.5* (1+1i*chirp0).*tau.^(2*mshape));
end
uu0 = uu;

% store dispersive phase shifts to speedup code
lindispersion = exp(0.5i*beta2*omega.^2*deltaz-0.5*alpha*deltaz); %phase factor, includes alpha
hhz= 1i*N^2 *deltaz; %nonlinear phase factor
temp = uu.*exp(abs(uu).^2*hhz/2); % note hhz/2 first step, hhz (non-linear) decay
tempFD = temp;

G = zeros(nt);
D = eye(nt);
for j=1:nt
    if j == 1
        G(j,nt) = 1;
        G(j,j) = - 2;
        G(j,j+1) = 1;
    elseif j == nt
            G(j,j-1) = 1;
            G(j,j) = -2;
            G(j,1) = 1;
    else
        G(j,j-1) = 1;
        G(j,j) = -2;
        G(j,j+1) = 1;
    end
end
G = G/(dtau*dtau);

LG =  (D - deltaz*0.5i*beta2*G);
[l,u] = lu(LG);

% uu2 = fft(ifft(temp));
% 
% duudt = -fft(ifft(temp).*omega.^2);
% duuFDdt = G*tempFD';


for n=1:step_num
    %calculating derivatives of wave in frequency and time domain
    
    
%     plot(abs(temp).^2);
%     hold
%     plot(abs(tempFD).^2,'r+');

    % Line Disp 
    f_temp = ifft(temp).*lindispersion; % step using dispersion decay (goes from wave 1-number; taking advantage of the way this loop works)
    uu = fft(f_temp); %amplitude updates the wave at that step (after step at line 73)
    % 
    uuFD = (u\(l\(tempFD')))';
    
%     plot(abs(uu).^2,'bo');
%     plot(abs(uuFD).^2,'g*');
%     
%     drawnow();
%     pause(0.1)
    
    %non linear dispersion
    temp = uu.*exp((abs(uu).^2).*hhz); 
    tempFD = uuFD.*exp((abs(uuFD).^2).*hhz); 
    
    f_temp = ifft(temp);% stores the frequency domain; used later to store final wave
        
    tempplot = fftshift(ifft(uu)).*(nt*dtau)/sqrt(2*pi);% plotted spectrum
    tempplotFD = fftshift(ifft(uuFD)).*(nt*dtau)/sqrt(2*pi);% plotted spectrum
    
    % logs data at every step (except last)so it may be plotted later
    wavespace(n+1,:)= abs(uu).^2;
    wavespaceFD(n+1,:)= abs(uuFD).^2;
    wavefreq(n+1,:)=abs(tempplot).^2;
    wavefreqFD(n+1,:)=abs(tempplotFD).^2;
    
end

figure;% to create new plot

Za = deltaz.*(0:step_num)
Toff = Za*beta1*Ld; %ps

t = 512*dtau*T0; 
TT = Toff + t;

% plotting wave envelope
subplot(2,1,1)
mesh(deltaz.*tau,deltaz.*(0:step_num),deltaz.*wavespace)
xlabel('T/T0');
ylabel('Z/Ld');
zlabel('Intensity');
view(37.5,30);

% plotting frequency envelope
subplot(2,1,2)
mesh(deltaz.*fftshift(omega)/(2*pi),deltaz.*(0:step_num), deltaz.*wavefreq)
xlabel('(v-v0)T0');
ylabel('Z/Ld');
zlabel('Intensity');
view(37.5,30);
