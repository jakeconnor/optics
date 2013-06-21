% This code solves the NLS equation with the split-step method
% idu/ dz - sgn(beta2) /2 d^2u/d(tau) ^2+N ^2*|u| ^2*|u| ^2u=0
% Written by Govind P. Agrawal in Marcg 2005 for the NLFO book
%---Specify imput parameters
clear all;
distance = input('Enter fiber length (in units of L_D) ='); % entire legth, graphed at the end
beta2 = input('dispersion: 1 for normal, -1 for anomalous = '); % ?
N = input('Nonlinear parameter N = '); %Solition order?
mshape = input('m = 0 for sech, m > 0 for super-Gaussian = '); % decided on input wave
alpha = input('alpha, positive for loss, negative for gain = ');
ss = input('self steepening factor = ');
% y1=input('y1= ');
% y2=input('y2= ');
chirp0 = 0; % input pulse chirp (default value) ???????????????????????????????????????????????????
%---set simulation parameters
nt = 1024; Tmax = 32; % FFT points and window size for graphs
step_num = round (20 * distance * N^2); % No. of z steps
deltaz = distance/step_num; % step size in z
dtau = (2 * Tmax)/nt; % step size in tau (radians from z)
%---tau and omega arrays
tau = (-nt/2:nt/2-1) * dtau; % temporal grid
omega = (pi/Tmax) * [(0:nt/2-1) (-nt/2:-1)]; % frequency grid
%---Input Field profile
if mshape==0 % soliton sech shape creationg of input wave
    uu = sech(tau).* exp(-0.5i* chirp0* tau.^2);
else % super-Gaussian only seems to be real numbers
    uu = exp(-0.5* (1+1i*chirp0).*tau.^(2*mshape));
end

% ---plot inpute pulse shape and spectrum

temp= fftshift(ifft(uu)).*(nt*dtau)/sqrt(2*pi);%spectrum

%---Store dispersive phase shifts to speedup code
dispersion = exp(0.5i*beta2*omega.^2*deltaz-0.5*alpha*deltaz); %phase factor, includes alpha
hhz= 1i*N^2 *deltaz; %nonlinear phase factor
%********[beginning of MAIN Loop] ********

temp = uu.*exp(abs(uu).^2*hhz/2); % note hhz/2 first step, hhz (non-linear) decay

for n=0:step_num

    %todo: split-step slope steepening (after making it work)
    
    f_temp = ifft(temp).*dispersion; % other step, using dispersion decay
    uu = fft(f_temp); %amplitude updates the wave at that step
    temp = uu.*exp(abs(uu).^2.*hhz);
     f_temp = ifft(temp);
    
     % self steepening
    duu_dt=ifft(1i.*omega.*temp);
    duu_dt_star=ifft(1i.*omega.*fft(conj(temp)));
    selfsteepen=exp((2.*conj(temp).*duu_dt + uu.*duu_dt_star).*deltaz.*ss);
    temp= uu.* selfsteepen;
    uu=temp;
        
 
    %% to plot at every step
    %             figure;
    %             subplot(2,1,1)
    %             plot(tau, abs(uu).^2, '-k')
    %             axis([-10 10 0 y1]);
    %             xlabel('Normalized Time');
    %             ylabel('Normalized Power');
    %             subplot(2,1,2)
    %             plot(fftshift(omega)/(2*pi), abs(tempplot).^2,'-k')
    %             axis([-2 2 0 y2]);
    %             xlabel('Normalized Frequency');
    %             ylabel('Spectral Power');
    %logs waveform over space
    tempplot = fftshift(ifft(uu)).*(nt*dtau)/sqrt(2*pi);
    wavespace(n+1,:)= abs(uu).^2;
    wavefreq(n+1,:)=abs(tempplot).^2;
    
end
figure;
uu = temp.*exp(-abs(uu).^2*hhz/2); % Final field
wavespace(step_num+1,:)= abs(uu).^2;
wavefreq(step_num+1,:)=abs(tempplot).^2;
tempplot = fftshift(ifft(uu)).*(nt*dtau)/sqrt(2*pi); % plotted spectrum
subplot(2,1,1)
mesh(deltaz.*tau,deltaz.*(0:step_num),deltaz.*wavespace)
view(37.5,30);
% axis([-10 10 0 distance 0 y1])
xlabel('T/T0');
ylabel('Z/Ld');
zlabel('Intesnsity');
subplot(2,1,2)
mesh(deltaz.*fftshift(omega)/(2*pi),deltaz.*(0:step_num), deltaz.*wavefreq)
xlabel('(v-v0)T0');
ylabel('Z/Ld');
zlabel('Intensity');
view(37.5,30);
% axis([-2 2 0 distance 0 y2])

