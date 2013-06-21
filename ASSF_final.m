% This code solves the NLS equation with the split-step method
% idu/ dz - sgn(beta2) /2 d^2u/d(tau) ^2+N ^2*|u| ^2*|u| ^2u=0
% Modified from version written by Govind P. Agrawal in Marcg 2005 for the NLFO book


%FFFFFFIIIIXXXXXXXX
%Does not work for non zero sigma value if only one wave
%sigma set to 0 if only one wave. should recheck math to see if multiwave
%approach is valid.
clear all;

% accomodating for multiple wave input
number = 0 + input('Number of waves = ');
distance= 0+ input('Enter fiber length (in units of L_D) = ');% based off of first wave; entire length, graphed at the end
if number >1; sigma = 0+ input('Sigma (for polarization) = ');
else sigma=0; end

% looping to get correct no. of waves
for j=1:number
    disp(sprintf('Wave %d' , j))
    
    % specify input parameters
    beta2{j} =  0 + input('dispersion: 1 for normal, -1 for anomalous = '); % dispersion of the group velocity (responsible for pulse broadening)
    N{j} = 0 + input('Nonlinear parameter N = '); %Solition order
    mshape{j} = 0;%input('m = 0 for sech, m > 0 for super-Gaussian = '); % decided on input wave
    alpha{j} = 0+ input('alpha, positive for loss, negative for gain = ');% loss or gain factor (same for each wave?)
    ss{j} = 0+ input('self steepening factor = ');
end

chirp0 = 0; % input pulse chirp (default value)

% set simulation parameters
nt = 1024; Tmax = 32; % FFT points and window size for graphs
step_num = round (20 * distance * N{1}^2); % No. of z steps
deltaz = distance/step_num; % step size in z
dtau = (2 * Tmax)/nt; % step size in tau (radians from z)

% tau and omega arrays
tau = (-nt/2:nt/2-1) * dtau; % temporal grid
omega = (pi/Tmax) * [(0:nt/2-1) (-nt/2:-1)]; % frequency grid

% input field profile
for j=1:number
    if mshape{j}==0 % soliton sech shape creationg of input wave
        uu{j} = sech(tau).* exp(-0.5i* chirp0* tau.^2);
    else % super-Gaussian (only seems to be real numbers)
        uu{j} = exp(-0.5* (1+1i*chirp0).*tau.^(2*mshape{j}));
    end
end


for j=1:number
    % store dispersive phase shifts to speedup code
    lindispersion{j} = exp(0.5i*beta2{j}*omega.^2*deltaz-0.5*alpha{j}*deltaz); %phase factor, includes alpha
    hhz{j}= 1i*N{j}^2 *deltaz; %nonlinear phase factor
    temp{j} = uu{j}.*exp(abs(uu{j}).^2*hhz{j}/2); % note hhz/2 first step, hhz (non-linear) decay
    
end


for n=1:step_num
    sumwaves=0; % initializing summation of square of every wave
    
    
    for j=1:number
        %summation of the square of the waves
        sumwaves=sumwaves+sigma.*(abs(uu{j}).^2 );
    end
    
    %calculating derivatives of wave in frequency and time domain
    for j=1:number
        duu_dt{j}=ifft(1i.*omega.*temp{j}); % frequency domain(as ifft and fft are symmetrical)
        duu_dt_star{j}=ifft(1i.*omega.*fft(conj(temp{j}))); % time domain
        
        % calculating self steepening
        selfsteepen=exp((2.*conj(temp{j}).*duu_dt{j} + uu{j}.*duu_dt_star{j}).*deltaz.*ss{j});
        
        f_temp{j} = ifft(temp{j}).*lindispersion{j}; % step using dispersion decay (goes from wave 1-number; taking advantage of the way this loop works)
        
    end
    
    %looping from last wave to first
    for j=number:-1:1
        uu{j} = fft(f_temp{j}); %amplitude updates the wave at that step (after step at line 73)
        
        % subtracting intial wave from sumwaves as in the math all waves
        % but the first contribute sigma their squares and first wave
        % contributes 1
        % (recall the typical case when sigma is 2; every wave but first contributes twice their squares)
        temp{j} = uu{j}.*exp((sumwaves-(sigma-1).*abs(uu{j}).^2).*hhz{j}).* selfsteepen; % applying self-steepening and non linear dispersion
        
        f_temp{j} = ifft(temp{j});% stores the frequency domain; used later to store final wave
        
        
        tempplot{j} = fftshift(ifft(uu{j})).*(nt*dtau)/sqrt(2*pi);% plotted spectrum
        
        % logs data at every step (except last)so it may be plotted later
        wavespace{j}(n+1,:)= abs(uu{j}).^2;
        wavefreq{j}(n+1,:)=abs(tempplot{j}).^2;
    end
    
end
for j=1:number
    figure;% to create new plot
    
    % plotting wave envelope
    subplot(2,1,1)
    mesh(deltaz.*tau,deltaz.*(0:step_num),deltaz.*wavespace{j})
    xlabel('T/T0');
    ylabel('Z/Ld');
    zlabel('Intensity');
    view(37.5,30);
    
    % plotting frequency envelope
    subplot(2,1,2)
    mesh(deltaz.*fftshift(omega)/(2*pi),deltaz.*(0:step_num), deltaz.*wavefreq{j})
    xlabel('(v-v0)T0');
    ylabel('Z/Ld');
    zlabel('Intensity');
    view(37.5,30);
    
end
