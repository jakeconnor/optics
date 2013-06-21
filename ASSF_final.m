% This code solves the NLS equation with the split-step method
% idu/ dz - sgn(beta2) /2 d^2u/d(tau) ^2+N ^2*|u| ^2*|u| ^2u=0
% Written by Govind P. Agrawal in Marcg 2005 for the NLFO book
%---Specify imput parameters
clear all;
number = input('Number of waves = ');
distance= input('Enter fiber length (in units of L_D) = ');
sigma = input('Sigma (for polarization) = ');
for j=1:number
    disp(sprintf('Wave %d' , j))
    % entire legth, graphed at the end
    beta2{j} = input('dispersion: 1 for normal, -1 for anomalous = '); % ?
    N(j) = input('Nonlinear parameter N = '); %Solition order?
    mshape{j} = input('m = 0 for sech, m > 0 for super-Gaussian = '); % decided on input wave
    alpha{j} = input('alpha, positive for loss, negative for gain = ');%maybe the same?
    ss{j} = input('self steepening factor = ');
end
chirp = 0; % input pulse chirp (default value)
%---set simulation parameters
nt = 1024; Tmax = 32; % FFT points and window size for graphs
step_num = round (20 * distance * max(N).*2); % No. of z steps
deltaz = distance/step_num; % step size in z
dtau = (2 * Tmax)/nt; % step size in tau (radians from z)
%---tau and omega arrays
tau = (-nt/2:nt/2-1) * dtau; % temporal grid
omega = (pi/Tmax) * [(0:nt/2-1) (-nt/2:-1)]; % frequency grid
%---Input Field profile
for j=1:number
    if mshape{j}==0 % soliton sech shape creationg of input wave
        uu{j} = sech(tau).* exp(-0.5i* chirp* tau.^2);
    else % super-Gaussian only seems to be real numbers
        uu{j} = exp(-0.5* (1+1i*chirp).*tau.^(2*mshape{j}));
    end
end
% ---plot inpute pulse shape and spectrum

%temp= fftshift(ifft(uu)).*(nt*dtau)/sqrt(2*pi);%spectrum
for j=1:number
    %---Store dispersive phase shifts to speedup code
    dispersion{j} = exp(0.5i*beta2{j}*omega.^2*deltaz-0.5*alpha{j}*deltaz); %phase factor, includes alpha
    hhz{j}= 1i*N(j)^2 *deltaz; %nonlinear phase factor
    temp{j} = uu{j}.*exp(abs(uu{j}).^2*hhz{j}/2); % note hhz/2 first step, hhz (non-linear) decay
    
end
for n=0:step_num
    sumwaves=0;
    for j=1:number
        sumwaves=sumwaves+sigma.*(abs(uu{j}).^2 );
    end
    for j=1:number
        
        duu_dt{j}=ifft(1i.*omega.*temp{j});
        duu_dt_star{j}=ifft(1i.*omega.*fft(conj(temp{j})));
        selfsteepen=exp((2.*conj(temp{j}).*duu_dt{j} + uu{j}.*duu_dt_star{j}).*deltaz.*ss{j});
        
        f_temp{j} = ifft(temp{j}).*dispersion{j}; % other step, using dispersion decay
        
    end
    for j=number:-1:1
        uu{j} = fft(f_temp{j}); %amplitude updates the wave at that step
        temp{j} = uu{j}.*exp((sumwaves-(sigma-1).*abs(uu{j}).^2).*hhz{j}).* selfsteepen;
        f_temp{j} = ifft(temp{j});
 
        %% to plot at every step
        %             figure;
        %             subplot(2,1,1)
        %             plot(tau, abs(uu).^2, '-k')
        %             axis([-10 10 0 y1]);
        %             xlabel('Normalized Time');
        %             ylabel('Normalized Power');
        %             subplot(2,1,2)
        %             plot(fftshift(omega)/(2*pi), abs(temp_1plot).^2,'-k')
        %             axis([-2 2 0 y2]);
        %             xlabel('Normalized Frequency');
        %             ylabel('Spectral Power');
        %logs waveform over space
        tempplot{j} = fftshift(ifft(uu{j})).*(nt*dtau)/sqrt(2*pi);
        wavespace{j}(n+1,:)= abs(uu{j}).^2;
        wavefreq{j}(n+1,:)=abs(tempplot{j}).^2;
    end
    
end
for j=1:number
    figure;
    uu{j} = temp{j}.*exp(-abs(uu{j}).^2*hhz{j}/2); % Final field
    wavespace{j}(step_num+1,:)= abs(uu{j}).^2;
    wavefreq{j}(step_num+1,:)=abs(tempplot{j}).^2;
    tempplot{j} = fftshift(ifft(uu{j})).*(nt*dtau)/sqrt(2*pi); % plotted spectrum
    subplot(2,1,1)
    mesh(deltaz.*tau,deltaz.*(0:step_num),deltaz.*wavespace{j})
    view(37.5,30);
    % axis([-10 10 0 distance 0 y1])
    xlabel('T/T0');
    ylabel('Z/Ld');
    zlabel('Intesnsity');
    subplot(2,1,2)
    mesh(deltaz.*fftshift(omega)/(2*pi),deltaz.*(0:step_num), deltaz.*wavefreq{j})
    xlabel('(v-v0)T0');
    ylabel('Z/Ld');
    zlabel('Intensity');
    view(37.5,30);
    % axis([-2 2 0 distance 0 y2])
end
