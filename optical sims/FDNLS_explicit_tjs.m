
%---Specify imput parameters
clear;

close all
set(0,'DefaultFigureWindowStyle','docked')

UserInput = input('User Input? 1 or 0');

if(UserInput)
    number = 1; %input('Number of waves = ');
    distance= input('Enter fiber length (in units of L_D) = ');
    time = input('Time to run simulation = ');
    % sigma = input('Sigma (for polarization) = ');
    
    % 2,2,.5,1,1,0,0 shows smooth until it returns to beinning
    % 2,2,.25,1,1,0,0 breaks down much faster, higher velocity means more
    % instability
    % 2,2,1,1,1,0,0 is very stable, but leaves noise bhind as it returns to the
    % beginning, which eventually breaks down. using more space steps or less time steps increases the
    % speed it breaks down also.
    
    
    
    for j=1:number
        disp(sprintf('Wave %d' , j))
        
        beta1{j} = -input('Beta 1 = '); % ?
        beta2{j} = input('dispersion: 1 for normal, -1 for anomalous = '); % ?
        N(j) = input('Nonlinear parameter N = '); %Solition order?
        mshape{j} = input('m = 0 for sech, m > 0 for super-Gaussian = '); % decided on input wave
        loss{j} = input('alpha, positive for loss, negative for gain = ');%maybe the same?
        %     ss{j} = input('self steepening factor = ');
    end
    %---set simulation parameters
    nt = 1024; Tmax = 32; % FFT points and window size for graphs
    MinDisp = 0;
    periodic = 0;
else
    number = 1;
    distance = 2;
    time = 20;
    j = 1;
    beta1{j} = 1;
    beta2{j} = -1e-5;
    N(j) = 0;
    loss{j} = 0e-1;
    %---set simulation parameters
    nt = 512; Tmax = 32; % FFT points and window size for graphs
    MinDisp = 1;
    periodic = 1;
end
chirp = 0; % input pulse chirp (default value)

if (MinDisp)
%     step_num_z = round (100 * distance * max(N).*2);
    step_num_z = 512;
    deltaz = distance/step_num_z;   
    deltat = deltaz*beta1{1};
    step_num_t = time/deltat;
else
    %deltat/deltaz * 1/beta1 ==1/2
    ratio=-beta1{j}/2;
    %  step_num_t = round (50000 * time * max(N).*2); % No. of z steps
    step_num_z = round (100 * distance * max(N).*2);
    deltaz= distance/step_num_z;
    step_num_t=step_num_z*ratio;
    deltat = time/step_num_t; % step size in z
    j = 1;
end


% beta1 dA/dt + dA/dz = i/2 beta2 dA^2/dt^2 + F(A)

% Gmat = dA/dz   


Gmat=zeros(step_num_z);
for n=1:step_num_z
    Gmat(n,n)=sign(beta1{j})/deltaz + loss{j}/2; %%for upwinding
    
    if n-sign(beta1{j})>0 & n-sign(beta1{j}) <= step_num_z;
        Gmat(n,n-sign(beta1{j}))=-sign(beta1{j})/deltaz;
    end
end


if (periodic)
    Gmat(1,step_num_z) = -sign(beta1{j})/deltaz;
end

% g1(1,:)=0;
% g1(1,1)=1;
% g1(1,step_num_z)=1/deltaz;

pulsewidth_n = .3*step_num_z;
% z = (-step_num_z/2:step_num_z/2-1) * deltaz; % space grid
z2=-4*distance:8*distance/(pulsewidth_n-1):4*distance; %input wave, hopefully
%---Input Field profile
for j=1:number %%
    uu{j}=zeros(1,step_num_z);
    uu{j}(1:pulsewidth_n)=sech(z2);
    temp{j}(1,:) = uu{j}';
    X=zeros(step_num_z,step_num_t);
    X(:,1)=[uu{j}'];
    Xp = X(:,1);
    Xpp = X(:,1);
    Xppp = X(:,1);
    
    dXdT = (X(:,1) - Xp)/deltat;
    dXdT2 = (2*X(:,1) - 5*Xp + 4*Xpp - Xppp)/(deltat*deltat);
    %plot(z,uu{j})
end

figure;
subplot(1,3,1), plot(abs(X(:,1)),'rx'), axis([0 step_num_z 0 1]);
% hold

subplot(1,3,2), plot(abs(dXdT(:)),'bx');
% hold

subplot(1,3,3), plot(abs(dXdT2(:)),'gx');
% hold


for n=1:step_num_t-1
    %     sumwaves=0; %for multiple waves in one fiber
    %     for j=1:number
    %         sumwaves=sumwaves+sigma.*(abs(uu{j}).^2 );.
    %     end
    for j=1:number
           X(:,n+1) = X(:,n) + ... 
               deltat/beta1{j}*(-Gmat*X(:,n) + i*0.5*beta2{j}*dXdT2);
        
           Xp = X(:,n);
           if (n > 2) Xpp = X(:,n-1); end
           if (n > 3) Xppp = X(:,n-2); end
           
           dXdT = (X(:,n+1) - Xp)/deltat;
           dXdT2 = (2*X(:,n+1) - 5*Xp +4*Xpp - Xppp)/(deltat*deltat);
        %          du_dt=(temp{j}(n,:)-temp{j}(n-1,:))/deltat;
        %          du2_dt2(n,m)=(du_dt(n,m)-du_dt(n-1,m)/deltat;
        
        
    end
    % subplot(2,3,1);
    
    %  title('Wave (A) in space');
    %  hold on
    % subplot(2,3,2);
    %   mesh(X(1:step_num_z-1,:))%to plot over time
    %  title('Wave (A) in space & time');
    %  subplot(2,3,3);
    %  plot(X(1+step_num_z:2*step_num_z,n))
    %  title('dA/dt in space');
    % subplot(2,3,4);
    %
    %  title('dA/dt in space & time');
    %  subplot(2,3,5);
    
     if mod(n,20)==0
        %          plot(X(:,n))
        %         mesh(abs(X(1+step_num_z:2*step_num_z,:)))%to plot over time
        
        subplot(1,3,1), plot(abs(X(:,n+1)),'r-x'), axis([0 step_num_z 0 1]);
        subplot(1,3,2), plot(abs(dXdT(:)),'b-x');
        subplot(1,3,3), plot(abs(dXdT2(:)),'g-x');
        
        drawnow;
        pause(0.000001);
    end
end

