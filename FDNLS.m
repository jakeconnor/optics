
%---Specify imput parameters
clear all;
number = 1; %input('Number of waves = ');
distance= input('Enter fiber length (in units of L_D) = ');
time = input('Time to run simulation = ');
% sigma = input('Sigma (for polarization) = ');
for j=1:number
    disp(sprintf('Wave %d' , j))
    
    beta1{j} = input('Beta 1 = '); % ?
    beta2{j} = input('dispersion: 1 for normal, -1 for anomalous = '); % ?
    N(j) = input('Nonlinear parameter N = '); %Solition order?
    mshape{j} = input('m = 0 for sech, m > 0 for super-Gaussian = '); % decided on input wave
    loss{j} = input('alpha, positive for loss, negative for gain = ');%maybe the same?
    %     ss{j} = input('self steepening factor = ');
end
chirp = 0; % input pulse chirp (default value)

%---set simulation parameters
nt = 1024; Tmax = 32; % FFT points and window size for graphs
step_num_t = round (50 * time * max(N).*2); % No. of z steps
step_num_z = round (200 * distance * max(N).*2);

deltat = time/step_num_t; % step size in z
deltaz= distance/step_num_z;

g1=zeros(step_num_z);
for n=1:step_num_z
    g1(n,n)=loss{j}/2;
    if n-1>0
        g1(n,n-1)=-1/(deltaz);
    end
    if n+1<=step_num_z
        g1(n,n+1)=1/(deltaz);
    end
    
end

zero=zeros(step_num_z);

c1= beta1{j}.*speye(step_num_z);
c3=-speye(step_num_z);
c4=-0.5i.*beta2{j}.*speye(step_num_z);
capmat=[c1 zero; c3 zero];

g1(1,:)=0;
g1(1,1)=1;
g1(step_num_z,:)=0;
g1(step_num_z,step_num_z)=1;

g3=zero;
g3(1,1)=1;
g4=eye(step_num_z);
indmat=sparse([zero zero; g3 zero]);
% z = (-step_num_z/2:step_num_z/2-1) * deltaz; % space grid
z2=-4*distance:8*distance/299:4*distance; %input wave, hopefully
% space grid
%---Input Field profile
for j=1:number %%
    uu{j}=zeros(1,step_num_z);
    uu{j}(1:300)=sech(z2);
    temp{j}(1,:) = uu{j}';
    X=zeros(2*step_num_z,step_num_t);
    X(:,1)=[uu{j}';zeros(step_num_z,1)];
    %plot(z,uu{j})
end
% ---plot inpute pulse shape and spectrum

%temp= fftshift(ifft(uu)).*(nt*dtau)/sqrt(2*pi);%spectrum

figure;
for n=1:step_num_t-1
    %     sumwaves=0; %for multiple waves in one fiber
    %     for j=1:number
    %         sumwaves=sumwaves+sigma.*(abs(uu{j}).^2 );.
    %     end
    for j=1:number
        
        if n==1 
            du_dt=(temp{j}(n,:)-temp{j}(n,:))/deltat;
            X(:,n+1)=vertcat(temp{j}(n,:)',du_dt');
        else
            %             temp{j}(n,:)=X(1:step_num_z,n)';
%             X(1,n)=0; X(step_num_z,:)=0;
            X(:,n+1)=((capmat/deltat)+indmat)\((capmat*(X(:,n-1))/deltat));
        end 
        
        %          du_dz(n,m)=(temp{j}(n,m+1)-temp{j}(n,m-1))./(2*deltaz);
        %          du_dt=(temp{j}(n,:)-temp{j}(n-1,:))/deltat;
        %          du2_dt2(n,m)=(du_dt(n,m)-du_dt(n-1,m)/deltat;
        
        
    end
    % subplot(2,3,1);
    plot(X(1:step_num_z,n))
    %  title('Wave (A) in space');
    %  hold on
    % subplot(2,3,2);
    %   mesh(X(1:step_num_z-1,:))%to plot over time
    %  title('Wave (A) in space & time');
    %  subplot(2,3,3);
    %  plot(X(1+step_num_z:2*step_num_z,n))
    %  title('dA/dt in space');
    % subplot(2,3,4);
    %  mesh(abs(X(1+step_num_z:2*step_num_z,:)))%to plot over time
    %  title('dA/dt in space & time');
    %  subplot(2,3,5);
    
    
    pause(30/step_num_t)
end

