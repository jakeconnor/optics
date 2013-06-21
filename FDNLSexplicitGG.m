
%---Specify imput parameters
clear;
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
    
    beta1{j} = input('Beta 1 = '); % ?
    beta2{j} = input('dispersion: 1 for normal, -1 for anomalous = '); % ?
    N(j) = input('Nonlinear parameter N = '); %Solition order?
    mshape{j} = input('m = 0 for sech, m > 0 for super-Gaussian = '); % decided on input wave
    loss{j} = input('alpha, positive for loss, negative for gain = ');%maybe the same?
    gamma = input('gamma = ');
    %     ss{j} = input('self steepening factor = ');
end
chirp = 0; % input pulse chirp (default value)
r=deltat/(2*beta1*deltaz)
%---set simulation parameters
nt = 1024; Tmax = 32; % FFT points and window size for graphs
step_num_t = round (2500 * time * max(N).*2); % No. of z steps
step_num_z = round (100 * distance * max(N).*2);

deltat = time/step_num_t; % step size in z
deltaz= distance/step_num_z;

g1=zeros(step_num_z);
for n=1:step_num_z
    g1(n,n)=-sign(beta1{j})/deltaz-loss{j}/2; %%for upwinding
    if n+sign(beta1{j})>0 && n+sign(beta1{j})<=step_num_z;
        g1(n,n+sign(beta1{j}))=sign(beta1{j})/deltaz;
    end
end

gg1=zeros(step_num_z);
for n=1:step_num_z
    gg1(n,n)=-(sign(beta1{j})*2)/(deltaz.^2); %%for upwinding
    if n<step_num_z
        gg1(n,n+1)=sign(beta1{j})/(deltaz.^2);
    end
    if n>1
        gg1(n,n-1)=sign(beta1{j})/(deltaz.^2);
    end
end
zero=zeros(step_num_z);
gg1(step_num_z,1)=1/(deltaz.^2);
gg1(1,step_num_z)=1/(deltaz.^2);
gg3=speye(step_num_z);
GG=[gg1 zero; zero gg3];

% g1(step_num_z,step_num_z)=-1/deltaz;
g1(step_num_z,1)=1/deltaz;

g4=eye(step_num_z);


c1=beta1{j}*speye(step_num_z);

c3=-speye(step_num_z);
c4=-0.5i.*beta2{j}.*speye(step_num_z);
c4=i*speye(step_num_z);
capmat=[c1 zero; c3 c4];


% g1(1,:)=0;
% g1(1,1)=1;

g3=zero;
g3(1,1)=1;
indmat=sparse([g1 zero; g3 g4]);
% z = (-step_num_z/2:step_num_z/2-1) * deltaz; % space grid
z2=-4*distance:8*distance/299:4*distance; %input wave, hopefully
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
        
        if n==1 %need to fix X so is actually a column vector, i think?
            du_dt=(temp{j}(n,:)-temp{j}(n,:))/deltat;
            X(:,n+1)=vertcat(temp{j}(n,:)',du_dt');
        else
            %             X(end-2:end,n)=0;
            %             temp{j}(n,:)=[X(1:step_num_z,n)';zeros(step_num_z,1)];
            %             du_dz=(temp{j}(n,[2:step_num_z,1])-temp{j}(n,[step_num_z,1:step_num_z-1]))./(2*deltaz);
            X(:,n+1)=X(:,n)-(deltat.*((indmat*X(:,n))+(gamma*(GG*X(:,n)))))/beta1{j};
            if n==500
                gamma=0;
            end
        end
        
        
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
    %  mesh(abs(X(1+step_num_z:2*step_num_z,:)))%to plot over time
    %  title('dA/dt in space & time');
    %  subplot(2,3,5);
    
        if mod(n,16)==0
    mesh(X(1:step_num_z,:))
    pause(30/step_num_t);
         end
end

