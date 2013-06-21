
%---Specify imput parameters
clear all;
number = 1; %input('Number of waves = ');
distance= input('Enter fiber length (in units of L_D) = ');
time = input('Time to run simulation = ');
% sigma = input('Sigma (for polarization) = ');
for j=1:number
    disp(sprintf('Wave %d' , j))
    
    beta1{j} = input('Beta 1 = '); % ?
    beta2{j} = input('dispersion: 1 for normal, -1 for anomalou s = '); % ?
    N(j) = input('Nonlinear parameter N = '); %Solition order?
    mshape{j} = input('m = 0 for sech, m > 0 for super-Gaussian = '); % decided on input wave
    loss{j} = input('alpha, positive for loss, negative for gain = ');%maybe the same?
    %     ss{j} = input('self steepening factor = ');
end
chirp = 0; % input pulse chirp (default value)

%---set simulation parameters
nt = 1024; Tmax = 32; % FFT points and window size for graphs
step_num_t = round (500 * time * max(N).*2); % No. of z steps
deltat = time/step_num_t; % step size in z
deltaz=deltat/beta1{j};
r=deltat/(2*beta1{j}*deltaz);
step_num_z=distance/deltaz;
g1=zeros(step_num_z);
%gc1=zeros(step_num_z);
for n=1:step_num_z 
    %setting up an 8th order accurate central difference, overkill, but best chance.
    %     gc1(n,n)=loss{j}/2;
    %     if n+4 > step_num_z; gc1(n,n+4-step_num_z) = -1/(280*deltaz);
    %     else gc1(n,n+4)=-1/(280*deltaz); end
    %     if n+3 > step_num_z; gc1(n,n+3-step_num_z) = 4/(105*deltaz);
    %     else gc1(n,n+3)=4/(105*deltaz); end
    %     if n+2 > step_num_z; gc1(n,n+2-step_num_z) = -1/(5*deltaz);
    %     else gc1(n,n+2)=-1/(5*deltaz); end
    %     if n+1 > step_num_z; gc1(n,n+1-step_num_z) = 4/(5*deltaz);
    %     else gc1(n,n+1)=4/(5*deltaz); end
    %     if n-1 < 1; gc1(n, step_num_z+n-1)= -4/(5*deltaz);
    %     else gc1(n,n-1)= -4/(5*deltaz); end
    %     if n-2 < 1; gc1(n, step_num_z+n-2)= 1/(5*deltaz);
    %     else gc1(n,n-2)= 1/(5*deltaz); end
    %     if n-3 < 1; gc1(n, step_num_z+n-3)= -4/(105*deltaz);
    %     else gc1(n,n-3)= -4/(105*deltaz); end
    %     if n-4 < 1; gc1(n, step_num_z+n-4)= 1/(280*deltaz);
    %     else gc1(n,n-4)= 1/(280*deltaz); end
    %
   
    
    
    %6th order accurate upwinding
%     g1(n,n)=49/(20*deltaz);
%     if n-1 > 0; g1(n,n-1)= -6/deltaz;
%     else g1(n,step_num_z+n-1)= -6/deltaz; end
%     if n-2 > 0 ; g1(n,n-2) = 15/(2*deltaz);
%     else g1(n,step_num_z+n-2)=15/(2*deltaz); end
%     if n-3 > 0; g1(n,n-3)= -20/(3*deltaz);
%     else g1(n,step_num_z+n-3)=-20/(3*deltaz); end
%     if n-4 > 0; g1(n,n-4)=15/(4*deltaz);
%     else g1(n,step_num_z+n-4)=15/(4*deltaz); end
%     if n-5 > 0; g1(n,n-5)=-6/(5*deltaz);
%     else g1(n, step_num_z+n-5) = -6/(5*deltaz); end
%     if n-6 > 0; g1(n,n-6)=1/(6*deltaz);
%     else g1(n,step_num_z+n-6)=1/(6*deltaz); end
    
    %1st order accurate upwinding
    g1(n,n)=sign(beta1{j})/deltaz + loss{j}/2; %%for upwinding
    if n-sign(beta1{j})>0 & n-sign(beta1{j})<=step_num_z; %#ok<AND2>
        g1(n,n-sign(beta1{j}))=-sign(beta1{j})/deltaz;
    end
    
end
%gc1(1,step_num_z)=-1/deltaz;
%gc1(step_num_z,1)=1/deltaz;


g1(step_num_z,step_num_z)=1/deltaz;
g1(1,step_num_z)=-1/deltaz;
zero=zeros(step_num_z);
g4=eye(step_num_z);


c1=beta1{j}*speye(step_num_z);

c3=-speye(step_num_z);
c4=-0.5i.*beta2{j}.*speye(step_num_z);
c4=1i*speye(step_num_z);
capmat=[c1 zero; c3 c4];

% g1(1,:)=0;
% g1(1,1)=1;

g3=zero;
g3(1,1)=1;
indmat=sparse([g1 zero; g3 g4]);
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
burr=diag(ones(step_num_z,1))+diag((-r/2)*ones(step_num_z-1,1),-1)+diag((r/2)*ones(step_num_z-1,1),1);
nburr=diag(ones(step_num_z,1))+diag((-r/2)*ones(step_num_z-1,1),1)+diag((r/2)*ones(step_num_z-1,1),-1);
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
            X(:,n)=vertcat(temp{j}(n,:)',du_dt');
            X(:,n+1)=((capmat)+indmat*deltat)\((capmat*(X(:,n))));
        else
            %             du_dt= (burr\(nburr*X(1:step_num_z,n)) -
            %             X(1:step_num_z,n))/deltat; detection of where the derivative
            %             changes to/from zero.
            %             waveder = (abs(du_dt) >= eps);
            %             waveder= diff([1; waveder; 1]);
            %             firstzero = find(waveder < 0);
            %             lastzero = find(waveder > 0) -1;
            
            %             temp{j}(n,:)=X(1:step_num_z,n)';
            %             X(1,n)=0; X(step_num_z,:)=0;
            X(:,n+1)=((capmat)+indmat*deltat)\((capmat*(X(:,n))));
        end
        %          du_dz(n,m)=(temp{j}(n,m+1)-temp{j}(n,m-1))./(2*deltaz);
        %          du_dt=(temp{j}(n,:)-temp{j}(n-1,:))/deltat;
        %          du2_dt2(n,m)=(du_dt(n,m)-du_dt(n-1,m)/deltat;
        
        
    end
    if mod(n,round(step_num_z/100))==0
        % subplot(2,3,1);
        %     plot(X(1:step_num_z,n))
        %  title('Wave (A) in space');
        %  hold on
        % subplot(2,3,2);
        mesh(X(1:step_num_z-1,:))%to plot over time
        %  title('Wave (A) in space & time');
        %  subplot(2,3,3);
        %  plot(X(1+step_num_z:2*step_num_z,n))
        %  title('dA/dt in space');
        % subplot(2,3,4);
        %       mesh(abs(X(1+step_num_z:2*step_num_z,:)))%to plot over time
        %  title('dA/dt in space & time');
        %  subplot(2,3,5);
        
        
        pause(30/step_num_t)
    end
end

