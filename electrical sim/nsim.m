%NSIM (not yet fully tested) Simulates Basic Electrical Circuits with netlist option 
%Naming Convention:
% V = Voltage Source
% R = Resistor
% C = Capacitor
% L = inductors
% D = Voltage Controlled and Delayed Source
% E = Voltage Controlled Voltage Source
% X = Optical Coupler
% P = Phase Shift


%read in netlist (basically)
%for large circuits we should pre-allocate .in .out and .mag matrices,
%but we won't be getting that large
clear variables;
timesteps = 5000; %this is embarrasingly arbitrary.
deltat=0.01; %same as above

thresh=0.01; %more numbers pulled from thin air
%these comments are embarrasing
if(input('netlist(1) or manual input (0) = '))
    %read in netlist here
    [name, n1, n2, a1, a2, a3, a4] = textread('netlist.txt',' %s %s %s %s %s %s %s');
    %note: returns cell arrays. can input a matrix for magnitudes of
    %sources as [#,#,#,...], or a command such as
    %sin(linspace(0,2*pi,100)) as long as there are no spaces.
    nnodes = 0; %initialize this
    matsize=0; %will be counter for sources needing a current eqn
    nvoltage = 0; nvcvs = 0; nmem = 0; nphase = 0; ncoupler = 0; ncap = 0; nind = 0; nresist=0;


    %this bit requires components magnitude or relation (for components
    %with inputs) to be the 3rd value. might rewrite later to allow
    %magnitude to be last number without being kludgy
    for j=1:length(name) %using j, leaving i for sqrt(-1)
        n1{j}=str2double(n1{j});
        n2{j}=str2double(n2{j});
        if exist('a2{j}');  a2{j}=str2double(a2{j}); end
        if exist('a3{j}');  a3{j}=str2double(a3{j}); end
        if exist('a4{j}');  a4{j}=str2double(a4{j}); end
    end
    
    
    %now (also kludgy) convert the netlist into the old way of doing things
    %probably worth rewriting this to clean up in the future
    
    for j=1:length(name)
        n = eval(name{j}(2)); %use second character for component number
       switch(name{j}(1)) % check the first character of component name to identify type
           case 'V'
               vsource.in(n) = n1{j};
               vsource.out(n) = n2{j};
               vsource.o(n)=1; %i can fix this up later
               vsource.mag(n,:) = eval(a1{j}); %eval() so we get the numbers, not the string
               nnodes = max([nnodes, vsource.in(n), vsource.out(n)]);
               matsize=matsize+1;
               nvoltage=nvoltage+1;
           case 'E'
               vcvs.readin(n) = a2{j};
               vcvs.readout(n)= a3{j};
               vcvs.in(n)=n1{j};
               vcvs.in(n)=n2{j};
               vcvs.alpha{n} = a1{j}; % we leave this on as a string, we eval() later
               nnodes = max([nnodes, vcvs.in(n), vcvs.out(n),vcvs.readout(n),vcvs.readin(n)]);
               matsize=matsize+1;
               nvcvs=nvcvs+1;
           case 'D'
               mem.readin(n) = a3{j};
               mem.readout(n) = a4{j};
               mem.in(n) = a1{j};
               mem.out(n) = a2{j};
               mem.alpha(n) = a1{j};
               mem.delay(n) = a2{j};
               nnodes = max([nnodes, mem.in(n), mem.out(n),mem.readout(n),mem.readin(n)]);
               matsize=matsize+1;
               nmem=nmem+1;
           case 'P'
               phase.in(n) = n1{j};
               phase.out(n) = n2{j};
               phase.alpha(n) = eval(a1{j}); %to make it a number
               phase.phi(n) = a2{j};
               phase.delay(n) = a3{j};
               nnodes = max([nnodes, phase.in(n), phase.out(n)]);
               matsize=matsize+1;
               nphase=nphase+1;
           case 'X'
               coupler.ina(n)=n1{j};
               coupler.outa(n)=n2{j};
               coupler.inb(n)=a2{j};
               coupler.outb(n)=a3{j};
               coupler.coeff(n)=eval(a1{j}); %need a number here too
               nnodes = max([nnodes, coupler.ina(n), coupler.outa(n),coupler.inb(n),coupler.outb(n)]);
               ncoupler=ncoupler+1;
           case 'R'
               resistor.in(n)=n1{j};
               resistor.out(n)=n2{j};
               resistor.mag(n)=eval(a1{j});
               nnodes = max([nnodes, resistor.in(n), resistor.out(n)]);
               nresist=nresist+1;
           case 'C'
               cap.in(n)=n1{j};
               cap.out(n)=n2{j};
               cap.mag(n)=eval(a1{j});
               nnodes = max([nnodes, cap.in(n), cap.out(n)]);
               ncap=ncap+1;
           case 'L'
               ind.in(n)=n1{j};
               ind.out(n)=n2{j};
               ind.mag(n)=eval(a1{j});
               nnodes = max([nnodes, ind.in(n), ind.out(n)]);
               matsize=matsize+1;
               nind=nind+1;
       end
    end


    
    
    
    
else %in case you want to input the old way
    nvoltage = input('Input number of independent voltage sources = ');
    if nvoltage==0; vsource.in=[]; vsource.out=[]; end
    for n=1:nvoltage
        fprintf('Voltage Source %g \n',n);
        vsource.in(n) = input('Input Node = ');
        vsource.out(n) = input('Output Node = ');
        vsource.o(n) = input('Is Voltage Constant/Repeating (1), or not (0) = ');
        vsource.mag(n,:) = input('Voltage (input a [matrix] to vary over time) = ');
    end
    
    
    nvcvs = input('Input number of VCVSs = ');
    
    if nvcvs==0; vcvs.in=[]; vcvs.out=[]; end
    for n=1:nvcvs
        fprintf('Voltage Source %g \n',n);
        vcvs.readin(n) = input('Control  -ve = ');
        vcvs.readout(n) = input('Control Voltage +ve = ');
        vcvs.in(n) = input('Input Node = ');
        vcvs.out(n) = input('Output Node = ');
        vcvs.alpha{n} = input('Relation (use ''vdrop'' for voltage drop over control nodes) = ','s');
        %this one uses a cell because its holding a string
    end
    
    nmem = input('Input number of delayed VCVSs = ');
    if nmem==0; mem.in=[]; mem.out=[]; end
    for n=1:nmem
        fprintf('Voltage Source %g \n',n);
        mem.readin(n) = input('Control  -ve = ');
        mem.readout(n) = input('Control Voltage +ve = ');
        mem.in(n) = input('Input Node = ');
        mem.out(n) = input('Output Node = ');
        mem.alpha(n) = input('Alpha Value = ');
        mem.delay(n) = input('Delay Time (seconds) = ')/deltat; %delay in steps
    end
    
    nphase = input('Input number of phase shifters = ');
    if nphase==0; phase.in=[]; phase.out=[]; end
    for n=1:nphase
        fprintf('Phase Shift %g \n',n);
        phase.in(n) = input('Input Node = ');
        phase.out(n) = input('Output Node = ');
        phase.alpha(n) = input('Alpha (<1 is loss) = '); %%this should be in exp(alpha) form for consistency
        phase.phi(n) = input('Phase Shift (radians) = ');
        phase.delay(n) = input('Delay Time (seconds) = ')/deltat; %delay in steps
    end
    
    ncoupler = input('Input number of fiber couplers = ');
    if ncoupler==0; coupler.ina=[]; coupler.outa=[]; coupler.inb=[]; coupler.outb=[]; end
    for n=1:ncoupler
        coupler.ina(n) = input('Input node of fiber 1 = ');
        coupler.outa(n) = input('Output node of fiber 1 = ');
        coupler.inb(n) = input('Input node of fiber 2 = ');
        coupler.outb(n) = input('Output node of fiber 2 = ');
        coupler.coeff(n) = input('Coupling Coefficient = ');
    end
    
    
    
    nresist = input('Input number of resistors = ');
    if nresist==0; resistor.in=[]; resistor.out=[]; end
    for n=1:nresist
        fprintf('Resistor %g \n',n);
        resistor.in(n) = input('Input Node = ');
        resistor.out(n) = input('Output Node = ');
        resistor.mag(n) = input('Resistance = ');
    end
    
    ncap = input('Input number of capacitors = ');
    if ncap==0; cap.in=[]; cap.out=[]; end
    for n=1:ncap
        fprintf('Capacitor %g \n',n);
        cap.in(n) = input('Input Node = ');
        cap.out(n) = input('Output Node = ');
        cap.mag(n) = input('Capacitance = ');
    end
    
    nind = input('Input number of inductors = ');
    if nind==0; ind.in=[]; ind.out=[]; end
    for n=1:nind
        fprintf('Capacitor %g \n',n);
        ind.in(n) = input('Input Node = ');
        ind.out(n) = input('Output Node = ');
        ind.mag(n) = input('Inductance = ');
    end
    
nnodes = max([ max([vsource.out]) max([resistor.out]) max([cap.out]) ...
    max([vsource.in]) max([resistor.in]) max([cap.in]) max([ind.in]) ...
    max([ind.out]) max([vcvs.in]) max([vcvs.out]) max([mem.in])  ...
    max([mem.out]) max([phase.in]) max([phase.out])    ]);
end
%calculate number of nodesand branches (components)


matsize = matsize + nnodes;
% + size(vsource.mag,1) +  size(ind.in,2) + size(vcvs.in,2) + size(mem.in,2) + size(phase.in,2);
%plus buffer columns for sources


%pre allocate all our matrices
condmat = spalloc(matsize,matsize,2*size(resistor.in,2));
capmat = spalloc(matsize,matsize,2*size(cap.in,2));
source = zeros(matsize,timesteps);
nlvsource = zeros(matsize,1);
nlvsource_old = zeros(matsize,1);
voltage = zeros(matsize,1);
%now we put all the components in their respective matrices
for n=1:nvoltage
    %for standard vector of voltage sources
    if vsource.in(n)~=0;
        condmat(nnodes+n,vsource.in(n)) = -1;
    end
    if vsource.out(n)~=0;
        condmat(nnodes+n,vsource.out(n)) = 1;
        condmat(vsource.out(n),nnodes+n) = 1;
    end
    if ~vsource.o(n)
        source(nnodes+n,1:size(vsource.mag(n),2))=vsource.mag(n,:);
    else
        for t=1:timesteps
            tmp = mod(t,size(vsource.mag,2));
            %             source(nnodes+n,t)=vsource.mag(n, mod( t ,size(vsource.mag(n),2))+1);
            source(nnodes+n,t)= vsource.mag(n, tmp+1);
        end
    end
    
end

for n=1:nvcvs
    if vcvs.in(n)~=0;
        condmat(vcvs.in(n),nnodes+nvoltage+n)=condmat(vcvs.in(n),nnodes+nvoltage+n)-1;
    end
    if vcvs.out(n)~=0;
        condmat(vcvs.out(n),nnodes+nvoltage+n)=condmat(vcvs.out(n),nnodes+nvoltage+n)+1;
        condmat(nnodes+nvoltage+n,vcvs.out(n)) = 1;
    end
end


for n=1:nmem
    %for standard vector of voltage sources
    if mem.in(n)~=0;
        condmat(nnodes+nvoltage+nvcvs+n,mem.in(n)) = -1;
    end
    if mem.out(n)~=0;
        condmat(nnodes+nvoltage+nvcvs+n,mem.out(n)) = 1;
        condmat(mem.out(n),nnodes+nvoltage+nvcvs+n) = 1;
    end
    % NOPE source(nnodes+nvoltage+nvcvs+n)=mem.mag(n);
end

for n=1:nphase
    %for standard vector of voltage sources
    if phase.out(n)~=0;
        condmat(nnodes+nvoltage+nvcvs+nmem+n,phase.out(n)) = 1;
        condmat(phase.out(n),nnodes+nvoltage+nvcvs+nmem+n) = 1;
    end
end
%this might be a coupler
for n=1:ncoupler
    if coupler.ina(n)~=0 && coupler.outa(n)~=0
        condmat(coupler.outa(n),coupler.ina(n))=-(1-coupler.coeff(n));
    end
    if coupler.ina(n)~=0 && coupler.outb(n)~=0
        condmat(coupler.outb(n),coupler.ina(n))=-(coupler.coeff(n));
    end
    if coupler.inb(n)~=0 && coupler.outa(n)~=0
        condmat(coupler.outa(n),coupler.inb(n))=-(coupler.coeff(n));
    end
    if coupler.inb(n)~=0 && coupler.outb(n)~=0
        condmat(coupler.outb(n),coupler.inb(n))=-(1-coupler.coeff(n));
    end
    
end


for n=1:nresist
    if resistor.in(n)~=0;
        condmat(resistor.in(n),resistor.in(n))= condmat(resistor.in(n),resistor.in(n)) + 1/resistor.mag(n);
    end
    if resistor.out(n)~=0;
        condmat(resistor.out(n),resistor.out(n))=condmat(resistor.out(n),resistor.out(n)) + 1/resistor.mag(n);
    end
    if resistor.in(n)~=0 && resistor.out(n)~=0
        condmat(resistor.in(n),resistor.out(n)) = condmat(resistor.in(n),resistor.out(n)) - 1/resistor.mag(n);
        condmat(resistor.out(n),resistor.in(n))=condmat(resistor.out(n),resistor.in(n)) - 1/resistor.mag(n);
    end
end
for n=1:ncap
    if cap.in(n)~=0;
        capmat(cap.in(n),cap.in(n))= capmat(cap.in(n),cap.in(n)) + cap.mag(n);
    end
    if cap.out(n)~=0;
        capmat(cap.out(n),cap.out(n))= capmat(cap.out(n),cap.out(n)) + cap.mag(n);
    end
    if cap.in(n)~=0 && cap.out(n)~=0;
        capmat(cap.in(n),cap.out(n)) = capmat(cap.in(n),cap.out(n)) - cap.mag(n);
        capmat(cap.out(n),cap.in(n)) = capmat(cap.out(n),cap.in(n)) - cap.mag(n);
        
    end
end

for n=1:nind
    if ind.in(n)~=0;
        %         capmat(cap.in(n),cap.in(n))= capmat(cap.in(n),cap.in(n)) + cap.mag(n);
        condmat(nnodes+nvoltage+nvcvs+n,ind.in(n)) = -1;
        condmat(ind.in(n),nnodes+nvoltage+nvcvs+n) = 1;
    end
    if ind.out(n)~=0;
        %         capmat(cap.out(n),cap.out(n))= capmat(cap.out(n),cap.out(n)) + cap.mag(n);
        condmat(ind.out(n),nnodes+nvoltage+nvcvs+n) = -1;
        condmat(nnodes+nvoltage+nvcvs+n,ind.out(n)) = 1;
    end
    if ind.in(n)~=0 && ind.out(n)~=0;
        %         indmat(ind.in(n),ind.out(n)) = indmat(cap.in(n),cap.out(n)) - cap.mag(n);
        %         indmat(ind.out(n),ind.in(n)) = indmat(cap.out(n),cap.in(n)) - cap.mag(n);
    end
    capmat(nnodes+nvoltage+n,nnodes+nvoltage+n)=ind.mag(n);
end




%do the LU decomposition
[LL,LU,LP] = lu((capmat + deltat*condmat));
% x= voltage(:,n+1
%A=(capmat + deltat*condmat)
%b=(capmat*voltage(:,n)+source*deltat)
%y=LL\LP*b
%x=LU\(LL\LP*b)

%time marching -- should have itcheck for convergence in solid-state system)
for n=1:timesteps
    d_nlvs = inf*ones(matsize,1);
    while max(abs(d_nlvs))>thresh
        %voltage(:,n+1)= ((capmat + deltat*condmat)\(capmat*voltage(:,n)+source*deltat));
        for j=1:nmem
            if n > mem.delay(j)
                vdrop=0;
                if mem.readin(j)~=0; vdrop = -voltage(mem.readin(j),n-mem.delay(j)); end
                if mem.readout(j) ~=0; vdrop = vdrop + voltage(mem.readout(j),n-mem.delay(j)); end
                source(nnodes+nvoltage+nvcvs+j,n)= mem.alpha(j)*vdrop;
            end
        end
        for j=1:nphase
            if n > phase.delay(j)
                source(nnodes+nvoltage+nvcvs+nmem+j,n)= phase.alpha(j)*exp(-1*i*phase.phi(j))*voltage(phase.in(j),n-phase.delay(j));
            end
        end
        
        
        voltage(:,n+1)= LU\(LL\LP*(capmat*voltage(:,n)+(source(:,n)+nlvsource)*deltat));
        for j=1:nvcvs
            vdrop =0; %calculating voltage drop to use in main calculations
            if vcvs.readin(j)~=0; vdrop = -voltage(vcvs.readin(j),n+1); end
            if vcvs.readout(j)~=0; vdrop = vdrop + voltage(vcvs.readout(j),n+1); end
            nlvsource(nnodes+nvoltage+j)=eval(vcvs.alpha{j});
        end
        d_nlvs= (nlvsource-nlvsource_old);
        %         for i=1:nvcvs
        %             nlvsource(nnodes+nvoltage+i)=nlvsource(nnodes+nvoltage+i)+d_nlvs(nnodes+nvoltage+i);
        %         end
        nlvsource_old=nlvsource;
        
    end
end
figure;
hold on

%plot nicely
for n=1:matsize
    subplot(matsize,1,n);
    if n <= nnodes
        plot(1:timesteps+1,real(voltage(n,:)),1:timesteps+1,imag(voltage(n,:)));
        ts =  strcat('Node 0',num2str(n));
        title(ts);
        ylabel('Voltage (V)');
    else
        plot(1:timesteps+1,-real(voltage(n,:)),1:timesteps+1,-imag(voltage(n,:)));
        ts =  strcat('Source 0',num2str(n));
        title(ts);
        ylabel('Current (A)');
    end
end


