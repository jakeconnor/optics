%ESIM Simulates Basic Electrical Circuits

%read in netlist (basically)
%for large circuits we should pre-allocate .in .out and .mag matrices,
%but we won't be getting that large
clear variables;
timesteps = 10000; %this is embarrasingly arbitrary.
deltat=0.01; %same as above
nvoltage = input('Input number of independent voltage sources = ');

if nvoltage==0; vsource.in=[]; vsource.out=[]; end
for n=1:nvoltage
    fprintf('Voltage Source %g \n',n);
    vsource.in(n) = input('Input Node = ');
    vsource.out(n) = input('Output Node = ');
    vsource.o(n) = input('Is Voltage Constant/Repeating (1), or not (0) = ');
    vsource.mag = input('Voltage (input a [matrix] to vary over time) = ');
end

nvcvs = input('Input number of VCVSs = ');

if nvcvs==0; vcvs.in=[]; vcvs.out=[]; end
for n=1:nvcvs
    fprintf('Voltage Source %g \n',n);
    vcvs.readin(n) = input('Control  -ve = ');
    vcvs.readout(n) = input('Control Voltage +ve = ');
    vcvs.in(n) = input('Input Node = ');
    vcvs.out(n) = input('Output Node = ');
    vcvs.alpha(n) = input('Alpha Value = ');
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
%calculate number of nodesand branches (components)
nnodes = max([ max([vsource.out]) max([resistor.out]) max([cap.out]) ...
    max([vsource.in]) max([resistor.in]) max([cap.in]) max([ind.in]) ...
    max([ind.out]) max([vcvs.in]) max([vcvs.out]) max([mem.in]) max([mem.out]) ]);

matsize= nnodes + size(vsource.mag,1) + ...
    size(ind.in,1) + size(vcvs.in,1) + size(mem.in,1);%plus buffer columns for sources

condmat = spalloc(matsize,matsize,2*size(resistor.in,2));
capmat = spalloc(matsize,matsize,2*size(cap.in,2));
source = zeros(matsize,timesteps);
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
        source(nnodes+n,1:size(vsource.mag(n)),2)=vsource.mag(n);
    else
        for t=1:timesteps
           source(nnodes+n,t)=vsource.mag(mod(t,size(vsource.mag,2))+1);
        end
    end
    
end

for n=1:nvcvs
    if vcvs.readin(n)~=0; condmat(nnodes+nvoltage+n,vcvs.readin(n))=1; end
    if vcvs.readout(n)~=0; condmat(nnodes+nvoltage+n,vcvs.readout(n))=-1; end
    if vcvs.in(n)~=0;
        condmat(nnodes+nvoltage+n,vcvs.in(n))=condmat(nnodes+nvoltage+n,vcvs.in(n))-vcvs.alpha(n);
        condmat(vcvs.in(n),nnodes+nvoltage+n)=condmat(vcvs.in(n),nnodes+nvoltage+n)-1;
    end
    if vcvs.out(n)~=0;
        condmat(nnodes+nvoltage+n,vcvs.out(n))=condmat(nnodes+nvoltage+n,vcvs.out(n))+vcvs.alpha(n);
        condmat(vcvs.out(n),nnodes+nvoltage+n)=condmat(vcvs.out(n),nnodes+nvoltage+n)+1;
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
        condmat(nnodes+nvoltage+n,ind.in(n)) = -1;
        condmat(ind.in(n),nnodes+nvoltage+n) = 1;
    end
    if ind.out(n)~=0;
        %         capmat(cap.out(n),cap.out(n))= capmat(cap.out(n),cap.out(n)) + cap.mag(n);
        condmat(ind.out(n),nnodes+nvoltage+n) = -1;
        condmat(nnodes+nvoltage+n,ind.out(n)) = 1;
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
    %voltage(:,n+1)= ((capmat + deltat*condmat)\(capmat*voltage(:,n)+source*deltat));
    for j=1:nmem
        if n > mem.delay(j)
            vdrop=0;
            if mem.readin(j)~=0; vdrop = voltage(mem.readin(j),n-mem.delay(j)); end
            if mem.readout(j) ~=0; vdrop = vdrop + voltage(mem.readout(j),n-mem.delay(j)); end
            source(nnodes+nvoltage+nvcvs+j,n)= mem.alpha(j)*vdrop;
        end
    end
    
    voltage(:,n+1)= LU\(LL\LP*(capmat*voltage(:,n)+source(:,n)*deltat));
    
end
figure;
hold on
for n=1:nnodes
    plot(voltage(n,:),n*deltat);
end


