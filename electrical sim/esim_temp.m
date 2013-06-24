%ESIM Simulates Basic Electrical Circuits

%read in netlist (basically)
%for large circuits we should pre-allocate .in .out and .mag matrices,
%but we won't be getting that large
clear variables;
timesteps = 100; %this is embarrasingly arbitrary.
deltat=0.01; %same as above
nvoltage = input('Input number of voltage sources = ');

if nvoltage==0; vsource.in=[]; vsource.out=[]; end
for n=1:nvoltage
    fprintf('Voltage Source %g \n',n);
    vsource.in(n) = input('Input Node = ');
    vsource.out(n) = input('Output Node = ');
    vsource.mag(n) = input('Output Voltage = ');
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
    fprintf('Capacitors %g \n',n);
    cap.in(n) = input('Input Node = ');
    cap.out(n) = input('Output Node = ');
    cap.mag(n) = input('Capacitance = ');
end
%calculate number of nodesand branches (components)
nnodes = max([ max([vsource.out]) max([resistor.out]) max([cap.out]) ...
    max([vsource.in]) max([resistor.in]) max([cap.in]) ]);
nbranches = size(vsource.in,2) + size(resistor.in,2) + size(cap.in,2);

condmat = spalloc(nnodes+1,nbranches,2*size(resistor.in,2));
capmat = spalloc(nnodes+1,nbranches,2*size(cap.in,2));
source = zeros(nnodes+1,1);
voltage = zeros(nnodes+1,1);
%now we put all the components in their respective matrices
for n=1:nvoltage
    source(vsource.in(n)+1,1)=-vsource.mag(n);
    source(vsource.out(n)+1,1)=vsource.mag(n);
end
for n=1:nresist
    condmat(resistor.in(n)+1,n)=-1/resistor.mag(n);
    condmat(resistor.out(n)+1,n)=1/resistor.mag(n);
end
for n=1:ncap
    capmat(cap.in(n)+1,nresist+n)=-1/cap.mag(n);
    capmat(cap.out(n)+1,nresist+n)=1/cap.mag(n);
end

%time marching -- should have itcheck for convergence in solid-state system)
for n=1:timesteps
    voltage(:,n+1)= ((capmat + deltat*indmat)\(V(:,n)+source*deltat);
end
figure;
hold on
for n=1:nnodes+1
    plot(voltage(n,:),n*deltat);
end


