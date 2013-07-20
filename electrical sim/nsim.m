%N2Sim Netlist Reading Circuit Simulator
%Names 		 						Syntax		
%  Voltage source								V# In Out Magnitude		(will repeat if this is a matrix)
%  Resistor										R# In Out Magnitude
%  Capacitor									C# In Out Magnitude
%  inductors									L# In Out Magnitude
%  Voltage Controlled and Delayed Source		D# In Out ControlIn ControlOut Alpha Delaytime
%  Voltage Controlled Voltage Source			E# In Out ControlIn ControlOut Alpha
%  Optical Coupler								X# In Out In 		Out 	   CouplingCoeff
%  Phase Shift									P# In Out Alpha		Phi		   Delay
%		input variables 					  name p1 p2  p3		p4			p5    p6
% i reserved for sqrt(-1)

clear variables;
timesteps = 5000; %how long to simulate for
deltat=0.01; % time in seconds for time step
thresh=0.01; %accuracy of newton-rhapson loops for NLVCVs 

%read in netlist here 
try
	[name, p1, p2, p3, p4, p5, p6] = textread('netlist.txt',' %s %s %s %s %s %s %s');
catch
	fprintf('netlist.txt not found');
	return
end
for n=1:size(name); p1{n}=str2num(p1{n}); p2{n}=str2num(p2{n}); end
%could just read in ints, but this keeps indexing consistent

%preallocate matrices (we'll sparse them later)
nnodes=max(max([ cell2mat(p1) cell2mat(p2) ])); capmat = zeros(nnodes); condmat = zeros(nnodes);

nirow = 0; %current row counter for voltage sources, to keep track from one component to the next
vlist=[]; rlist=[]; clist=[]; llist=[]; dlist=[]; elist=[]; xlist=[]; plist=[];
%make lists to keep track of which components exist where, will be useful for components needing extra-matrix math

%add all our components
for j=1:length(name)
	switch(name{j}(1))
	case 'V'
		nirow=nirow+1;
		vlist=[vlist j;nirow];
		%appends location of coponents in netlist, and their assoscated current row for later use
		Inames(nirow,:)=name{j}; %name for graphing later
		p3{j}=str2num(p3{j}); %needs to be a number (or matrix of)
		if p1{j}~=0; condmat(nnodes+nirow,p1{j}) = -1; end
		if p2{j}~=0
			condmat(nnodes+nirow,p2{j}) = 1;
			condmat(p2{j},nnodes+nirow) = 1;
		end
		for t=1:timesteps
			tmp=mod(t,size(p3{j},2));
			source(nnodes+nirow,t)=p3{j}(tmp+1);
		end
	case 'E'
		nirow=nirow+1;
		elist=[elist j;nirow];
		Inames(nirow,:)=name{j}; %name for graphing later
		p5{j}=str2num(p5{j}); 
		if p1{j}~=0; condmat(p1{j},nnodes+nirow)=1; end 
		if p2{j}~=0
			condmat(p2{j},nnodes+nirow)=1;
			condmat(nnodes+nirow,p2{j})=1;
		end
	case 'D'
		nirow=nirow+1;
		dlist=[dlist j;nirow];
		Inames(nirow,:)=name{j}; 
		p5{j}=str2num(p5{j}); 
		p6{j}=str2num(p6{j})/deltat; %make it a number and convert it to steps for later
		if p1{j}~=0; condmat(nnodes+nirow,p1{j})=-1; end
		if p2{j}~=0
			condmat(nnodes+nirow,p2{j})=1;
			condmat(nnodes+nirow,p2{j})=1;
		end
		p3{j}=str2num(p3{j}); 
	case 'P'
		nirow=nirow+1;
		plist=[plist j];
		Inames(nirow,:)=name{j};
		if p2{j}~=0
			condmat(nnodes+nirow,p2{j})=1;
			condmat(p2{j},nnodes+nirow)=1;
		end
	case 'X'
		nirow=nirow+1;
		xlist=[xlist j;nirow];
		Inames(nirow)=name{j}; 
		p5{j}=str2num(p5{j}); 
		if p1{j}~=0 && p2{j}~=0
			condmat(p2{j},p1{j})=sqrt(1-p5{j});
		end
		if p1{j}~=0 && p4{j}~=0
			condmat(p4{j},p1{j})=sqrt(p5{j});
		end
		if p2{j}~=0 && p3{j} ~=0
			condmat(p3{j},p2{j})=sqrt(p5{j});
		end
		if p2{j}~=0 && p4{j} ~=0
			condmat(p4{j},p2{j})=sqrt(1-p5{j});
		end
	case 'R'
		rlist=[rlist j];
		p3{j}=str2num(p3{j}); 
		if p1{j}~=0
			condmat(p1{j},p1{j})=condmat(p1{j},p1{j}) + 1/p3{j};
		end
		if p2{j}~=0
			condmat(p2{j},p2{j})=condmat(p2{j},p2{j}) + 1/p3{j};
		end
		if p1{j}~=0 && p2{j}~=0
			condmat(p1{j},p2{j})=condmat(p1{j},p2{j}) + 1/p3{j};
			condmat(p2{j},p1{j})=condmat(p2{j},p1{j}) + 1/p3{j};
		end
	case 'C'
		clist=[clist j];
		p3{j}=str2num(p3{j}); 
		if p1{j}~=0
			capmat(p1{j},p1{j})=capmat(p1{j},p1{j}) + p3{j};
		end
		if p2{j}~=0
			capmat(p2{j},p2{j})=capmat(p2{j},p2{j}) + p3{j};
		end
		if p1{j}~=0 && p2{j}~=0
			capmat(p1{j},p2{j})=capmat(p1{j},p2{j}) + p3{j};
			capmat(p2{j},p1{j})=capmat(p2{j},p1{j}) + p3{j};
		end
	case 'L'
		p3{j}=str2num(p3{j}); 
		nirow=nirow+1;
		llist=[llist j;nirow];
		Inames(nirow,:)=name{j}; 
		if p1{j}~=0;
			condmat(nnodes+nirow,ind.in(n)) = -1;
			condmat(ind.in(n),nnodes+nirow) = 1;
		end
		if p2{j}~=0;
			condmat(ind.out(n),nnodes+nirow) = -1;
			condmat(nnodes+nirow,ind.out(n)) = 1;
		end
		if p1{j}~=0 && p2{j}~=0; %%I have no idea what I did here. I'm sure something belongs here
		end  	
		capmat(nnodes+nirow,nnodes+nirow)=p3{j};
	end
end
%now having added lines to condmat, we need to make capmat the same size
capmat=[capmat 						zeros(size(capmat,1),nirow);...
		zeros(nirow,size(capmat,2))	zeros(nirow)				]

%that done, define our matrices for actual math using LU decomp.
[LL,LU,LP]=lu(sparse(capmat)+deltat*sparse(condmat)); %leaves cap&cond in full form for analysis
													%but completes the LU facotrization for speed
%more pre-allocations													
matsize=size(condmat,1);
nlvsource = zeros(matsize,1);
nlvsource_old = zeros(matsize,1);
voltage = zeros(matsize,timesteps);

for t=1:timesteps
	d_nlvs=inf*ones(matsize,1); %start with infinite derivative to prevent premature loop exit
	while(max(d_nlvs))>thresh 
		for nn=1:size(dlist,2) %for each memory delay
			j=dlist(1,nn); %grab the row that the delay component is in
			if t>p6{j} %if time is > than delay
				vdrop=0; %so vdrop is accurate when p4{j}==0
				if p4{j}~=0; vdrop=-voltage(p4{j},t-p6{j}); end
				if p5{j}~=0; vdrop= vdrop + voltage(p5{j},t-p6{j}); end
				source(nnodes+dlist(2,nn),t)= p5{j}*vdrop; 
				%TODO (maybe) support nonlinear here also, depreciate either E or D component
			end
		end
		for nn=1:size(plist,2)
			j=plist(1,nn);
			if t>p5{j}
				source(nnodes+plist(1,nn),t)=voltage(p1{j},t-p5{j})*exp(-1*i*p4{j})*exp(-0.5*p3{j});
			end
		end
		voltage(:,t+1)= LU\(LL\LP*(capmat*voltage(:,t)+(source(:,t)+nlvsource)*deltat)); %MATH!
		%now reconsider nonlinear voltage courses 
		%(at t>1 the above calc uses previous steps vopltages as first guess)
		for nn=1:size(elist,2)
			vdrop=0; %same as in D-component
			if p3{j}~=0; vdrop = vdrop - voltage(p3{j},t+1); end
			if p4{j}~=0; vdrop = vdrop + voltage(p4{j},t+1); end
			nlvsource(nnodes+dlist(2,nn))=eval(p5{j}); %this one's a bit different, runs user-supplied command
		end
		d_nlvs= (nlvsource- nlvsource_old); %calculate change, to compare against accuracy threshhold
		nlvsource_old=nlvsource; %set current to old for next round's comparison
	end
end


figure; %plot the results
for n=1:(nnodes+nirow)
	subplot(nnodes+nirow,1,n)
	if n<=nnodes
		plot(1:timesteps+1,real(voltage(n,:)),1:timesteps+1,imag(voltage(n,:)));
		title(['Node ' num2str(n)]);
		ylabel('Voltage (V)');
	else
		plot(1:timesteps+1,-real(voltage(n,:)),1:timesteps+1,-imag(voltage(n,:)));
		title(['Source ' Inames(n-nnodes)]);
		ylabel('Current (A)');
	end
end
return