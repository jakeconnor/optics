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
%  Fiber or Waveguide							F# In Out Length    Beta1      N     Loss
%  Splitter                                     S# In Out Out
%  Joiner                                       J# In In  Out
%		input variables 					  name p1 p2  p3		p4	p5      p6 	     p7		p8			p9
%  Bidirectional Fiber or Waveguide				B# In In  Out       Out	Length  Beta1A   NA    LossA	  Beta1B
%														  	p10			p11	  p12
%												   		  	NB 			LossB Sigma

clear variables;
% todo: adjust timesteps and deltat based on something other than random numbers
timesteps = 1000; %how long to simulate for
deltat=0.000002; % time in seconds for time step
thresh=0.01; %accuracy of newton-rhapson loops for NLVCVs

%read in netlist here
try
    [name, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12] = textread('netlist.txt',' %s %s %s %s %s %s %s %s %s %s %s %s %s');
catch
    fprintf('netlist.txt not found');
    return
end
for n=1:size(name); p1{n}=str2num(p1{n}); p2{n}=str2num(p2{n}); end
%could just read in ints, but this keeps indexing consistent (otherwise matlab turns these into [])

%preallocate matrices (we'll sparse them later)
nnodes=max(max([ cell2mat(p1) cell2mat(p2) ])); capmat = zeros(nnodes); condmat = zeros(nnodes);

nirow = 0; %current row counter for voltage sources, to keep track from one component to the next
vlist=[]; rlist=[]; clist=[]; llist=[]; dlist=[]; elist=[]; xlist=[];
plist=[]; flist=[]; fdt=[]; slist=[]; jlist=[];
%make lists to keep track of which components exist where, will be useful for components needing extra-matrix math

%add all our components
for j=1:length(name)
    switch(name{j}(1)) %checks first character of name (handles non-unique names fine)
        case 'V'
            nirow=nirow+1; %add a current row for each component that needs one.
            vlist=[vlist [j;nirow]];
            %appends location of components in netlist, and their assoscated current row for later use
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
            dlist=[dlist [j;nirow]];
            Inames(nirow,:)=name{j};
            p5{j}=str2num(p5{j});
            p6{j}=str2num(p6{j});
            if p1{j}~=0; condmat(nnodes+nirow,p1{j})=-1; end
            if p2{j}~=0
                condmat(nnodes+nirow,p2{j})=1;
                condmat(nnodes+nirow,p2{j})=1;
            end
            p3{j}=str2num(p3{j});
        case 'P'
            nirow=nirow+1;
            p5{j}=str2num(p5{j});
            p3{j}=str2num(p3{j});
            plist=[plist [j; nirow]];
            Inames(nirow,:)=name{j};
            if p2{j}~=0
                condmat(nnodes+nirow,p2{j})=1;
                condmat(p2{j},nnodes+nirow)=1;
            end
        case 'X'
            nirow=nirow+1;
            xlist=[xlist [j]];
            Inames(nirow,:)=name{j};
            p5{j}=str2num(p5{j});
            p4{j}=str2num(p4{j});
            p3{j}=str2num(p3{j});
            %first current row
            if p1{j}~=0 && p2{j}~=0
                condmat(nnodes+nirow,p1{j})=sqrt(1-p5{j});
            end
            if p2{j}~=0 && p3{j} ~=0
                condmat(nnodes+nirow,p3{j})=sqrt(p5{j});
            end
            condmat(nnodes+nirow,p2{j})=-1;
            condmat(p2{j},nnodes+nirow)=1;
            %for second current row
            nirow=nirow+1;
            condmat(p4{j},p4{j})=-1;
            if p1{j}~=0 && p4{j}~=0
                condmat(nnodes+nirow,p1{j})=sqrt(p5{j});
            end
            
            if p2{j}~=0 && p4{j} ~=0
                condmat(nnodes+nirow,p2{j})=sqrt(1-p5{j});
            end
            condmat(nnodes+nirow,p4{j})=-1;
            condmat(p4{j},nnodes+nirow)=1;
        case 'S'
            slist=[slist j];
            nirow=nirow+1;
            p3{j}=str2num(p3{j});
            Inames(nirow,:)=name{j};
            condmat(nnodes+nirow,p1{j})=0.5;
            condmat(nnodes+nirow,p2{j})=-1/sqrt(2);
            condmat(p2{j},nnodes+nirow)=1;
            
            nirow=nirow+1;
            condmat(nnodes+nirow,p1{j})=0.5;
            condmat(nnodes+nirow,p3{j})=-1/sqrt(2);
            condmat(p3{j},nnodes+nirow)=1;
            
        case 'J'
            nirow=nirow+1;
            Inames(nirow,:)=name{j};
            jlist=[jlist [j;nirow]];
            p3{j}=str2num(p3{j});
            condmat(nnodes+nirow,p1{j})=-sqrt(2);
            condmat(nnodes+nirow,p2{j})=-sqrt(2);
            condmat(nnodes+nirow,p3{j})=2;
            condmat(p3{j},nnodes+nirow)=1;
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
                condmat(p1{j},p2{j})=condmat(p1{j},p2{j}) - 1/p3{j};
                condmat(p2{j},p1{j})=condmat(p2{j},p1{j}) - 1/p3{j};
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
                capmat(p1{j},p2{j})=capmat(p1{j},p2{j}) - p3{j};
                capmat(p2{j},p1{j})=capmat(p2{j},p1{j}) - p3{j};
            end
        case 'L'
            p3{j}=str2num(p3{j});
            nirow=nirow+1;
            llist=[llist j;nirow];
            Inames(nirow,:)=name{j};
            if p1{j}~=0;
                condmat(nnodes+nirow,p1{j}) = -1;
                condmat(p1{j},nnodes+nirow) = 1;
            end
            if p2{j}~=0;
                condmat(p2{j},nnodes+nirow) = -1;
                condmat(nnodes+nirow,p2{j}) = 1;
            end
            
            capmat(nnodes+nirow,nnodes+nirow)=p3{j};
        case 'F'
            nirow=nirow+1;
            flist=[flist [j; nirow]];
            Inames(nirow,:)=name{j};
            p3{j}=str2num(p3{j});
            p4{j}=str2num(p4{j});
            p5{j}=str2num(p5{j});
            p6{j}=str2num(p6{j});
            fiber{j} = fdfiber(p3{j},p4{j},p5{j},p6{j});
            if p2{j}~=0
                condmat(nnodes+nirow,p2{j})=1;
                condmat(p2{j},nnodes+nirow)=1;
            end
            fiber{j}.plotting=0;
            fdt = [fdt fiber{j}.deltat];
        case 'B'
            p3{j}=str2num(p3{j});
            p4{j}=str2num(p4{j});
            p5{j}=str2num(p5{j});
            p6{j}=str2num(p6{j});
            p7{j}=str2num(p7{j});
            p8{j}=str2num(p8{j});
            p9{j}=str2num(p9{j});
            p10{j}=str2num(p10{j});
            p11{j}=str2num(p11{j});
            p12{j}=str2num(p12{j});
            
            
            
            nirow=nirow+1;
            flist=[flist [j; nirow]];
            Inames(nirow,:)=name{j};
            
            fiber{j} = bfib4(p5{j},p6{j},p7{j},p8{j},p9{j},p10{j},p11{j},p12{j});
            if p4{j}~=0
                condmat(nnodes+nirow,p4{j})=1;
                condmat(p4{j},nnodes+nirow)=1;
            end
            
            
            nirow=nirow+1;
            if p3{j}~=0
                condmat(nnodes+nirow,p3{j})=1;
                condmat(p3{j},nnodes+nirow)=1;
            end
            
            fiber{j}.sigma=p7{j};
            fiber{j}.plotting= 1;
            fdt = [fdt fiber{j}.deltat];
    end
end

condmat=[condmat zeros(size(condmat,1),nnodes+nirow-size(condmat,2));
    zeros(nnodes+nirow-size(condmat,1),size(condmat,1))...
    zeros(nnodes+nirow-size(condmat,2),nnodes+nirow-size(condmat,1)) ];
%now having added lines to condmat, we need to make capmat the same size. source sometimes too
capmat=[capmat 						zeros(size(capmat,1),size(condmat,2)-size(capmat,2));...
    zeros(size(condmat,1)-size(capmat,1),size(capmat,2))	zeros(size(condmat,1)-size(capmat,1))	];

%that done, define our matrices for actual math using LU decomp.
[LL,LU,LP]=lu(sparse(capmat)+deltat*sparse(condmat)); %leaves cap&cond in full form for analysis
%but completes the LU facotrization for speed
%pre-allocations
matsize=size(condmat,1);
source=[source; zeros(nnodes+nirow-size(source,1),size(source,2))];


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
                %TODO: support nonlinear here also, depreciate either E or D component
            end
        end
        for nn=1:size(plist,2)
            j=plist(1,nn);
            if t>p5{j}
                source(nnodes+plist(2,nn),t)=voltage(p1{j},t-p5{j})*exp(-1i*eval(p4{j}))*exp(-0.5*p3{j});
            end
        end
        
        for nn=1:size(flist,2) %for optical fiber
            j=flist(1,nn); %call the object relevant to the fiber, drop output in voltage source vector
            fiber{j}=fiber{j}.simulateFiber(voltage(p1{j},t),voltage(p2{j},t),t,(t-1),deltat);
            source(nnodes+flist(2,nn),t)=fiber{j}.Vouta;
            try    
            source(nnodes+flist(2,nn)+1,t)=fiber{j}.Voutb;
            catch
            end
            
        end
        
        voltage(:,t+1)= LU\(LL\LP*(capmat*voltage(:,t)+(source(:,t)+nlvsource)*deltat)); %MATH!
        % recalculate nonlinear voltage courses
        %(at t>1 the above calc uses previous steps vopltages as first guess)
        for nn=1:size(elist,2)
            vdrop=0; %same as in D-component
            if p3{j}~=0; vdrop = vdrop - voltage(p3{j},t+1); end
            if p4{j}~=0; vdrop = vdrop + voltage(p4{j},t+1); end
            try
            nlvsource(nnodes+dlist(2,nn))=eval(p5{j}); %runs user-supplied command
            catch
                error('Nonlinear voltage source relation is invalid')
            end
        end
        d_nlvs= (nlvsource - nlvsource_old); %calculate change, to compare against accuracy threshhold
        nlvsource_old=nlvsource; %set current to old for next round's comparison
    end
end

figure; %plot the results
for n=1:(nnodes+nirow)
    if n<=nnodes
        subplot(2,max(nnodes,nirow),n)
        plot(1:timesteps+1,real(voltage(n,:)),1:timesteps+1,imag(voltage(n,:)));
        title(['Node ' num2str(n)]);
        ylabel('Voltage (V)');
    else
        subplot(2,max(nnodes,nirow),n)
        plot(1:timesteps+1,-real(voltage(n,:)),1:timesteps+1,-imag(voltage(n,:)));
        title(['Source ' Inames(n-nnodes)]);
        ylabel('Current (A)');
    end
end
return
