classdef fdfiber %An object that simulates fibers
    properties
        fiberLength
        beta1
        N
        alpha
        plotting
        Vouta
        deltat
        movie
    end
    properties 
        step_num_z
        deltaz
        uu
        uu_old
        Gmat
        vars_created = 0
    end
    methods
        function obj = fdfiber(length, beta1, N, alpha )
            %construct the object, build our G matrix, calculate steps sizes and so on.
            
            %store our parameters in the object
            obj.fiberLength = length;
            obj.beta1 = beta1;
            obj.N = N;
            obj.alpha=alpha;
            obj.step_num_z=max([obj.fiberLength*100*(obj.N) 200]);
            obj.deltaz=obj.fiberLength/obj.step_num_z;
            obj.deltat=obj.beta1*obj.deltaz;

            %preallocate matrices
            obj.uu{1}=0; %define it as a cell array
             [obj.uu{1:(obj.fiberLength/obj.deltat)}]=deal(((1+1i)*zeros(obj.step_num_z,1)));  
             obj.Gmat=zeros(obj.step_num_z);
            %and form the G matrix
            for n=1:obj.step_num_z
                obj.Gmat(n,n)=sign(obj.beta1)/obj.deltaz + obj.alpha/2; %%for upwinding
                if n-sign(obj.beta1)>0 && n-sign(obj.beta1) <= obj.step_num_z;
                    obj.Gmat(n,n-sign(obj.beta1))=-sign(obj.beta1)/obj.deltaz;
                end
            end
        end
        
        function obj = simulateFiber(obj, Vin,~,currenttstep,lasttstep,circ_deltat)
            tstep = ceil((currenttstep*circ_deltat-lasttstep*circ_deltat)/obj.deltat); %decide on number of time steps
            %we round up, and can interpolate later
            
            %actually simulate for the time steps
            for n=1:tstep
                lasttstep=lasttstep+~lasttstep; % if laststep is 0 sets it to 1
                obj.uu_old=obj.uu{lasttstep}(:); 
                obj.uu{currenttstep}(:)=obj.uu{lasttstep}(:);
                obj.uu{currenttstep}(1)=Vin; %bring in the most current part ofthe input wave
                %propogates the fiber
                obj.uu{currenttstep}(:)= obj.uu{currenttstep}(:)+((obj.deltat/obj.beta1)* ...
                    (-obj.Gmat*obj.uu{currenttstep}(:)) + ...
                    (obj.deltat/obj.beta1)*1i*obj.N.*(abs(obj.uu{currenttstep}(:)).^2).*obj.uu{currenttstep}(:));
                dXdT= (obj.uu{currenttstep}-obj.uu_old)/obj.deltat; %calculate the derivative
            end
            
            if obj.plotting
                %plot real and imaginary components
                plot(1:size(obj.uu,1),real(obj.uu(:,currenttstep)),'r-')
                hold on;
                plot(1:size(obj.uu,1),imag(obj.uu(:,currenttstep)),'g-')
                hold off;
                drawnow;
                pause(0.000001);
            end
            
            weight=mod(circ_deltat*(currenttstep-lasttstep)/obj.deltat,1); %interpolate if we had to make an extra step
            obj.Vouta=(weight* obj.uu{currenttstep}(end)+(1-weight)*obj.uu_old(end)); %and output a value to the circuit simulator
        end
        
    end
end

