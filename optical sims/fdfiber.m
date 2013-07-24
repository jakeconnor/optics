classdef fdfiber %An object that simulates fibers
    properties
        fiberLength
        beta1
        N
        alpha
        plotting
        Vout
    end
    properties
        step_num_z
        deltaz
        deltat
        uu
        uu_old
        Gmat
        vars_created  = 0
    end
    methods
        function obj = fdfiber(length, beta1, N, alpha )
            if nargin>0
                obj.fiberLength = length;
                obj.beta1 = beta1;
                obj.N = N;
                obj.alpha = alpha;
            end
        end
        
        function [obj] = simulateFiber(obj, Vin, time)
            if ~obj.vars_created
                obj.step_num_z=obj.fiberLength*1000*(1+0.1*obj.N);
                obj.deltaz=obj.fiberLength/obj.step_num_z;
                obj.deltat=obj.beta1*obj.deltaz;
                obj.uu=((1+1i)*zeros(1,obj.step_num_z))';
                obj.Gmat=zeros(obj.step_num_z);
                for n=1:obj.step_num_z
                    obj.Gmat(n,n)=sign(obj.beta1)/obj.deltaz + obj.alpha/2; %%for upwinding
                    if n-sign(obj.beta1)>0 & n-sign(obj.beta1) <= obj.step_num_z;
                        obj.Gmat(n,n-sign(obj.beta1))=-sign(obj.beta1)/obj.deltaz;
                    end
                end
                obj.vars_created=1;
            end
            tstep = ceil(time/obj.deltat); %decide on number of time steps
            %we round up, and can interpolate later
            
            %actually simulate for the time steps
            
            for n=1:tstep
                obj.uu_old=obj.uu;
                obj.uu(1)=Vin; %bring in the most current part ofthe input wave
                obj.uu= obj.uu + (obj.deltat/obj.beta1*(-obj.Gmat*obj.uu) + (obj.deltat/obj.beta1)*1i*obj.N.*(abs(obj.uu).^2).*obj.uu);
                dXdT= (obj.uu-obj.uu_old)/obj.deltat;
                dXdT2=dXdT/obj.deltat;
                if obj.plotting
                    plot(1:length(obj.uu),real(obj.uu),'r-')%, axis([0 step_num_z 0 1]);
                    hold on;
                    plot(1:length(obj.uu),imag(obj.uu),'g-')
                    hold off;
                    drawnow;
                    pause(0.000001);
                end
            end
            
            
            weight=mod(time/obj.deltat,1); %interpolate if we had to make an extra step
            obj.Vout=(weight*obj.uu(end)-(1-weight)*obj.uu_old(end))/2;
            
            
            
        end
        
    end
end

