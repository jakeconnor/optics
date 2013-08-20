classdef bfib3
    %bfib3
    
    
    properties
        fiberLength
        beta1
        N
        loss
        plotting
        Vouta
        Voutb
        deltat
        step_num_z
        deltaz
        uu
        uu_old
        Gmat
        sigma
        interp
        currentstep
        laststep
        time
        oldstep
    end
    
    methods
        function obj=bfib3(length,beta1a,na,lossa,beta1b,nb,lossb,sigma)
            obj.fiberLength=length;
            obj.beta1{1}=beta1a;
            obj.beta1{2}=beta1b;
            obj.N{1}=na;
            obj.N{2}=nb;
            obj.loss{1}=lossa;
            obj.loss{2}=lossb;
            obj.sigma=sigma;
            
            %build our two fibers
            for nfib=1:2
                obj.step_num_z{nfib}=max([obj.fiberLength*100*obj.N{nfib} 200]);
                obj.deltaz{nfib}=obj.fiberLength/obj.step_num_z{nfib};
                obj.deltat{nfib}=obj.beta1{nfib}*obj.deltaz{nfib};
                obj.uu{nfib}=((1+1i)*zeros(obj.step_num_z{nfib},obj.fiberLength/obj.deltat{nfib}));%%
                obj.Gmat{nfib}=zeros(obj.step_num_z{nfib});
                for n=1:obj.step_num_z{nfib}
                    obj.Gmat{nfib}(n,n)=sign(obj.beta1{nfib})/obj.deltaz{nfib} + obj.loss{nfib}/2; %%for upwinding
                    if n-sign(obj.beta1{nfib})>0 & n-sign(obj.beta1{nfib}) <= obj.step_num_z{nfib};
                        obj.Gmat{nfib}(n,n-sign(obj.beta1{nfib}))=-sign(obj.beta1{nfib})/obj.deltaz{nfib};
                    end
                end
                spacing{nfib}=obj.deltaz{nfib}*(1:size(obj.uu{nfib},1));
            end
            %next bit requires both spacing vectors
            for nfib=1:2
                oppfib=mod(nfib,2)+1;%quick var for what the number of the other fiber is
                obj.interp{nfib}= zeros(size(obj.uu{nfib},1),size(obj.uu{oppfib},1)); %interp{nfib} is mult by oppfib to get size of nfib
                for ni=1:size(obj.uu{nfib},1)
                    for nj=2:size(oppfib,1)
                        if spacing{nfib}(ni)<spacing{oppfib}(nj) %build our interpolation matrix
                            zweight=1 - (spacing{oppfib}(nj)-spacing{nfib}(ni))/obj.deltaz{oppfib};
                            obj.interp{nfib}(ni,nj)=zweight;
                            obj.interp{nfib}(ni,nj-1)=1-zweight;
                        end
                    end
                end
                obj.currentstep{nfib}=1;
                obj.laststep{nfib}=0;
                obj.time{nfib}=0;
            end
        end
        
        
        
        %% simulateFiber: propogates the fiber for some time
        function obj = simulateFiber(obj, Vina, Vinb, circstep, ~, circdeltat)
            Vin{1}=Vina; Vin{2}=Vinb;
            
            %this time around the waves keep track of real time, so if the
            
            while obj.time{1} < circstep*circdeltat || obj.time{2} < circstep*circdeltat
                nfib=1+(obj.time{1}>obj.time{2}); %sets actively simming fiber to whichever has lower time
                oppfib=mod(nfib,2)+1; %set opposite fiber
                obj.laststep{nfib}=obj.laststep{nfib}+~obj.laststep{nfib}; % if lasstep is 0 sets it to 1
                obj.uu_old{nfib}=obj.uu{nfib}(:,obj.laststep{nfib});
                obj.uu{nfib}(:,obj.currentstep{nfib})=obj.uu{nfib}(:,obj.laststep{nfib});
                obj.uu{nfib}(1,obj.currentstep{nfib})=Vin{nfib}; %bring in the most current part ofthe input wave
                
                
                obj.uu{nfib}(:,obj.currentstep{nfib})= obj.uu{nfib}(:,obj.currentstep{nfib})+((obj.deltat{nfib}/obj.beta1{nfib})* ...
                    (-obj.Gmat{nfib}*obj.uu{nfib}(:,obj.currentstep{nfib})) + ...
                    (obj.deltat{nfib}/obj.beta1{nfib})*1i ...
                    *(obj.N{nfib}.*(abs(obj.uu{nfib}(:,obj.currentstep{nfib})).^2)...
                    + obj.sigma.*(abs(obj.uu{oppfib}(:,obj.currentstep{oppfib})).^2).* ...
                    obj.uu{nfib}(:,obj.currentstep{nfib})));

                dXdT{nfib}= (obj.uu{nfib}(obj.currentstep{nfib})-obj.uu_old{nfib})/obj.deltat{nfib};
                
                obj.laststep{nfib}=obj.currentstep{nfib};
                obj.time{nfib}=obj.time{nfib}+obj.deltat{nfib};
                obj.currentstep{nfib}=obj.currentstep{nfib}+1;
            end
            
            for nfib=1:2
                weight{nfib}=mod(circdeltat*(obj.currentstep{nfib}-obj.laststep{nfib})/obj.deltat{nfib},1);
            end
            
            obj.Vouta=(weight{1}*obj.uu{1}(end,obj.currentstep{1}-1)+(1-weight{1})*obj.uu_old{1}(end));
            obj.Voutb=(weight{2}*obj.uu{2}(end,obj.currentstep{2}-1)+(1-weight{2})*obj.uu_old{2}(end));
            
            
        end
        
        
        
    end
    
end

