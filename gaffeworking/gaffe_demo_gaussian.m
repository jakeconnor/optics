%GAFFE_DEMO_GAUSSIAN Demonstrate the collapse threshold for a Gaussian.
%   Given an initial Gaussian field and a power relative to the critical
%   power given in [1] this propagates the Gaussian two nominal diffraction
%   lengths or until five minutes have passed.  
%
%   It then plots the beam spot-size v.s. the propagation distance.
%
%   Case study:
%
%      From equation 4 where the complex field u is u(x,y,z):
%
%            d      d^2     d^2
%         i --- u + ---  u + ---  u + |u|^2 u = 0
%           dz      dx^2    dy^2
%
%     a) Identify the integration direction z:
%     b) The linear operator is then terms involving d^2
%     c) The nonlinear operator is the term involving |u|^2
%     e) Transform linear operator to Fourier domain, replace d/dx with 
%        i kx and d/dy with iky so:
%
%     L = -i.*(KX{1}.^2 + KX{2}.^2);
%     N = i.*abs(u).^2;
%
%  Where KX{1} is the first spatial frequency (kx) and KX{2} is the second 
%  transverse spatial frequency (ky).
%      
%  Ncr for Gaussian = 1.8962 (table 1)
%
% When modelling cases involving or near blowup or collapse small numerical
% errors are important. For this reason be prepared to wait as you may need
% to adjust TolGlobal, TolAliased1 and TolAliased2.
%
% GAFFE_DEMO_GAUSSIAN(POWER) Scales the power of the Gaussian beam to Power
% (0.5) by default.  Here 1.0 is the critical power.  For a value close to
% 1 the run-time can be very large.
%
% GAFFE_DEMO_GAUSSIAN(POWER,MESH) Selects a different mesh resizing
% technique, for valid options see EVOLVEMESHSET.
% 
%  
% [1] Fibich "Critical power for self-focusing in bulk media
%     and in hollow waveguides", Optics Letters 25:5 p335
%
%See also: GAFFE_DEMO_STEEPEN, GAFFE_DEMO_SOLITON, EVOLVE

% $Author: graceej $ $Date: 2010/05/23 14:37:24 $
% $Revision: 1.6 $


% Copyright (c) 2009, Edward J. Grace
% All rights reserved.
 
% Redistribution and use in source and binary forms, with or 
% without modification, are permitted provided that the following 
% conditions are met:
 
%     * Redistributions of source code must retain the above 
%       copyright notice, this list of conditions and the following 
%       disclaimer.
%     * Redistributions in binary form must reproduce the above 
%       copyright notice, this list of conditions and the following 
%       disclaimer in the documentation and/or other materials 
%       provided with the distribution.
%     * Neither the name of the Imperial College London nor the 
%       names of its contributors may be used to endorse or 
%       promote products derived  this software without specific 
%       prior written permission.
 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND 
% CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, 
% INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF 
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS 
% BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED 
% TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
% DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON 
% ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR 
% TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF 
% THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF 
% SUCH DAMAGE.

function gaffe_demo_gaussian(varargin)
Ncr=1.86225; % Townes boundary (Blowup)
alpha_Gaussian=1.8962;  % collapse boundary for Gaussian.
factor=0.5; 
mesh='defaults';
for n=0:nargin-1
    if isnumeric(varargin{n+1})
        factor=varargin{n+1};
    else
        mesh=varargin{n+1};
    end
end


%% Sample the profile over a large dense grid. The propagation engine is
% smart enough to downsample or upsample appropriately.
[x,dx,kx,dkx]= fftspace(12*pi,240);
[X,Y]=ndgrid(x,x); [KX,KY]=ndgrid(kx,kx);

Pc=2*pi*alpha_Gaussian;

% Initial field
u0_axis = sqrt(factor*Pc/pi);
u = u0_axis*exp(-0.5*(X.^2+Y.^2));
% Report power and critial power
Pm=sum(sum(abs(u).^2))*dx*dx;
expectation = 'Global existance (this should be quick)';
if (Pm >= 2*pi*1.86225) % Townes limit
    expectation='Blowup (be patient)';
    if Pm >= 2*pi*1.8962 % Gaussian collapse
        expectation='Collapse (be very patient)'
    end
end
       
    
fprintf('Power in model (Pm)          : %f\n',Pm);
fprintf('Critical power (Pc) (2*pi*Nc): %f\n',Pc);
fprintf('P/Pc                         : %f\n',Pm/Pc);
fprintf('Expect                       : %s\n\n\n',expectation);

%% Get default options and display them.
Options = evolve('defaults');
Options.Callback.Mesh = evolvemeshset(mesh);

%% Setup propagation.
Ld=0.5;
z0=2*Ld;
Options.MaxZ = z0;


%% Nonlinear forward propagation.
Options.MaxZ = z0;
Options.Callback.OperatorLinear =    @DefaultDiffractionFibich;
Options.Callback.OperatorNonlinear = @DefaultKerr;
Options.Callback.EndIteration =      @SampleField;
Options.Trace = 0;


% Solve the problem.
[zfinal,u1,X1,KX1] = evolve(u,x,x,Options);
disp('Finished nonlinear forward propagation.');

% Plot results
figure(1);
clf;
[IsQuit, n, z, wx, wy, Nx, u0] = SampleField('get');
[ax,h1,h2]=plotyy(z/Ld,wx./sqrt(2),z/Ld,abs(u0).^2/abs(u0_axis).^2);
set(ax(1),'Fontsize',15);
set(ax(2),'Fontsize',15);
xlabel(ax(1),'z/Ld','FontSize',18);
xlabel(ax(2),'z/Ld','FontSize',18);
ylabel(ax(1),'w_x/w_x(0)','FontSize',18);
ylabel(ax(2),'I_{axis}/I_{axis}(0)','FontSize',18);
grid on;
set(h2,'Marker','.');
end

%Sample the field at each iteration.  This is executed after each iteration
%and accumulates the samples in internal variables.  These are returned
%when an argument is supplied.
function [IsQuit,n_,z_,wx_,wy_,Nx_, u0_] = SampleField(varargin)

% By default do not quit.
IsQuit = 0;

% Maintain internal state between callbacks.  This is where we log the 
% various values of interest.
persistent n z wx wy Nx T u0;

% If the argument in is 'begin' then initialise everything and return.
if nargin == 1
    switch lower(varargin{1})
        case 'begin'
            tic;
            T = toc;
            
            n = 0;
        
            z = zeros(100,1);
            wx=z; wy=z; Nx=z; u0=z;
            
            fprintf('\n\nProgress:     ');
        case 'end'
            % Finalise the state of the variables.
            z = z(1:n); wx=wx(1:n); wy=wy(1:n); Nx=Nx(1:n); u0=u0(1:n);
            
            fprintf('\n\nFinished!\n');
        case 'get'
            % Relay persistent internal variables to external function.
            n_=n; z_=z; wx_=wx; wy_=wy; Nx_=Nx; u0_=u0;
            
        case 'set'
            % A blank argument should be considered the same as no
            % argument.  This allows other values to be set.
        otherwise
            throw(MException('SanityCheck:EndIteration',sprintf('The end of iteration callback should be called with either ''begin'', ''end'' or ''get'', currently ''%s'' is not supported.',varargin{1})));            
    end
    return
end



% Carry on as normal, collecting data.  Here we collect the current
% iteration count in the caller (evolve) and various aspects of the
% recently propagated field (u1).
n = evalin('caller','n')+1;
z(n)  = evalin('caller','z');
wx(n) = evalin('caller','width(u1,X{1})/2');
wy(n) = evalin('caller','width(u1,X{2})/2');
Nx(n) = evalin('caller','N(1)');
u0(n) = evalin('caller','u1(1,1)');

% Here we update the progress.
fprintf('\b\b\b\b%0.3d%%',floor(100*(z(n)/evalin('caller','o.MaxZ'))));

% If the timeout has been exceeded (5 minutes) then we set IsQuit to 1.
if toc-T > 60*5
    fprintf('Quitting as run time has exceeded 5 minutes\n');
    IsQuit = 1;
end

end
%% CVS Log
%
% $Log: gaffe_demo_gaussian.m,v $
% Revision 1.6  2010/05/23 14:37:24  graceej
% * Cross reference incorrect.
%
% Revision 1.5  2010/05/23 13:17:15  graceej
% * Added new demonstration and improved performance.  Corrected bug in callback and rationalised forms of callbacks.
%
% Revision 1.4  2010/05/22 20:57:59  graceej
% * Initial checkin of second order soliton demo.
%
% Revision 1.3  2009/10/24 11:08:05  graceej
% * Modified license to make use of the current BSD style open source initiative license.
%
% Revision 1.2  2009/05/06 20:09:45  graceej
% * Updated make package target to include distribution of the license file and README.
%
% Revision 1.1  2009/05/06 16:30:34  graceej
% * Moved demos to main branch.
%
% Revision 1.2  2009/05/06 16:23:02  graceej
% * Updated documentation to reflect options.
%
% Revision 1.1  2009/05/06 09:11:06  graceej
% * Rename demo with gaffe_demo_ prefix.
%
% Revision 1.2  2009/04/19 19:03:49  graceej
% * Reworked gaussian ray transfer matrix code.
% * Corrected minor bugs in demos, using the old field for sampling rather than the newly propagated field.
%
% Revision 1.1  2009/04/17 11:42:32  graceej
% * Brought over auto_checking/BPC-LIB tag=end and checked in.
% * This old repository is now defunct.
%
% Revision 1.2  2009/04/17 11:37:27  graceej
% * Wholescale modification of the entire library.
%
% Revision 1.1.2.3  2009/01/20 21:06:13  graceej
% * Fixed, now modeling what we thought we were modleling.
%
% Revision 1.1.2.2  2009/01/20 21:01:33  graceej
% * Corrected to use appropriate initial condition.
%
% Revision 1.1.2.1  2009/01/20 14:31:31  graceej
% * Added diffraction operator appropriate for Fibich paper.
%
% Revision 1.1.2.2  2008/09/18 21:33:50  graceej
% * Modified to work more smoothly in octave.
%
% Revision 1.1.2.1  2008/09/16 16:21:50  graceej
% * Initial checkin.
%
%