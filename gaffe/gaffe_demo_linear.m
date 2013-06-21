%GAFFE_DEMO_LINEAR Demonstrate propagation of a linear pulse in 1D.
%
%  This function demonstrates the linear propagation of a simple 1D
%  Gaussian pulse.  This is a simple example that makes use of
%  the default end of iteration function to demonstrate how to form a basic
%  simulation. 
%
%  The Gaussian pulse disperses over 4 diffraction lengths.  Since the
%  problem is linear we must set a maximum z step that should be allowed
%  since without this the z step will become so large that it may only have
%  one intermediate sample point.
%
%  An important point to understand is that the transverse mesh resizes
%  automatically. 
%
%See also: DefaultEndIteration, GAFFE_DEMO_GAUSSIAN, GAFFE_DEMO_STEEPEN, GAFFE_DEMO_SOLITON,
%EVOLVE 

% $Author: graceej $ $Date: 2010/05/23 13:17:15 $
% $Revision: 1.1 $


% Copyright (c) 2010, Edward J. Grace
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

function gaffe_demo_linear(varargin)

% Use the default end of iteration callback.
Options = evolve('defaults');
Options.Callback.EndIteration = @DefaultEndIteration;
Options.MaxZ  = 4;
Options.MaxDZ = Options.MaxZ / 100;

% Build the coordinate system with FFT ordering.
[t, dt, omega, domega] = fftspace(10,221);

% Initial condition, a Gaussian pulse.
u0 = exp(-t.^2);

% Setup the number of snapshots we want.
N_Snap = 5;
DefaultEndIteration('set',N_Snap,Options);

Options.Callback.OperatorLinear    = @DefaultDispersion;
Options.Callback.OperatorNonlinear = @DefaultIdentity;
% Run the simulation.
evolve(u0,t,Options);

% Now get the snapshot fields.
[IsQuit, z, ts, us] = DefaultEndIteration('get');

% Plot the absolute field at each z slice.
plot(fftshift(ts{1}{1}),fftshift(abs(us{1})),'.-',...
     fftshift(ts{2}{1}),fftshift(abs(us{2})),'.-',...
     fftshift(ts{3}{1}),fftshift(abs(us{3})),'.-',...
     fftshift(ts{4}{1}),fftshift(abs(us{4})),'.-');
axis([-10 10 0 1]);
grid on;
for n=1:length(z)
    l{n} = sprintf('z=%0.2f (N=%i)',z(n),length(us{n}));
end
legend(l);
xlabel('x','FontSize',15);
ylabel('|u|','FontSize',15);
end

% $Log: gaffe_demo_linear.m,v $
% Revision 1.1  2010/05/23 13:17:15  graceej
% * Added new demonstration and improved performance.  Corrected bug in callback and rationalised forms of callbacks.
%
%
