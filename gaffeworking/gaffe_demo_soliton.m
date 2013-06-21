%GAFFE_DEMO_SOLITON Demonstrate soliton propagation.
%
%  This function demonstrates the propagation of a second order soliton in
%  a nonlinear fibre.  The 1D solution is propagated with a variable mesh
%  and variable step length and compared to the known analytical value at
%  each point.
%
%  Graphs of the solution, maximum field error, step length and grid size
%  are produced and updated regularly.  At the end of the simulation the
%  maximum observed error between the analytical field and the numerical
%  simulation is reported.
%
%  Calling with
%
%    GAFFE_DEMO_SOLITON(MESH);
%
%  Allows a different mesh scheme to be used.  For more details see
%  EVOLVEMESHSET
%
%  Additional parameters 
%
%    GAFFE_DEMO_SOLITON(MESH,TOLGLOBAL,TOLALIASED)
%
%  Control the global adaptive setp lentgth toleranace TOLGLOBAL and the
%  mesh resizing tolerance TOLALIASED.
%
%See also: GAFFE_DEMO_GAUSSIAN, GAFFE_DEMO_STEEPEN, SECOND_ORDER_SOLITON,
%EVOLVE, EVOLVEMESHSET

% $Author: graceej $ $Date: 2010/05/24 12:38:16 $
% $Revision: 1.5 $


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

function gaffe_demo_soliton(varargin)

% Define defaults.
mesh='fibonacci&humble';
TolGlobal  = 1E-9;
TolAliased = 1E-16;

if nargin > 0
    mesh=varargin{1};
end
if nargin > 1
  TolGlobal=varargin{2};
end
if nargin > 2
  TolAliased=varargin{3};
end

% Initial number of samples.
N = nearest_hamming(fibonacci(16));

% Size of real space (multiples of pi)
W = 20*pi;

% Propagate a maximum of 1 
Options.MaxZ =(2.5)*pi/2;

% Set mesh resizing to choice supplied by argument -- defaults.
Options.Callback.Mesh = evolvemeshset(mesh);

% Register this function handle as the nonlinear operator to use.
%Options.Callback.OperatorNonlinear = @DefaultKerr;
% Switch off any dispersion by using the identity linear operator.
%Options.Callback.OperatorLinear = @DefaultDispersion;
% Execute plot_progress after each iteration.
Options.Callback.EndIteration = @plot_progress;
% Evolve the field.  This does the heavy lifting.
Options.TolAliased1=TolAliased;
Options.TolAliased2=TolAliased;
Options.TolGlobal=TolGlobal;
Options.MaxIterations=1E6;
N = Options.Callback.Mesh.SetInitialSize(N);

% Set up the space.
[tau,dtau,omega,domega]= fftspace(W,N);

% Generate the initial second order soliton pulse.
u0 = second_order_soliton(tau,0);
u0 = ifft(fft(u0));

% Run the simulation.
[z,u1,tau1,omega1] = evolve(u0,tau,Options);

% Get internal values in plot_progress function.
[IsQuit, err] = plot_progress('get');

% Finally report the maximum observed infinity norm error.
fprintf('Maximum error between analytical solution and numerical method is %d\n',max(err));

end


%PLOT_PROGRESS A callback that plots progress.  Modelled in the same vain
%as DefaultEndIteration it updates a display of the solution, local and
%global errors.
function  [IsQuit, err] = plot_progress(varargin)

% This should always be set.
IsQuit = 0;

% Internal state variables.
persistent  dxi_ xi_  n_  N_ local_error_ iter_ err_ tlast tfirst;

% If we have passed parameters, then see what they are.
if nargin > 0
    switch lower(varargin{1})
        case 'set'
            throw(MException('SanityCheck:plot_progress','No internal variables to set'));
        case 'begin'
            % Called at the beginning of an iteration.
            n_ = 0;

            tic;
            tfirst = toc;
            tlast  = tfirst;
            dxi_   = zeros(100,1);
            xi_    = zeros(100,1);
            err_   = zeros(100,1);
            iter_  = zeros(100,1);
            N_     = zeros(100,1);
            local_error_ = zeros(100,1);
        case 'end'
            dxi_ = dxi_(1:n_); xi_ = xi_(1:n_); err_=err_(1:n_); 
            iter_=iter_(1:n_); N_ = N_(1:n_);   local_error_=local_error_(1:n_);
        case 'get'
            err = err_;
        otherwise
            throw(MException('SanityCheck:plot_progress','Argument to plot_progress should be ''get'', ''set'', ''begin'', ''end'' or nothing.'));
    end
    return
end

n_ = evalin('caller','n')+1;

% Get the position, step length and mesh size. 
xi_(n_)  = evalin('caller','z');
dxi_(n_) = evalin('caller','dz_try');
N_(n_)   = evalin('caller','size(u0,1)');

% Grab the current coordinate system and field.
t = fftshift(evalin('caller','x{1}'));
u = fftshift(evalin('caller','u1(:,1)'));

% Calculate local infinity norm error and get local fine / course solution
% error.
err_(n_) = max(abs(u - second_order_soliton(t,xi_(n_))));
local_error_(n_) = evalin('caller','local_error');

tnow=toc;
DT=tnow-tlast;
iter_(n_) = n_;

% Update the display every second.
if DT >= 1 || xi_(n_) >= evalin('caller','o.MaxZ - dz - dz')
    
    % Update figure 1 with the current propagation step.
    figure(1);
    clf;
    
    % Soliton profile.
    subplot(3,1,1);
    plot(t,abs(u).^2,'.',t,abs(second_order_soliton(t,xi_(n_))).^2,'o-');
    xlabel('t','FontSize',15); ylabel('Intensity','FontSize',15);
    legend('u_{gaffe}','u_{analy}');
    axis([-3 3 0 18]);
    grid on;
    
    % Step size and mesh size.
    subplot(3,1,2);
    plotyy(2*xi_(1:n_)/pi,N_(1:n_),2*xi_(1:n_)/pi,dxi_(1:n_));
    xlabel('z','FontSize',15);
    legend('N','dz');
    grid on;
    
    % Plot of local error and absolute error w.r.t. analytical soln.
    subplot(3,1,3);
    plotyy(iter_,local_error_,iter_,err_);
    grid on;
    xlabel('n','FontSize',15); 
    legend('|epsilon_{local}|','|epsilon_{abs}|');
    
    
    drawnow;
    tlast=tnow;
end

% Force the propagation engine to quit after 30 mins
if  abs(tnow-tfirst) > 60*30
    disp('Quitting since we have been going for over 30 minutes.');
    IsQuit=1;
end

end

% $Log: gaffe_demo_soliton.m,v $
% Revision 1.5  2010/05/24 12:38:16  graceej
% * Fixed yet another trivial bug with gaffe_demo_soliton.
%
% Revision 1.4  2010/05/24 12:22:33  graceej
% * Grabbing step size from wrong variable
%
% Revision 1.3  2010/05/24 08:23:52  graceej
% * Improved documentation.
% * Updated callback to follow the same as the DefaultEndIteration example.
%
% Revision 1.2  2010/05/22 21:42:47  graceej
% * Update with more reasonable display of progress.
%
% Revision 1.1  2010/05/22 20:57:59  graceej
% * Initial checkin of second order soliton demo.
%
