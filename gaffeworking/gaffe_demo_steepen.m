%GAFFE_DEMO_STEEPEN Demonstrate the effects of self-steepening.
%  This demo uses GAFFE to demonstrate the effects of self-steepening on an
%  optical pulse. It numerically reproduces figure 4.15 from Agrawal
%  "Nonlinear Fiber Optics" 1989 ISBN 0-12-045140-9
%
%  This example pushes the limits of this technique as the pulse
%  ultimately develops a shock, hence is poorly represented by a
%  Fourier decomposition as it becomes less smooth.
%
%  By adjusting the TolGlobal parameter and TolAliased parameter
%  this nicely demonstrates the requirement for increasing the
%  (spectral) mesh size as nonlinearity broadens the spectrum.
%  Setting TolGlobal to a very small value is needed to maintain a
%  close solution to the actual (analytic) solution, this is highly
%  computationally demanding as the step-size must be very small.
%
%See also: GAFFE_DEMO_GAUSSIAN, GAFFE_DEMO_SOLITON, EVOLVE

% $Author: graceej $ $Date: 2010/05/24 08:24:19 $
% $Revision: 1.8 $


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

function gaffe_demo_steepen(varargin)
% Self steepening factor
S=0.02; 
% Global tolerance
TolGlobal=1E-9;
% Aliasing tolerance
TolAliased=1E-13;
% Number of samples
N = nearest_humble(fibonacci(10));
% Size of real space (multiples of pi)
W = 6*pi;


mesh='fibonacci&humble';
if nargin > 0
    mesh=varargin{1};
end
if nargin > 1
  TolGlobal=varargin{2};
end
if nargin > 2
  TolAliased=varargin{3};
end


% Initialise the plot_progress function - set up plots and internal data.
plot_progress('set',S);

% Propagate a maximum of 1 
Options.MaxZ =0.27/S;

% Set mesh resizing to choice supplied by argument -- defaults.
Options.Callback.Mesh = evolvemeshset(mesh);

% Define a function handle that describes the effects of self-steepening
% controlled by the parameter S and Kerr self-phase modulation. 
SelfSteepeningAndKerr = @(dxi, omega, tau, u, fftu) ...
    DefaultSelfSteepening(S*dxi,omega,tau,u,fftu)% .* ...
    %DefaultKerr(dxi,[],[],u,[]);

% Register this function handle as the nonlinear operator to use.
Options.Callback.OperatorNonlinear = SelfSteepeningAndKerr;
% Switch off any dispersion by using the identity linear operator.
Options.Callback.OperatorLinear = @DefaultIdentity;
% Execute plot_progress after each iteration.
Options.Callback.EndIteration = @plot_progress;
% Evolve the field.  This does the heavy lifting.
Options.TolAliased1=TolAliased;
Options.TolAliased2=TolAliased;
Options.TolGlobal=TolGlobal;
Options.MaxIterations=1E9;

N = Options.Callback.Mesh.SetInitialSize(N);
Options.Callback.Mesh.IsAliasDown1 = @IsAliasDown1_Custom;


[tau,dtau,omega,domega]= fftspace(W,N);

% Generate the initial Gaussian pulse.
u0 = sqrt(exp(-tau.^2));
u0=ifft(fft(u0));

[z,u1,tau1,omega1] = evolve(u0,tau,Options);
end

%PLOT_PROGRESS A callback that plots progress.  This shou§ld be registered
%to be called at the end of each iteration, at which point it will plot the
%fields and spot sizes etc.
function  IsQuit = plot_progress(varargin)

% By default do not quit.
IsQuit = 0;

% Maintain internal state between callbacks.
persistent S n N T z tlast;

if nargin > 0
    switch lower(varargin{1})
        case 'begin'
            % This is called once before we begin iterating.
            tic
            T = toc;
            tlast = T;
            
            n = 0;
            
            z  = zeros(100,1); 
            N  = zeros(100,1);
        case 'end'
            % This is called once after we have finished iterating.
            z = z(1:n); N = N(1:n);
            
        case 'get'
            % This should be used to return values from the internal state.
        case 'set'
            % This should be used to set values to the internal state.
            S = varargin{2};
        otherwise
            throw(MException('SanityCheck:EndIteration',sprintf('The end of iteration callback should be called with either ''begin'', ''end'' or ''get'', currently ''%s'' is not supported.',varargin{1})));            
    end
    return
end

n = evalin('caller','n') + 1;
N(n) = evalin('caller','N(1)');
z(n) = evalin('caller','z');

%t0_(n_) = width(evalin('caller','u0'),evalin('caller','X{1}'))/2;
%N_(n_) = evalin('caller','size(u0,1)');
%tnow=toc;
DT=toc - tlast;

% Update progress plot at most every second.
if DT >= 1 || z(n) >= evalin('caller','o.MaxZ - dz - dz')
    
    % Pull in the current field and fourier transform.
    t     = fftshift(evalin('caller','x{1}'));
    omega = fftshift(evalin('caller','kx{1}'));
    u     = fftshift(evalin('caller','u1(:,1)'));
    U     = fftshift(evalin('caller','U1(:,1)'));

    % Implicit function for intensity as a function of xi at a given time t.
    % this function should be zero when the intensity is correct.
    zerofn = @(Intens) Intens - exp(-(t-3.*S.*z(n).*Intens).^2);

    figure(2);
    clf;
    subplot(3,1,1);
    ss=sseqn(t,z(n),S);
    
    plot(t,abs(ifft(fft(sqrt(exp(-t.^2))))).^2,'k--',t,abs(ifft(fft(u))).^2,t,abs(ifft(fft(ss))),'r--');
    xlabel('t/t_0','FontSize',15); ylabel('Intensity','FontSize',15);
    grid on;
    title(sprintf('S z/L_{nl} = %0.2f,  N_{samples} = %i',S*z(n),N(n)),'FontSize',15);
    axis([-3 3 0 1.1]);
    drawnow;
    subplot(3,1,2);
    semilogy(omega,abs(U).^2./max(abs(U).^2));
    xlabel('\omega-\omega_0','FontSize',15); ylabel('Relative spectral power','FontSize',15);
    grid on;
    title('Spectrum of self-steepening pulse (log scale)');
    axis([-100 100 1e-30 1]);
    subplot(3,1,3);
    plot(t,zerofn(abs(u).^2),t,zerofn(ss.'));
    axis on;
    axis([-3 3 -1E-6 1E-6]);
    legend('FFT','Semi-analytic');
    grid on;
    drawnow;
    
    tlast=toc;
end

if abs(toc-T) > 60*30
    fprintf('Quitting as run time has exceeded 30 minutes.\n');
    IsQuit=1;
end
end

%SSEQN Self-Steepening Equation this returns the self-steepening
%profile determined from the method of characteristics for a
%Gaussian initial profile.
function result = sseqn(t,z,s)
Intens = @(t,z,s) fzero(@(Izt) Izt - exp(-(t-3.*s.*Izt.*z).^2),0,optimset('TolFun',1E-22,'TolX',1E-22));
[t,z,s]=ndgrid(t,z,s);
for n=1:length(t)
  result(n)=Intens(t(n),z(n),s(n));
end
end

function N = IsAliasDown1_Custom(u,Options,dim) 
edge = floor(size(u,dim)/2)-fibonacci(ifibonacci(size(u,dim))-5);
N = isfftnalias(u,Options.TolAliased1,edge,dim);
end

% 
% $Log: gaffe_demo_steepen.m,v $
% Revision 1.8  2010/05/24 08:24:19  graceej
% * Improved documentation.
%
% Revision 1.7  2010/05/23 14:55:30  graceej
% * Did not mark z as persistent!
%
% Revision 1.6  2010/05/23 13:17:15  graceej
% * Added new demonstration and improved performance.  Corrected bug in callback and rationalised forms of callbacks.
%
% Revision 1.5  2010/05/22 20:58:27  graceej
% * Link help to gaffe_demo_soliton
%
% Revision 1.4  2009/10/24 11:08:05  graceej
% * Modified license to make use of the current BSD style open source initiative license.
%
% Revision 1.3  2009/05/06 20:09:45  graceej
% * Updated make package target to include distribution of the license file and README.
%
% Revision 1.2  2009/05/06 17:24:32  graceej
% * Update release date to some standard style.
% * Correct error in gaffe_demo_steepen
%
% Revision 1.1  2009/05/06 16:30:34  graceej
% * Moved demos to main branch.
%
% Revision 1.1  2009/05/06 16:27:28  graceej
% * Renamed self-steepening demo.
% * Removed unwanted alternative demos.
%
% Revision 1.1  2009/04/17 11:42:32  graceej
% * Brought over auto_checking/BPC-LIB tag=end and checked in.
% * This old repository is now defunct.
%
% Revision 1.2  2009/04/17 11:37:27  graceej
% * Wholescale modification of the entire library.
%
% Revision 1.1.2.6  2008/11/17 17:55:30  graceej
% * Improved error tolerance.
%
% Revision 1.1.2.5  2008/11/06 17:58:39  graceej
% * Improved documentation
% * Added command-line arguments.
%
% Revision 1.1.2.4  2008/11/06 16:59:11  graceej
% * Self-steepen default procedure.
%
% Revision 1.1.2.3  2008/09/16 13:52:17  graceej
% * Modified to select mesh resizing choice.
%
% Revision 1.1.2.2  2008/09/16 11:24:32  graceej
% * In addition plot spectrum.
%
% Revision 1.1.2.1  2008/09/16 09:43:01  graceej
% * Copied over from demo_self_steepening and renamed.
%
% Revision 1.1.2.2  2008/09/11 17:32:49  graceej
% * Fixed reporting of mesh size.
%
% Revision 1.1.2.1  2008/09/08 14:28:46  graceej
% * Modified to generalise the number of dimensions.
%
% Revision 1.1  2008/08/21 18:03:22  graceej
% * Comprehensively revised to reproduce Agrawal example.
%
% Revision 1.2  2008/08/18 16:30:23  graceej
% * Modified initial condition to alter power rather than height.
%
% Revision 1.1  2008/08/15 15:54:50  graceej
% * Renamed test_propation_engine.m to test_solvenlse2.m
%
% Revision 1.11  2008/08/15 15:51:09  graceej
% * Renamed propagation_engine.m to solvenlse2.m
%
% Revision 1.10  2008/08/15 15:39:15  graceej
% * Testing a symmetric initial condition on an asymmetrical grid.
%
% Revision 1.9  2008/08/15 13:33:19  graceej
% * Square symmetry still present, but should not be implicitly assumed in propagation_engine.
%
% Revision 1.8  2008/08/15 00:20:39  graceej
% * Return KX as well for convenience.
%
% Revision 1.7  2008/08/14 17:00:54  graceej
% * Setup sensible default propagation problem.
% * Add a measure of how fast the propagation is being carried out, we can see the slowdown as it approaches the singularity very clearly for height > 2.
%
% Revision 1.6  2008/08/13 19:09:12  graceej
% * Modified to use a more challenging initial condition and to plot the progress of a slice through the beam.
%
% Revision 1.5  2008/08/12 22:56:13  graceej
% * Modified propagation engine to resample the field before propagation to sit in the centre of a nyquist optimal window.
% * Check for aliasing after propagation and, if necessary, upsample the appropriate spaces.
% * Now behaves highly reliably and more robustly over wide conditions of accuracy.
%
% Revision 1.4  2008/07/23 18:50:27  graceej
% * Now demonstrates the propagation engine for a canonical Gaussian near collapse.
%
% Revision 1.3  2008/07/15 07:32:26  graceej
% * Demonstrate with a beam passing through focus
% * Closer to blowup
%
% Revision 1.2  2008/07/14 18:37:11  graceej
% * Correct minor timeout bug.
%
% Revision 1.1  2008/07/13 15:21:30  graceej
% * Commit of propagation engine and supporting code.
%
%
