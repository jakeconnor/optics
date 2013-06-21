%EVOLVE Evolve a generalised nonlinear wave equation.
%
%   Given an initial (sampled) field and its coordinate system, this code
%   efficiently evolves the signal (usually a pulse or beam) according to
%   the specified operators up to the maximum evolutionary value.
%
%   The initial field, assumed to be fft zero center ordered U0 and a
%   corresponding vector describing the coordinate system are (at minimum)
%   all that are required.  The propagation engine will then propagate this
%   field over the default distance returning the final field in U1, its
%   coordinate system vector (which may differ due to sampling changes) and the
%   Z distance.
%
%   If the initial field is a 1D vector the algorithm assumes you are
%   modeling the 1D NLSE.  This is typically used to model fibre
%   propagation.
%
%     [Z, U1, T1, BETA1] = EVOLVE(U0, T);
%
%   If the initial condition is a 2D matrix the algorithm assumes you wish
%   to model the 2D NLSE.  This is typically used to model nonlinear (Kerr
%   type) beam propagation in the paraxial approximation. If just a single
%   coordinate vector is supplied the field is assumed to have a
%   symmetrical coordinate system.
%
%     [Z,U1,X1,KX1,Y1,KY1] = EVOLVE(U0,X);
%
%   To propagate fields that are asymmetrical with respect to X and Y,
%   either because they are inherently asymmetrical or they are sampled on
%   an asymmetric grid you should supply the coordinate vector for the
%   other dimension (Y).
%
%     [Z,U1,X1,KX1,Y1,KY1] = EVOLVE(U0,X,Y);
%
%   You can alter the internal behaviour by modifying the options
%   structure, described below and passing it as the last argument, either
%   by calling with.
%
%     [Z,U1,X1,KX1,Y1,KY1] = EVOLVE(U0,X,OPTIONS);
%
%   or, for asymmetrical fields.
%
%     [Z,U1,X1,KX1,Y1,KY1] = EVOLVE(U0,X,Y,OPTIONS);
%
%
%OPTIONS STRUCTURE
%
%   The default options for the engine can be determined by calling
%   SOLVENLSE2 with the single argument 'defaults', which will return a
%   structure containing the default options.
%
%     OPTIONS = EVOLVE('defaults');
%
%   The option fields are described below.
%
%     TolGlobal          The global error tolerance to aim for.  This is a
%                        measure of how good the 4th order split-step
%                        method approximates commutation of the linear and
%                        nonlinear operators.
%
%     TolAliased         A measure of how closely we approximate the
%                        finitely represented function before we consider
%                        it to be aliased, either in real space or in
%                        frequency space.
%
%     FracRAliased       When more a fraction of amplitude TolAliased
%                        greater than the peak field amplitude or a
%                        fraction of power greater than TolAliased is found
%                        outside of the normalised radius R, the field is
%                        considered to be aliased.
%
%                        When the field is aliased the grid size is
%                        increased.
%
%     FracROverSpecified When the field is well contained (not considered
%                        to be aliased) within a small fraction of the
%                        space it is considered to be over specified.  When
%                        this occurs the grid is reduced to accomodate the
%                        new field.
%
%     MaxZ               The maximum propagation distance to allow.
%
%     MinDZ              The minimum z step to allow.
%
%     MaxDZ              The maximum z step to allow. If the z step that
%                        is permitted by the global error exceeds this
%                        value it dz is clamped to MaxDZ.
%
%     InitialDZ          The first z step to try.
%
%     MaxIterations      The maximum number of iterations to allow.
%
%     MaxN               The maximum size of a grid in samples along one
%                        edge, i.e. MaxN^2 total points.
%
%     Trace              If set on prints a verbose trace of the progress
%                        of the propagation engine.
%
%     Callback           A structure consisting of function handles that
%                        are called at certain points during the
%                        propagation process.
%
% CALLBACK STRUCTURE
%
%    The callback option is a structure of callback functions which should
%    take no arguments.  By default the called functions should use
%    'evalin' to access the variables in the caller workspace (inside
%    EVOLVE).
%
%      BeginIteration    Called at the beginning of an iteration.
%
%      EndIteration      Called at the end of an iteration.
%
%
%See also: EVOLVEMESHSET, DefaultEndIteration, GAFFE_DEMO_LINEAR,
%GAFFE_DEMO_SOLITON, GAFFE_DEMO_STEEPEN, GAFFE_DEMO_GAUSSIAN

% $Author: graceej $ $Date: 2010/05/24 08:23:06 $
% $Revision: 1.10 $
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

function [z,u1,X,KX]=evolve(varargin)

% Initial fudge value.  This should be expunged properly.  Its purpose is
% to prevent trying to down-size a field that has been upscaled too soon.
%
% In practice we should have a more careful set of checks that defaults to
% the virgin field (before being propagated).
%Fudge1 = 2; % Moved to options.

%% Setup options
% Define the default options, these are options which appear reasonable for
% general problems that I have encountered. They are by no means
% authoritative.
% FracROverSpecified is (1 - 1/G)^2 where G is the golden ratio.
G=(1+sqrt(5))/2;
defaultopt = struct(...
    'TolGlobal',1E-6,...
    'TolAliased1',1E-8,...
    'TolAliased2',1E-9,...
    'Fudge1',2,...
    'FracRAliased',1/G,...
    'FracROverSpecified',1/G/G/G,...
    'MaxZ',100,...
    'MinDZ',10^ceil(log10(eps*1E4)),...
    'MaxDZ',10,...
    'InitialDZ',[],...
    'MaxIterations',1E4,...
    'MaxN',6718464,...
    'Trace',0,...
    'User1',[],...
    'Callback',struct(...
    'BeginIteration',[],...
    'EndIteration',[],...
    'OperatorLinear',@DefaultDispersion,...
    'OperatorNonlinear',@DefaultKerr,...
    'EvolveStep',@DefaultEvolveStep,...
    'Mesh',evolvemeshset('default')...
    )...
    );
% If we are called with just "defaults" and are expected to return just one argument,
% return the default options.  This is the format expected by 'optimget'.
if nargin==1 && nargout <= 1 && isequal(varargin{1},'defaults')
    z = defaultopt;
    return
end


%% Parse calling arguments
options = struct();
% Obtain initial field.
u0 = varargin{1};
% Make sure that the field is square (same number of samples along each
% edge) and set the size variable N.
N = size(u0);
% Check that there is no leading singleton dimension.
if N(1) == 1
    throw(MException('SanityCheck:LeadingSingletonDimension','Please reshape the input field and coordinates to avoid leading singleton dimensions'));
end
% Determine number of dimensions (a 1D signal has two dimenions, one is a
% singleton).
NDim = ndims(u0);
% Get a list of non singleton dimensions
NonSingletonDimensions = find(N>1);
x = cell(NDim,1);
for idx = 2:nargin
    if isstruct(varargin{idx})
        options = varargin{idx};
    else
        x{NonSingletonDimensions(idx-1)}=varargin{idx};
    end
end
for d = 1:length(x)
    if isempty(x{d})
        x{d} = [0];
    end
end




%% Merge options with existing options.
o=structmerge(defaultopt,options);

%% Sanity check options.
if o.FracROverSpecified > 1 || o.FracROverSpecified < 0
    throw(MException('SanityCheck:FracROverSpecified','FracROverSpecified should be < 1.0 and > 0'));
end
if o.FracRAliased > 1 || o.FracRAliased < 0
    throw(MException('SanityCheck:FracRAliased','FracRAliased should be < 1 and > 0.'));
end
if o.FracROverSpecified > o.FracRAliased
    throw(MException('SanityCheck:FracROverSpecified','FracROverSpecified should be less than FracRAliased'));
end
if ~isempty(o.InitialDZ) && o.InitialDZ < o.MinDZ
    throw(MException('SanityCheck:InitialDZ','InitialDZ should be greater than MinDZ'));
end

% Check callbacks.  If they are not empty check if they are callback
% functions.  If not, then complain!
if ~isempty(o.Callback.BeginIteration) && ~isa(o.Callback.BeginIteration,'function_handle')
    throw(MException('SanityCheck:CallbackNotFunction','Callback.BeginIteration should be a function handle.'));
end
if ~isempty(o.Callback.EndIteration) && ~isa(o.Callback.EndIteration,'function_handle')
    throw(MException('SanityCheck:CallbackNotFunction','Callback.EndIteration should be a function handle.'));
end




if prod(N)^(1/NDim) > o.MaxN
    throw(MException('SanityCheck:FieldTooBig',sprintf('The field size exceeds the maximum dimension MaxN=%i',MaxN')));
end

%% Setup supporting coordinate systems and transforms.
% U0 should always contain the most recent fft of the field.
U0 = fftn(u0);


%% Initialise the field according to the initialisation callback.
N_new = o.Callback.Mesh.SetInitialSize(N);
if any(N_new ~= N)
    u0 = fftnpad(u0,N_new);
    U0 = fftn(u0);
    if o.Trace
        fprintf('Resizing initial field from %s to %s.\n',sprintsize(N),sprintsize(N_new));
    end
end
for d = 1:NDim
    [x{d},dx{d},kx{d},dkx{d}] = fftspace(x{d}*(N_new(d)/N(d)),N_new(d));
end
N = N_new;

%% Grid the spaces
X = cell(size(x)); KX = cell(size(x));
[X{:}]  = ndgrid(x{:});
[KX{:}] = ndgrid(kx{:});
N_mesh = N;

%% Initialise the callback to be executed at the beginning and end of each iteration.



% Call the begin iteration callback, if it is defined.
if ~isempty(o.Callback.BeginIteration)
    is_terminated_by_callback = feval(o.Callback.BeginIteration,'begin');
end

% The callback should, by default, expect no arguments.  If there is an
% argument called 'init' the callback should internally clear all its
% persistent variables.
if ~isempty(o.Callback.EndIteration)
    feval(o.Callback.EndIteration,'begin');
end

%% Select initial dz.
if isempty(o.InitialDZ)
    % The step size should be the minimum of a step size that offers a
    % maximum nonlinear phase shift of (say) pi/8 and a step size that offers a
    % maximum linear phase shift in the spectrum of (say) pi/4 or a tenth of the
    % total maximum propagation distance (to make it worthwhile simulating
    % something).
    
    % Step size based on sampling maximum z distance.
    dz_max_Z = o.MaxZ/10.0;
    
    % Step size if we assume a maximum nonlinear phase shift.
    dz_max_NL = (pi/8)./abs(maximumvalue(u0)).^2;
    
    % Step size based on maximum spectral phase rotations.
    dz_max_L = dz_max_Z*ones(1,NDim);
    for d = NonSingletonDimensions
        dz_max_L(d) = 2*(pi/4)./(width(U0,KX{d})/2);
    end
    
    % Initial dz should be the minumum of these.
    dz = min([dz_max_NL dz_max_L dz_max_Z]);
    
    % If this dz is too small then the initial dz should be set to several
    % orders of magnitude more than the minimum allowed step size.
    dz=max([dz o.MinDZ*100]);
    
    % Since we know how the error should behave as we scale dz, assuming we
    % are in a sensible region, it should be possible to rapidly find the
    % ideal step size by simple root finding.
    % @bug Following function should probably involve a reentrant call to
    % solvenlse.
    if o.Trace
        fprintf('Initial estimate for dz=%e, refining.\n',dz);
    end
    % Root find to determine the dz for an error that is half the tolerated
    % global error.
    dz = estimate_dz(o.TolGlobal*0.5,u0,U0,dz,X,KX,o.Callback.EvolveStep,o.Callback.OperatorLinear,o.Callback.OperatorNonlinear);
    if o.Trace
        fprintf('Initial stepsize refined to dz=%e.\n',dz);
    end
else
    dz = o.InitialDZ;
end


%% Main loop.
%
% We iterate until any of the conditions for exit have been met, for
% example exceeding the maximum number of iterations, maximum sample size
% etc.
n = 0; % Iteration counter.
z = 0; % Current z position @bug Calculated by addition of floating point.
% Size of the resized, new field.
is_terminated_by_callback = 0; % Flag to indicate that a callback function requests exit.
% True for each dimension of the field that requires enlarging in real
% space.
is_too_big_x = zeros(1,NDim);
% True for each dimension of the field that requires enlarging in real
% space.
is_too_big_kx = zeros(1,NDim);
% True for each dimension of the field that requires shrinking in
% real space.
is_too_small_x = zeros(1,NDim);
% True for each dimension of th efield that requires shrinking in frequency
% space.
is_too_small_kx = zeros(1,NDim);
% The step length of the last propagation.
dz_try=dz;
dz_next=dz;
dz_last=dz_try;


% Briefly the algorithm works like this.
%
% Check to see if the field is over represented in real space or frequency
% space and did not end up aliased because of the last attempted propagation.
% If it is over specified, then reduce the number of samples appropriately and
% recalculate the associated coordinate systems.
%
% Propagate the field over the current chosen test step length dz_try.
%
% If the field ends up aliased in real or k space, i.e. it's expanded too
% much because of propagation or too much in k space due to large
% nonlinearity then enlarge the mesh and retry the propagation.
%

% Force N_mesh to be zero so that the mesh is made.
N_mesh = 0*N;
n_timeout=0;
% Used to control the restarting of an iteration.
is_retry_iteration=0;
is_final_step=0;
is_ended=0;
while n < o.MaxIterations && ...
        z < o.MaxZ && ...
        ~is_terminated_by_callback && ...
        ~is_ended
    % Start this iteration by assuming we won't retry the iteration.
    is_retry_iteration=0;
    
    if(prod(N) > o.MaxN)
        throw(MException('SolveNLSE:Bailout','Grid size exceeds maximum'));
    end
    if (dz_try < o.MinDZ & ~is_final_step)
        throw(MException('SolveNLSE:Bailout','Step size is smaller than minimum'));
    end
    
    % Call the begin iteration callback, if it is defined.
    if ~isempty(o.Callback.BeginIteration)
        is_terminated_by_callback = feval(o.Callback.BeginIteration);
    end
    
    if o.Trace
        fprintf('Start attempt of step %d.\n',n);
    end
    
    % Only allow attempted signal truncation after a certain number of
    % iterations since the last signal truncation.
    % @bug This is a bug -- there is clearly a complex and subtle interplay
    % between the noise floor and signal resizing when up/down scaling.
    if n >= n_timeout
        % Check to see if the field is too narrow in real space and truncate it.
        % This is only attempted if previously it was not too big in real
        % space.
        if any(~is_too_big_x)
            for d = NonSingletonDimensions
                if ~is_too_big_x(d)
                    [u0, N_new(d)] = truncate_field(u0, o, d);
                    [kx{d},dkx{d},x{d},dx{d}] = fftspace(kx{d},N_new(d));
                end
            end
            if any(N_new ~= N)
                U0 = fftn(u0);
                is_too_small_x = N_new ~= N;
                N=N_new;
            end
        end
        % Repeat for spatial frequency representation.
        if any(~is_too_big_kx)
            for d = NonSingletonDimensions
                if ~is_too_big_kx(d)
                    [U0, N_new(d)] = truncate_field(U0, o, d);
                    [x{d},dx{d},kx{d},dkx{d}] = fftspace(x{d},N_new(d));
                end
            end
            if any(N_new ~= N)
                % N.B. MATLAB fft does not scale by 1/N.  By default we must scale
                % the spectral power to conserve energy.
                U0 = U0*prod(N_new./N);
                u0 = ifftn(U0);
                is_too_small_kx = N_new ~= N;
                N=N_new;
            end
        end
        n_timeout = n + o.Fudge1;
    end
    
    % Check to see if the new field size differs from the previous field
    % size.  If so, regenerate the coordinate meshes.
    if any(N ~= N_mesh)
        if o.Trace
            fprintf('Field size has changed from %s to %s, recalculating coordinate meshes.\n',sprintsize(N_mesh),sprintsize(N));
        end
        [X{:}] = ndgrid(x{:});
        [KX{:}] = ndgrid(kx{:});
        N_mesh=N;
    end
    
    % Check to see if the trial propagation distance dz_try exceeds the
    % maximum dz we have permitted (if any).  If so then cap it to the
    % upper limit. length is borked too
    if dz_try > o.MaxDZ
        if o.Trace
            fprintf('Trial dz step (%e) exceeds MaxDZ (%e) capping to upper limit. \n',dz_try,o.MaxDZ);
        end
        dz_try = min(dz_try,o.MaxDZ);
    end
    
    % Propagate forwards one step of length dz and determine error
    % information using the supplied forward propagator.
    [local_error,u1,U1,uc,UC,uf,UF,U0] = ...
        o.Callback.EvolveStep(u0, U0, dz_try, X, KX, ...
        o.Callback.OperatorLinear, ...
        o.Callback.OperatorNonlinear);
    
    
   % if mod(n,50)==0 %samples every 50 steps
        global j waves initialsize;
        if n ==0
            waves.u= zeros(size(u0),100000);
            waves.tau=x{1,1};
            waves.z= zeros(1,100000);
            initialsize=size(u0,1);
        end
        if size(u1,1)==initialsize
            waves.u(:,(n)+1)=u1;
            waves.z(:,(n)+1)=z;
        else
            waves.u(:,(n)+1)=interp1(u1,1:size(u1,1)/initialsize:size(u1,1),'spline');
            waves.z(:,(n)+1)=z;
        end
        %for j=1:n
  %   end
%         try %%this will save the waveforms at their original size
                % too much computation time, so saving them at lower
                % resolution instead
%             waves{j}.u(:,n+1)=u1;
%             waves{j}.tau(:,n+1)=x{1,1};
%             waves{j}.z(:,n+1)=z;
%         catch
%             %truncate matrices
%             waves{j}.u(:,~any(waves{j}.u,1))=[];
%             waves{j}.tau(:,~any(waves{j}.tau,1))=[];
%             waves{j}.z=waves{j}.z(1:size(waves{j}.u,2));
%             j=j+1;
%             try 
%                 size(waves{j}.u); %check if new matrices need pre-allocation
%             catch
%                 pack; %frees up memory, kind of.
%                 waves{j}.u= zeros(size(u0),10000);%and pre allocate
%                 waves{j}.tau= zeros(size(u0),10000);
%                 waves{j}.z= zeros(1,10000);    
%             end
%             waves{j}.u(:,n+1)=u1;
%             waves{j}.tau(:,n+1)=x{1,1};
%             waves{j}.z(:,n+1)=z;
%         end

    
    if o.Trace
        fprintf('Trial propagation over a distance dz=%0.4e \n',dz_try);
    end
    
    % If TolGlobal is NaN then the step size is not allowed to dynamically
    % change.
    if ~isnan(o.TolGlobal)
        % Determine the next step size to use depending on the current error
        % and global error tolerance.
        [dz_next, isok] = select_dz(dz_try, local_error, o.TolGlobal);
        % Check to see if we are within the error bound.  If not, discard this
        % step.
        % @bug The error could also be because the mesh is too small.
        if ~isok
            if o.Trace
                fprintf('Propagation step exceeds error limit, local error (%e) > TolGlobal (%e), discarding step.\n',local_error,o.TolGlobal);
            end
            % Try the new (smaller) step size.
            if dz_next > dz_try
                throw(MException('SanityCheck:StepSize','Next expected propagation step (%e) is larger than previous trial (%e).'),dz_next,dz_try);
            end
            dz_try = dz_next;
            % Make sure that we retry the iteration.  Don't do it immediately
            % however, first check to see if there are any other problems.
            is_retry_iteration = 1;
        else
            % The error bound was acceptable, report on new step size if it
            % differs.
            if o.Trace && dz_next ~= dz_try
                if dz_next > dz_last
                    fprintf('Increasing next step size to %0.4e .\n',dz_next);
                else
                    fprintf('Decreasing next step size to %0.4e .\n',dz_next);
                end
            end
        end
    end
    
    % Test to see if the field is aliased as a result of the propagation
    % step. If it is then we wish to rescale the initial field.
    for d = NonSingletonDimensions
        is_too_big_x(d)  = o.Callback.Mesh.IsAliasDown1(u1, o, d);
        is_too_big_kx(d) = o.Callback.Mesh.IsAliasDown1(U1, o, d);
    end
    
    
    % If the field is aliased in real space AND k space, nothing can be
    % done, the field is rubbish.  Since the field before starting
    % propagation should be fine the reason the field is now aliased in
    % real space and k space is because the propagation has been too
    % extreme.  An example would be something that has experienced
    % significant linear and nonlinear effects.  The only way to avoid this
    % without enlarging the spaces is to propagate over a smaller step.
    %
    % @bug Perhaps we should also increase the spectral and spatial
    % representation.
    if any(is_too_big_x & is_too_big_kx)
        if o.Trace
            fprintf('Field aliased in both real space and k space, reducing step size.\n');
        end
        dz_try = dz_try/2;
        is_retry_iteration = 1;
        continue
    end
    
    
    %
    big_small_x = is_too_big_x & is_too_small_x;
    big_small_kx = is_too_big_kx & is_too_small_kx;
    
    
    % If the field has, for example, propagated a long way it will have expanded outside
    % the boundary of the computational window.  If this is the case we
    % should enlarge the real space (interpolate k space) and try
    % propagating again.
    if any(is_too_big_x) % Only test if there is a dimension that's aliased.
        for d = NonSingletonDimensions
            if is_too_big_x(d)
                % Get the size of the newly upscaled field along this dimension.
                N_new(d) = o.Callback.Mesh.GetUpSize(N(d),o);
                % Adjust the coordinate systems for this dimension based on
                % the new sample size.
                [kx{d},dkx{d},x{d},dx{d}] = fftspace(kx{d},N_new(d));
                if o.Trace
                    fprintf('Propagation caused aliasing in real space along dimension %d.  Increasing x space from %s to %s samples.\n',d,sprintsize(N),sprintsize(N_new));
                end
            end
        end
        if any(N_new ~= N) % @bug this should be redundent.
            % Pad the signal based upon the new dimensions.
            u0 = fftnpad(u0,N_new);
            % Make sure the current transform does not lie!
            U0 = fftn(u0); % @bug This should not be strictly needed.
            % Retry the propagation.
            N = N_new;
            is_retry_iteration = 1;
        end
    end
    
    % If the field has, for example, experienced a lot of nonlinearity its spectrum can
    % expand outside the computational window.  As with the real space
    % window this means we should expand the spectral representation
    % (interpolate real space).
    if any(is_too_big_kx)
        for d = NonSingletonDimensions
            if is_too_big_kx(d)
                N_new(d) = o.Callback.Mesh.GetUpSize(N(d),o);
                [x{d},dx{d},kx{d},dkx{d}] = fftspace(x{d},N_new(d));
                if o.Trace
                    fprintf('Propagation caused aliasing in frequency space along dimension %d.  Increasing kx space from %s to %s samples.\n',d,sprintsize(N),sprintsize(N_new));
                end
            end
        end
        if any(N_new ~= N) % @bug this should be redundent.
            % N.B. Scale by the number of samples to conserve power.
            U0 = fftnpad(U0,N_new)*prod(N_new./N);
            u0 = ifftn(U0);
            N = N_new;
            is_retry_iteration=1;
        end
    end
    
    % Determine if we should retry the iteration.  If so, use the samller
    % of the current step and the suggested next step.
    if is_retry_iteration
        dz_try = min([dz_try dz_next]);
        continue;
    end
    
    
    % At this point we have successfully completed an iteration, stepping
    % forward a distance dz_try.
    z = z + dz_try;
    
    
    % Call the begin iteration callback, if it is defined.
    if ~isempty(o.Callback.EndIteration)
        is_terminated_by_callback = feval(o.Callback.EndIteration);
    end
    
    % The propagated field u1 should now be propagated forward.  So the
    % input field becomes the previous output field.
    u0 = u1;
    U0 = U1;
    
    
    % Check to see if the next step over a distance of dz_next will take us
    % beyond the maximum propagation distance of the simulation. If so we
    % override the next step length to make sure it takes us exactly to the
    % end.
    if is_final_step
        is_ended=1;
    elseif z + dz_next > o.MaxZ
        % If we are close to machine precision o.MaxZ - z could reverse
        % sign and unwittingly cause us to try and propagate in reverse.
        % Formally make certain it has a consistent sign.
        sn = sign(o.MaxZ);
        if ~sn; sn=1; end
        dz_try = sn*abs(o.MaxZ - z);
        is_final_step=1;
        if o.Trace
            fprintf('Overriding step length for final step to %e.\n',dz_try);
        end
    else
        dz_last = dz_try;
        dz_try = dz_next;
    end
    
    if o.Trace
        fprintf('Iteration %i, z=%f complete.\n\n',n,z);
    end
    
    
    % Increment iteration counter.
    n = n + 1;
end

%% Finalise the begin / end of iteration callback if it exists.


% Call the begin iteration callback, if it is defined.
if ~isempty(o.Callback.BeginIteration)
    is_terminated_by_callback = feval(o.Callback.BeginIteration,'end');
end

% Call the begin iteration callback, if it is defined.
%
% This gives the callback an opportunity to finalise its internal state so
% that when it is called with 'get' it can return what is wanted.  For
% example the callback may have preallocated over-sized buffers -- this
% gives it an opportunity to resize them correctly.

if ~isempty(o.Callback.EndIteration)
    feval(o.Callback.EndIteration,'end');
end

end



%TRUNCATE_FIELD Optionally truncate a representation of a field.
%
%   This helper function checks to see if a field is aliased, if not it
%   attempts to truncate the field so that it fits in between being too
%   small (over specified) and too large (aliased).
%
function [u,N_new] = truncate_field(u,o,dim)
N=size(u,dim);
N_new=N;

% Check to see if we are working along a trivial (singleton) dimension.  If
% so just exit.
if N == 1
    too_big = 0; too_small = 0; u=u;
    return
end


% If the field before we try to propagate it forwards is aliased
% (by this measure) we have irrevocably lost information and
% subsequent propagation will produce rubbish -- so we *must* bail
% out at this point.
is_too_big = o.Callback.Mesh.IsAliasDown1(u, o, dim);
if is_too_big
    throw(MException('SanityCheck:FieldAliased',sprintf('The initial representation appears to be aliased along dimension %d.',dim)));
end

% Check to see if the field is over represented.  This is equivalent to it
% not being aliased on the smaller grid.
is_too_small  = ~o.Callback.Mesh.IsAliasDown2(u, o, dim);
if is_too_small
    if o.Trace
        fprintf('Field representation over specified, attempting truncation.\n');
    end
    %Determine the new number of samples to use for the function u.
    N_new = o.Callback.Mesh.GetDownSize(size(u,dim),o);
    
    % Under certain circumstances it is possible for N_new > N i.e.
    % due to rounding error when we have a small number of samples. It is
    % dependent on the (unknown) callback so anything is possible.
    %
    % We check to see if the proposed number of samples *is* smaller.  If
    % so we down sample.  If not we leave things as they are.
    if N_new < N
        % Pad u1, since N < size(u) it will "unpad", remove a section from
        % the centre.
        N_all = size(u);
        N_all(dim) = N_new;
        u1 = fftnpad(u, N_all);
        
        % Truncating the field, effectively filtering it with a
        % rectangular window, will raise the noise floor (wings)
        % of the associated domain.  We should therefore check that
        % truncating the field as we intend does not end up causing
        % aliasing in the other domain
        %
        % N.B. it could be either
        % fft2 or ifft2, it doesn't matter as we are only interested in the
        % magnitude of the result.
        U1 = fftn(u1);
        if o.Callback.Mesh.IsAliasDown1(U1, o, dim)
            if o.Trace
                fprintf('Truncated representation is aliased in complement space along dimension %d, reverting field.\n',dim);
            end
            N_new = N;
        else
            u = u1;
        end
    else
        N_new = N;
        if o.Trace
            fprintf('Truncating representation would not significantly reduce field size, reverted.\n');
        end
    end
end
end

%SPRINTSIZE Given a size vector generate a string representing the size.
function S = sprintsize(SZ)
S = sprintf('%dx',SZ(1:end-1));
S = [S sprintf('%d',SZ(end))];
end

%MAXIMUMVALUE Given an N dimensional array, determine the maximum value.
function V = maximumvalue(A)
ND = ndims(A);
V = A;
for dim = 1:ND
    V = max(V);
end
end


%SELECT_DZ Determine the stepsize and if it is within the error bound.
%
%    Given the prevailing step size DZ, the local error DELTA and the
%    required target global error DELTA_G determine the new stepsize and if
%    it is suitable for use in the next propagation step.
%
%    [DZ_NEW,ISOK] = SELECT_DZ(DZ,DELTA,DELTA_G);
%
%    When ISOK is zero (false) the returned step size is too large and any
%    result obtained propagating over that distance should be discarded.
function [dz,isok] = select_dz(dz,delta,delta_G)
isok=1;
% Determine new step size based on current step size and error targets.
if delta >= 2*delta_G
    dz=dz/2.0;
    isok=0; % Not ok, should discard solution.
elseif delta > delta_G && delta < 2*delta_G
    dz=dz/2.0^(1/3);
elseif delta < 1/2*delta_G
    dz=dz*2.0^(1/3);
else
    dz=dz;
end
end


%ESTIMATE_DZ Determine an initial dz based on 4th order error behaviour.
%
%   Given an initial field, its Fourier representation and coordinate
%   system and a (reasonable) estimate for dz functions describing the
%   linear and nonlinear operators and the associated coordinate systems,
%   this function attempts to determine a value of dz that satisfies the
%   target error.
%
%   Providing the initial guess is within the range of dz values over which
%   we expect to see cubic error behavior this works well.  If we are
%   unwittingly outside this area the gradient of the line is likely to be
%   significantly different from 3 and the goodness of fit (standard
%   deviation of residuals) will be poor.
%
%   We then pick the dz which produces an error that is closest to the
%   target error as the staring point.
%
%     DZ = ESTIMATE_DZ(TARGET_ERROR,U0, FFT_U0, DZ_0,X,KX, funLinear, funNonlinear);
function [dz, actual_error] = estimate_dz(target_error, u0, U0, dz, ...
    X, ...
    KX, ...
    funStep, ...
    funLinear, ...
    funNonlinear)

%% Secant root finding method to obtain a good estimate of dz.

% Since we expect the error properties of the propagation technique to
% follow a cubic relationship with respect to dz (linear on log log plot)
% by root finding for the log(dz) that gives the log(target_error) we can
% expect to find the appropriate dz very quickly.  The secant method is
% highly effective for functions with smooth first and second derivatives
% and a simple single root.  We expect the error curve to meet both these
% criteria.
%
% Typically this determines dz to better than 0.1% within 2 iterations.
%
% When the method fails (which can occur in the case of no nonlinearity) the initial dz is used instead.
x_prior = log(dz);
x_now   = log(dz./2);

f = @(zeta)  log(funStep(u0, U0, exp(zeta), X, KX, funLinear, funNonlinear)) - log(target_error);

% Bootstrap method
f_prior = f(x_prior);
f_now   = f(x_now);

% First test
err_first = f_prior.^2;

% We expect the root to be determined to the precision below within typically 3 iterations.
% If it requires more than 6 iterations we should bail out -- we are probably in
% a pathalogical region.
TolRoot = 1E-9;
for n=1:6
    x_next = x_now - (x_now - x_prior)./(f_now - f_prior) * f_now;
    f_prior = f_now;
    x_prior = x_now;
    x_now = x_next;
    f_now = f(x_now);
    % Convergence test, if we are calculating rubbish (e.g. NaN Inf) then exit.
    if ~isfinite(f_now)
        break
        % If the error is less than the tolerance and better than the
        % initial guess then we have found the root.
    elseif f_now.^2 < TolRoot && f_now.^2 < err_first
        dz = exp(x_now);
        break
    end
end
end



%% CVS log information
%
% $Log: evolve.m,v $
% Revision 1.10  2010/05/24 08:23:06  graceej
% * Improved documentation.
%
% Revision 1.9  2010/05/23 13:17:15  graceej
% * Added new demonstration and improved performance.  Corrected bug in callback and rationalised forms of callbacks.
%
% Revision 1.8  2010/05/22 20:50:06  graceej
% * Fixed a bug where the spectrum wasn't being updated from the previous run.
%
% Revision 1.7  2009/10/24 11:04:53  graceej
% * Modified to BSD license.
% * Modified so that a value of NaN for the TolGlobal prevents the dynamic step sizing from being carried out.
%
% Revision 1.6  2009/05/06 20:09:45  graceej
% * Updated make package target to include distribution of the license file and README.
%
% Revision 1.5  2009/05/06 16:30:34  graceej
% * Moved demos to main branch.
%
% Revision 1.4  2009/05/06 09:13:57  graceej
% * Use fftspace instead of makefftspace.
%
% Revision 1.3  2009/04/21 12:36:09  graceej
% * Corrected minor typo in initial grid sizing.
% * TODO Modify so that the input coordinate systems are simple linear vectors rather than of the same dimension as u.
%
% Revision 1.2  2009/04/20 17:38:09  graceej
% * Simplify treatment of final step.
%
% Revision 1.1  2009/04/17 11:42:32  graceej
% * Brought over auto_checking/BPC-LIB tag=end and checked in.
% * This old repository is now defunct.
%
% Revision 1.2  2009/04/17 11:37:25  graceej
% * Wholescale modification of the entire library.
%
% Revision 1.1.2.8  2008/11/06 17:41:12  graceej
% * Removed breakpoint line.
%
% Revision 1.1.2.7  2008/11/06 15:31:29  graceej
% * Added (fudge) User1 option.  This is undocumented.
%
% Revision 1.1.2.6  2008/09/16 15:35:10  graceej
% * Reduce default maximum number of iterations.
%
% Revision 1.1.2.5  2008/09/16 14:53:08  graceej
% * Forward integration carried out by callback.
% * Callback is by default the standard forward propagation.
%
% Revision 1.1.2.4  2008/09/16 14:37:10  graceej
% * Changed function names to reflect new names.
% * Fixed function end error.
%
% Revision 1.1.2.3  2008/09/16 14:32:18  graceej
% * Moved estimate_dz to local function.
%
% Revision 1.1.2.2  2008/09/16 14:30:59  graceej
% * Moved select_dz to local function.
%
% Revision 1.1.2.1  2008/09/16 14:28:47  graceej
% * Initial checkin.
%
% Revision 1.5.2.10  2008/09/15 18:18:42  graceej
% * Made mesh resizing / stepsize selection more robust.
% * Prevented bailout if minimum stepsize reached on last iteration.
% * Increased mesh maximum.
%
% Revision 1.5.2.9  2008/09/11 17:30:32  graceej
% * Asymmetric precision requirements for up/down scaling.
% * Uses humble fibonacci number mesh as default.
% * Fudge1 used to prevent mesh oscillation for at least 2 iterations.  Non ideal -- but appears to work well.
%
% Revision 1.5.2.8  2008/09/11 16:45:33  graceej
% * Added Fudge1 option to options.  Not ideal, but need to think more carefully about how to handle the case of oscillating mesh sizes.
%
% Revision 1.5.2.7  2008/09/11 14:39:03  graceej
% * Intermediate modification - intend to solve problem of oscillating mesh size.
%
% Revision 1.5.2.6  2008/09/09 14:17:17  graceej
% * Removed legacy code and improved error messages.
%
% Revision 1.5.2.5  2008/09/08 14:28:46  graceej
% * Modified to generalise the number of dimensions.
%
% Revision 1.5.2.4  2008/09/05 20:41:58  graceej
% * Use fftpadn instead of obsolete pad2.  TODO - modify pad2, to operate along a given dimension.
%
% Revision 1.5.2.3  2008/09/02 13:52:37  graceej
% * Modified to use power of two length signals by default.
%
% Revision 1.5.2.2  2008/08/29 14:00:37  graceej
% * Using Fibonacci length signals.
%
% Revision 1.5.2.1  2008/08/28 11:17:42  graceej
% * Initial changes, use Fibonacci number (golden) ratios for determining if things are or are not aliased.
%
% Revision 1.5  2008/08/22 13:38:20  graceej
% * Fixed minor logical error in final step.
%
% Revision 1.4  2008/08/21 17:21:46  graceej
% * By default assume we are considering the 1D NLSE and use dispersion as the default linear operator.
%
% Revision 1.3  2008/08/20 17:18:42  graceej
% * Removed explicit IsNonlinear check.  This is managed by setting Callback.OperatorNonlinear to @DefaultIdentity.
% * Replaced default propagation kernels with arbitrary linear and nonlinear function callbacks.
% * Modified to use improved initial dz function for arbitrary operators.
% * Call "step" instead of "propagate4" to carry out evolution with arbitrary operator functions.
% * Moved caching of linear operators to "step" function.
% * Corrected minor bug in final step at z = o.MaxZ
% * Corrected bug when aliasing occurs, reverting step size to the previous step size.
%
% Revision 1.2  2008/08/19 15:23:54  graceej
% * Corrected output arguments.
% * Corrected help documentation.
%
% Revision 1.1  2008/08/19 14:09:13  graceej
% * Modification of solvenlse2 to allow 1D signals.
% * Caught a few minor bugs.
%
% Revision 1.3  2008/08/16 16:33:28  graceej
% * Bug fix when aborting signal truncation.
%
% Revision 1.2  2008/08/15 21:50:02  graceej
% * Modified to use new isaliased which relies purely on power.
% * Modified to use alias_edge which determines the edge of aliasing using purely power.
% * Changed alias tolerance to reflect useage of power (approximately squares requried tolearance).
%
% Revision 1.1  2008/08/15 15:51:08  graceej
% * Renamed propagation_engine.m to solvenlse2.m
%
% Revision 1.22  2008/08/15 15:48:56  graceej
% * Changed function name from propagation_engine to solvenlse2.
% * Updated help.
%
% Revision 1.21  2008/08/15 15:40:28  graceej
% * Make command line backward compatible.
% * Fix bug in grid down sizing.
% * Make alias tolerance a little tighter by default.
%
% Revision 1.20  2008/08/15 13:33:19  graceej
% * Square symmetry still present, but should not be implicitly assumed in propagation_engine.
%
% Revision 1.19  2008/08/15 13:28:07  graceej
% * Modified to break assumption of a square grid.
%
% Revision 1.18  2008/08/15 13:09:40  graceej
% * First set of changes for extension to asymmetric mesh.
% * Explicitly add other dimension, but use Ny everywhere.
% * Currently do not consider x, will change at next checkin.
%
% Revision 1.17  2008/08/15 12:01:34  graceej
% * Operate assuming the N we are refering to is Ny, to make it easier to convert to  the asymmetric case.
% * Moved helper functions out to seperate M files as they will be the same for the asymmetric case.
%
% Revision 1.16  2008/08/15 00:21:15  graceej
% * Return KX as well for convenience.
%
% Revision 1.15  2008/08/14 23:57:24  graceej
% * Aliasing or unaliasing limits are set to be related to the golden ratio.
%
% Revision 1.14  2008/08/14 23:29:08  graceej
% * Fixed edge case in final dz when close to ZMax.
% * Added function for improving the estimate of the step size to use.
%
% Revision 1.13  2008/08/14 21:14:22  graceej
% * Made default initial conditions more conservative.
% * Fixed minor bug in final step selection for completion of the step to Zmax
%
% Revision 1.12  2008/08/14 19:49:00  graceej
% * Better estimate of initial step size.
% * Fixed minor tracing bug.
%
% Revision 1.11  2008/08/14 18:46:37  graceej
% * Fixed bug in truncation of over specified signal.
% * Improved comments.
%
% Revision 1.10  2008/08/14 17:19:53  graceej
% * Describe what the options do and how to set them.
%
% Revision 1.9  2008/08/14 17:02:46  graceej
% * Merged changes to head that delt with refactoring of the code.
%
% Revision 1.8.2.5  2008/08/14 16:59:40  graceej
% * Moved code that determines upscaling to seperate function.
% * Added lengthy commentary justifying the upscaling factor.
% * Removed redundent argument from downscaling function.
%
% Revision 1.8.2.4  2008/08/14 16:28:32  graceej
% * Moved repeated code to a function that attempts to down size the grid if it is not too big.
%
% Revision 1.8.2.3  2008/08/14 15:38:51  graceej
% * Inlined structmerge to manage the merger of options structures.
%
%