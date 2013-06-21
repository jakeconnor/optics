%DefaultEndIteration Default EndIteration callback.
%
%  This function is designed to be called at the end of every iteration in
%  the main function (evolve).  It can be used as canonical boiler-plate
%  code for devising ones own end of iteration callbacks but should be
%  generally valid for all simluations.
%
%  The callback should contain a 'set', 'begin', 'end' and 'get' 
%  section described below:
%
%    'set'   - When passed 'set' as the first argument user specific
%              initialiseation code is called, this is generally controlled
%              by the user code.
%
%    'begin' - When the argument 'begin' is passed in by 'evolve' 
%              then the function should initialise all its internal 
%              variables.
%
%    'end'   - When the argument 'end' is passed in by 'evolve' then the
%              function should finalise all its internal variables, for
%              example by resizing large buffers to contain just the data
%              of interest.
%
%    'get'   - When passed 'get', usually by an external piece of user code
%              the internal persistent variables should be passed out of
%              the function to the user code.
%
%  By default this function stores a snapshot of the field and its 
%  associated coordinates at (non uniform) intervals along its evolution.
%  The nominal number of snapshots can be controlled by calling the
%  function in the following manner.
%
%    DEFAULTENDITERATION('set',N_SNAP);
%
%See also: GAFFE_DEMO_GAUSSIAN, GAFFE_DEMO_STEEPEN, GAFFE_DEMO_SOLITON,
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
function [IsQuit, z_, x_, u_] =  DefaultEndIteration(varargin)

% This should always be set.
IsQuit = 0;

% Persistent internal variable should be setup, here we are going to have a
% pair of cell arrays, one for the coordinate system, the other for the
% field snapshot.
persistent n z x u N_SNAP Options;

% If we have passed parameters, then see what they are.
if nargin > 0
    switch lower(varargin{1})
        case 'set'
            % Call user specific setup here.  By default we want 'N_SNAP'
            % snapshots, over MAX_Z distance so N_SNAP is setup here.
            N_SNAP =varargin{2};
            Options  = varargin{3};
        case 'begin'
            n = 0;
            z = zeros(1);
            u = cell(N_SNAP,1);
            x = cell(N_SNAP,1);
        case 'end'
            z = z(1:n);
            u = u(1:n);
            x = x(1:n);
        case 'get'
            z_ = z; x_=x; u_=u;
        otherwise
            throw(MException('SanityCheck:EndIteration',sprintf('The end of iteration callback should be called with either ''begin'', ''end'' or ''get'', currently ''%s'' is not supported.',varargin{1})));            
    end
    return
end

% Normal service is calling with no arguments.  Here we simply accrue
% samples.
z_now = evalin('caller','z');
if z_now - z(end) >= Options.MaxZ/N_SNAP
    % Get a snapshot.
    n=n+1;
    z(n) = evalin('caller','z');
    u{n} = evalin('caller','u1');
    x{n} = evalin('caller','x');
    z_now = z(n);
end
end
    
