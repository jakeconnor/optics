%EVOLVEMESHSET Choose a set of callbacks for resizing the mesh in EVOLVE.
%
%   Generate an option structure of callbacks for resizing the mesh. The
%   default choice for resizing the mesh of a simulation is to use the
%   nearest hamming number to a Fibonacci number.  This is a good compromise
%   between a natural nested hierachy of mesh sizes (Fibonacci lengths) and
%   the speed of a 7-smooth optimised (split radix) FFT.
%
%   Other options are to use a strict Fibonacci sequence or simple power of
%   two ordering.  Adaption of the mesh can be switched off by choosing the
%   identity option.
%
%     Options.Mesh = EVOLVEMESHSET(MESH_OPTION)
%
%   Options for the mesh are listed below in the section MESH SCHEMES
%
%EXAMPLE
%
%   Using the default options as returned by evolve('defaults') modify the
%   mesh resizing scheme from the default to a dyadic mesh.
%
%     Options = evolve('defaults');
%     Options.Callback.Mesh = evolvemeshset('dyadic');
%
%
%
%MESH SCHEMES
%   
%    DYADIC            Mesh sizes are conventional power of two (2^n). 
%
%    FIBONACCI         Mesh sizes are Fibonacci lengths F(n).
%
%    FIBONACCI&HUMBLE  Mesh sizes are Fibonacci length, padded to the next
%                      highest humble (7-smooth) number.
%
%    FIBONACCI&HAMMING Mesh sizes are Fibonacci length, padded to the next 
%                      highest hamming (5-smooth) number.
%
%    DOUBLEFIBONACCI   Mesh sizes are twice the Fibonacci length F(n).
%
%    LINEAR            Mesh is resized in a linear fashion.
% 
%    IDENTITY          The mesh is not resized.  
%
%See also: EVOLVE

% $Author: graceej $ $Date: 2010/05/24 08:22:45 $
% $Revision: 1.5 $


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

function OptionsMesh = evolvemeshset(varargin)
choice='defaults';
choices = {'dyadic','fibonacci','fibonacci&humble','fibonacci&hamming','identity','doublefibonacci','linear'};
if nargin > 0
    choice = lower(varargin{1});
else
    choices
    return
end

OptionsMesh = struct(...
    'SetInitialSize',[],...
    'IsAliasDown1',[],...
    'IsAliasDown2',[],...
    'GetDownSize',[],...
    'GetUpSize',[]...
    );
switch choice
    case 'dyadic',
        OptionsMesh.SetInitialSize = @SetInitialSize_Dyadic;
        OptionsMesh.IsAliasDown1   = @IsAliasDown1_Dyadic;
        OptionsMesh.IsAliasDown2   = @IsAliasDown2_Dyadic;
        OptionsMesh.GetDownSize    = @GetDownSize_Dyadic;
        OptionsMesh.GetUpSize      = @GetUpSize_Dyadic;
    case 'fibonacci',
        OptionsMesh.SetInitialSize = @SetInitialSize_Fibonacci;
        OptionsMesh.IsAliasDown1   = @IsAliasDown1_Fibonacci;
        OptionsMesh.IsAliasDown2   = @IsAliasDown2_Fibonacci;
        OptionsMesh.GetDownSize    = @GetDownSize_Fibonacci;
        OptionsMesh.GetUpSize      = @GetUpSize_Fibonacci;
    case 'doublefibonacci'
        OptionsMesh.SetInitialSize = @SetInitialSize_DoubleFibonacci;
        OptionsMesh.IsAliasDown1   = @IsAliasDown1_DoubleFibonacci;
        OptionsMesh.IsAliasDown2   = @IsAliasDown2_DoubleFibonacci;
        OptionsMesh.GetDownSize    = @GetDownSize_DoubleFibonacci;
        OptionsMesh.GetUpSize      = @GetUpSize_DoubleFibonacci;
    case 'doublefibonacci&humble'
        OptionsMesh.SetInitialSize = @SetInitialSize_HumbleDoubleFibonacci;
        OptionsMesh.IsAliasDown1   = @IsAliasDown1_HumbleDoubleFibonacci;
        OptionsMesh.IsAliasDown2   = @IsAliasDown2_HumbleDoubleFibonacci;
        OptionsMesh.GetDownSize    = @GetDownSize_HumbleDoubleFibonacci;
        OptionsMesh.GetUpSize      = @GetUpSize_HumbleDoubleFibonacci;
    case {'default','defaults','fibonacci&humble'},
        OptionsMesh.SetInitialSize = @SetInitialSize_HumbleFibonacci;
        OptionsMesh.IsAliasDown1   = @IsAliasDown1_HumbleFibonacci;
        OptionsMesh.IsAliasDown2   = @IsAliasDown2_HumbleFibonacci;
        OptionsMesh.GetDownSize    = @GetDownSize_HumbleFibonacci;
        OptionsMesh.GetUpSize      = @GetUpSize_HumbleFibonacci;
    case {'fibonacci&hamming'},
        OptionsMesh.SetInitialSize = @SetInitialSize_HammingFibonacci;
        OptionsMesh.IsAliasDown1   = @IsAliasDown1_HammingFibonacci;
        OptionsMesh.IsAliasDown2   = @IsAliasDown2_HammingFibonacci;
        OptionsMesh.GetDownSize    = @GetDownSize_HammingFibonacci;
        OptionsMesh.GetUpSize      = @GetUpSize_HammingFibonacci;
    case 'identity',
        OptionsMesh.SetInitialSize = @SetInitialSize_Identity;
        OptionsMesh.IsAliasDown1   = @IsAliasDown1_Identity;
        OptionsMesh.IsAliasDown2   = @IsAliasDown2_Identity;
        OptionsMesh.GetDownSize    = @GetDownSize_Identity;
        OptionsMesh.GetUpSize      = @GetUpSize_Identity;
    case 'linear',
        OptionsMesh.SetInitialSize = @SetInitialSize_Linear;
        OptionsMesh.IsAliasDown1   = @IsAliasDown1_Linear;
        OptionsMesh.IsAliasDown2   = @IsAliasDown2_Linear;
        OptionsMesh.GetDownSize    = @GetDownSize_Linear;
        OptionsMesh.GetUpSize      = @GetUpSize_Linear;
    otherwise
        throw(MException('SanityCheck:MeshResizing',['Option for mesh resizing callback "' choice '" unknown.']));
end
        
        
        