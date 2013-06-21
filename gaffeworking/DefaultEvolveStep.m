%DefaultEvolveStep Integrate a field u0 forwards a distance dz.
%
%    Given a field U0(X,Y), its Fourier representation FFT_U0(X,Y) and the
%    corresponding frequency space KX,KY this will propagate the field
%    forwards one step of length DZ applying the linear operator given by a
%    function L(KX,KY) in Fourier space and the nonlinear operator given by
%    a function N(X,Y) in real space.   
%
%    [LOCAL_ERR, U1] = DefaultEvolveStep(U0, FFT_U0, DZ, X, KX, @L, @N);
%
%    Initial field U0, step length DZ, coordinate systems and associated
%    Fourier space X,Y; KX,KY.
%
%    [LOCAL_ERROR, U1, FFT_U1, ...
%     UC, FFT_UC, ...
%     UF, FFT_UF, ...
%     FFT_U0] = DefaultEvolveStep(U0, FFT_U0, DZ, X, KX, @L, @N);
%
%    This third-order accurate default integrator is a generalised
%    implementation of the technique outlined in:
%
%    O. V. Sinkin, R. Holzl\"ohner, J. Zweck, and
%    C. R. Menyuk. "Optimization of the split-step Fourier method
%    in modeling optical-fiber communications systems."
%    J. Lightwave Tech., 21(1):61--€“68, 2000.
%
%See also: EVOLVE

% $Author: graceej $ $Date: 2010/05/23 13:17:15 $
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

function [local_error, u1,U1,uc,UC,uf,UF,U0] = DefaultEvolveStep(u0, ...
    U0, ...
    dz, ...
    X, ...
    KX, ...
    OperatorLinear, ...
    OperatorNonlinear)
persistent L2 L4 dz_last

% Only recalculate L2 or L4 (the linear operators for dz/2 and dz/4
% respectively) if the values on which they depend (dz, mesh) differ from
% the previous call.
% @bug Need N dimensional check for difference in L2.
if isempty(dz_last) || dz ~= dz_last || ...
        any(size(u0) ~= size(L2)) || ...
        any(size(u0) ~= size(L4)) 
    L2 = OperatorLinear(0.50*dz,KX,X,u0,U0);
    L4 = OperatorLinear(0.25*dz,KX,X,u0,U0);
    dz_last = dz;
end

% Propagate over two half steps to obtain the fine solution.
UF = U0.*L4; uf = ifftn(UF);
uf = uf.*OperatorNonlinear(0.5*dz,KX,X,uf,UF);
UF = fftn(uf); uf = ifftn(UF.*L2);
uf = uf.*OperatorNonlinear(0.5*dz,KX,X,uf,UF);
UF = fftn(uf).*L4; uf = ifftn(UF);

% Propagate over a single coarse step to obtain the course solution.
UC = U0.*L2; uc = ifftn(UC);
uc = uc.*OperatorNonlinear(dz,KX,X,uc,UC);
UC = fftn(uc).*L2; uc = ifftn(UC);

% Obtain the higher order solution and its FT (which may not be more accurate than uf)
u1 = uf;
U1 = UF;

% Determine the error between the course and fine solutions.
local_error = sqrt(l2norm(uf-uc)./l2norm(uf));
end
%% CVS Log
%
% $Log%
%
