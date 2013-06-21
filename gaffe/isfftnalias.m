%ISFFTNALIAS Check to see if an FFT ordered signal is aliased.
%
%   This function checks to see if an N dimensional FFT ordered signal is
%   aliased along a given dimension by comparing the amount of power
%   outside a given edge to the amount of power.
%
%   If this fraction exceeds a threshold it is considered to be aliased.
%
%      YN = ISFFTNALIAS(U,TOL,WN,DIM);
%
%   The signal is U, and the tolerance is TOL.  The index marking the edge
%   of the region to consider as the guard band, within which the power is
%   measured is WN and the dimension along which to check for aliasing is
%   DIM.
%
%See also: FFTNPAD

% $Author: graceej $ $Date: 2009/10/24 11:08:05 $
% $Revision: 1.4 $


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

function [ISIT] = isfftnalias(U,TOL,WN,DIM)
if nargin < 4 || isempty(DIM)
    DIM=1;
end
if nargin < 3 || isempty(WN)
    WN = floor(0.8*(size(U,DIM)/2));
end
if nargin < 2 || isempty(TOL)
    TOL = 1E-6;
end
% Manage trivial case of WN=0.  If WN=0 simply return true
if WN==0
    ISIT=1;
    return
end
% Shift signal so the dimension to consider is first.
[U,SHIFT] = shiftdim(U,DIM-1);
% Determine total power.
PTotal = abs(U).^2;
for d=1:ndims(PTotal)
    PTotal = cumsum(PTotal,d);
end
% Determine power in guard band.
N=size(U,1);
idx_Guard{1} = [WN:(N-WN+1)];
for d=2:ndims(PTotal)
    idx_Guard{d} = 1:size(U,d);
end
    
PGuard   = abs(U(idx_Guard{:})).^2;
for d=1:ndims(PGuard)
    PGuard = cumsum(PGuard,d);
end
% It is alaised if the power in the guard band exceeds the tolerance.
ISIT = PGuard(end) / PTotal(end) > TOL;

%% Sumall
sumall = @(A) sum(A(1:numel(A)));

%% CVS Log
% 
% $Log: isfftnalias.m,v $
% Revision 1.4  2009/10/24 11:08:05  graceej
% * Modified license to make use of the current BSD style open source initiative license.
%
% Revision 1.3  2009/05/06 20:09:45  graceej
% * Updated make package target to include distribution of the license file and README.
%
% Revision 1.2  2009/05/06 17:57:10  graceej
% * Added CVS Revision information.
%
% Revision 1.1  2009/04/17 11:42:32  graceej
% * Brought over auto_checking/BPC-LIB tag=end and checked in.
% * This old repository is now defunct.
%
% Revision 1.2  2009/04/17 11:37:26  graceej
% * Wholescale modification of the entire library.
%
% Revision 1.1.2.2  2008/11/06 15:34:10  graceej
% * Fix documentation.
%
% Revision 1.1.2.1  2008/09/16 14:08:56  graceej
% * Changed name from isalias.
% * Checked in initial version.
%



