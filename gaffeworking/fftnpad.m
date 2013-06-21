%FFTNPAD Pad an fft ordered array.
%   B = FFTNPAD(A) pads an array A such that size(B) == 2*size(A).
%
%   B = FFTNPAD(A,PADSIZE) pads an array A with PADSIZE(k) zeros along the
%   k-th dimension of A.  
%
%   B = FFTNPAD(A,PADSIZE,PADVAL) pads array A with PADVAL (a scalar)
%   instead of with zeros. 
%
%   If the dimensions specified in PADSIZE are smaller than the initial
%   dimensions of the signal it us un-padded.  
%
%See also: FFTSPACE

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

function B=fftnpad(A,varargin)
%% Manage arguments.
if nargin < 2
    PADSIZE=size(A)*2;
else
    PADSIZE=varargin{1};
end
if nargin < 3
    VALUE=0;
else
    VALUE=varargin{2};
end

SZ1=size(A);
ND=numel(SZ1);
%% Determine if padding or truncation is required for each dimension.
pad = zeros(ND,1);
for d=1:ND
    if PADSIZE(d) > SZ1(d)
        pad(d) = 1;
    else
        pad(d) = 0;
    end
end


% Setup output array.
B = ones(PADSIZE)*VALUE;

% Swap source and target depending if we are padding or unpadding.
for d=1:ND
    if ~pad(d)
        tmp=SZ1(d); SZ1(d) = PADSIZE(d); PADSIZE(d) = tmp;
    end
end
%% Determine source and target ranges for each dimension.
for d=1:ND
    InHalf=floor(SZ1(d)/2);
    Offset = mod(SZ1(d),2);
    
    % Determine source and target ranges for M
    src{d} = [1:(InHalf+Offset) (InHalf+Offset+1):SZ1(d)];
    tgt{d} = [1:(InHalf+Offset) (PADSIZE(d)-InHalf+1):PADSIZE(d)];
    
    % If we are unpadding, swap the source and target ranges.
    if ~pad(d)
        tmp=src{d}; src{d}=tgt{d}; tgt{d}=tmp;
    end
   
end
%% Assign source to target using cell {:} notation of indices (see fftshift)
B(tgt{:}) = A(src{:});
%% Log of changes
%
% $Log: fftnpad.m,v $
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
% Revision 1.2  2009/04/17 11:37:25  graceej
% * Wholescale modification of the entire library.
%
% Revision 1.1.2.1  2008/09/16 11:28:29  graceej
% * Name change from fftpadn.
%
% Revision 1.1.2.1  2008/09/05 20:38:57  graceej
% * Initial checkin of fftnpad.  It appears to work as advertised.
%
%