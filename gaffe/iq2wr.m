%IQ2WR Convert an inverse Q parameter to W and R.
%   [W,R] = IQ2WR(IQ,K) Converts a given inverse Q parameter with an assumed
%   value of K to the spot size and radius of curvature W and R
%   respectively.
%
%   [W,R] = IQ2WR(IQ) Converts the given inverse Q parameter assuming K=1.
%
%See also WR2IQ, IQ2IQ

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

function [w,r]=iq2wr(iq,varargin)
k=1;
if nargin>1
    k=varargin{1};
end
r=1./real(iq);
w=sqrt(-2./(k.*imag(iq)));

%
% $Log: iq2wr.m,v $
% Revision 1.4  2009/10/24 11:08:05  graceej
% * Modified license to make use of the current BSD style open source initiative license.
%
% Revision 1.3  2009/05/06 20:09:45  graceej
% * Updated make package target to include distribution of the license file and README.
%
% Revision 1.2  2009/04/19 19:03:49  graceej
% * Reworked gaussian ray transfer matrix code.
% * Corrected minor bugs in demos, using the old field for sampling rather than the newly propagated field.
%
% Revision 1.1  2009/04/18 18:24:45  graceej
% * Rationalised q parameter and gaussian related stuff.
%
%