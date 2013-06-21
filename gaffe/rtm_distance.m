%RTM_DISTANCE Generate a ray transfer matrix for distance L
%   RTM_DISTANCE(L) returns a 2x2 ray transfer matrix describing the
%   optical distance L.
%
%See also RTM_LENS

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

function m=rtm_distance(l)
m=[1 l; 0 1];