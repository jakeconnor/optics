
%SECOND_ORDER_SOLITON Return a second-order soliton of the nlse.
%
%   Return the field distribution of a second-order soliton of the NLSE, as
%   a function of local time and distance.  Such solitons are often known
%   as 'breathers'.
%
%     U = SECOND_ORDER_SOLITON(TAU,Z);
%
%   See also: GAFFE_DEMO_SOLITON


% $Author: graceej $ $Date: 2010/05/22 20:58:36 $
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

function u = second_order_soliton(tau,zeta)
soliton_numerator = (4.*(cosh(3.*tau) + 3.*exp(4.*i.*zeta).*cosh(tau)).*exp(i.*zeta./2.0));
soliton_denominator = (cosh(4.*tau) + 4.*cosh(2.*tau) + 3.*cos(4.*zeta));
% @bug Fix for limits formally using l'Hopital etc. this is a quick and
% dirty.
u=soliton_numerator./soliton_denominator;
u(isnan(u))=0;
end