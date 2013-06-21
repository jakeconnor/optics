% $Author: graceej $ $Date: 2009/10/24 11:08:04 $
% $Revision: 1.3 $


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

function N = IsAliasDown2_DoubleFibonacci(u,Options,dim)
edge=fibonacci(ifibonacci(size(u,dim)/2)-2);
N = isfftnalias(u,Options.TolAliased2,edge,dim);
end

%% CVS Log
%
% $Log: IsAliasDown2_DoubleFibonacci.m,v $
% Revision 1.3  2009/10/24 11:08:04  graceej
% * Modified license to make use of the current BSD style open source initiative license.
%
% Revision 1.2  2009/05/06 20:09:44  graceej
% * Updated make package target to include distribution of the license file and README.
%
% Revision 1.1  2009/05/06 16:42:13  graceej
% * Moved from evolvemesh.
%
% Revision 1.1  2009/04/17 11:42:33  graceej
% * Brought over auto_checking/BPC-LIB tag=end and checked in.
% * This old repository is now defunct.
%
% Revision 1.2  2009/04/17 11:37:28  graceej
% * Wholescale modification of the entire library.
%
% Revision 1.1.2.1  2008/09/17 16:35:00  graceej
% * Include functions for signals that are 2*Fn in length.  Ensures even signals.
%
% Revision 1.1.2.3  2008/09/16 14:11:49  graceej
% * isalias changed to isfftnalias.
%
% Revision 1.1.2.2  2008/09/16 12:56:02  graceej
% * Initial checkin.
%
%