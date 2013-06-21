%L2NORM Determine the L2 norm of a field.
%
%    The function l2norm determines the L2norm squared of the field UN assuming
%    that the area element deltaX deltaY == 1.
%
%    Example
%       
%        [x,dx] = fftspace(8,256); dy=dx;
%        [X,Y] = meshgrid(x);
%        z = exp(-X.^2-Y.^2);
%        % Should be pi/2
%        halfpi=l2norm(z).*dx.*dy
%
%    See also SUM, MAKEFFTSPACE

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

function [norm]=l2norm(un)
L2 = abs(un).^2;
% Operate over N dimensions.
for d = 1:ndims(un)
    L2 = cumsum(L2,d);
end
norm=L2(end);

%% CVS log
%
% $Log: l2norm.m,v $
% Revision 1.4  2009/10/24 11:08:05  graceej
% * Modified license to make use of the current BSD style open source initiative license.
%
% Revision 1.3  2009/05/06 20:09:45  graceej
% * Updated make package target to include distribution of the license file and README.
%
% Revision 1.2  2009/05/06 17:24:33  graceej
% * Update release date to some standard style.
% * Correct error in gaffe_demo_steepen
%
% Revision 1.1  2009/04/17 11:42:32  graceej
% * Brought over auto_checking/BPC-LIB tag=end and checked in.
% * This old repository is now defunct.
%
% Revision 1.2  2009/04/17 11:37:26  graceej
% * Wholescale modification of the entire library.
%
% Revision 1.1.2.1  2008/09/16 14:05:17  graceej
% * Add CVS log.
%
% Revision 1.6.2.1  2008/09/09 14:15:29  graceej
% * Modified for N dimensional data.
%
% Revision 1.6  2008/06/13 19:29:26  af604
% * Substituted makespace to makefftspace
%
% Revision 1.5  2008/06/10 15:45:39  af604
% * Added makefftspace to see also list
%
% Revision 1.4  2008/06/10 15:33:13  af604
% * Modified example to explain how the function can be used together with makefftspace
%
% Revision 1.3  2008/06/09 15:35:31  af604
% *I am trying to get rid of it...
%
% Revision 1.2  2008/06/09 11:15:18  graceej
% * Simplified calling structure
% * Use sum instead of trapz
% * No longer pass x and y coordinate meshes
%
%