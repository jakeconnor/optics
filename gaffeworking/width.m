%WIDTH Calculate full width of a complex distribution.
%
%  When applied to a field such as a 2D beam this function returns the
%  width of the second moment of intensity in the given direction
%  corresponding to ISO 11146.
%
%  Example:
%          [x] = makefftspace(8,256);
%          [X,Y]=meshgrid(x);
%          un=exp(-X.^2-Y.^2);
%          fullwidth=width(un, X)
%          % The result should be 2
%
%@bug Incorrectly accounts for even mesh length values. The aliased point
%should be double counted.
%
% See also L2NORM, MESHGRID, MAKEFFTSPACE.

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

function [fullwidth]=width(un, X)
% Calculate the centroid of the beam
L2 = abs(un).^2;
for d = 1:ndims(L2)
    L2 = cumsum(L2,d);
end
L2 = L2(end);
XCEN = X.*abs(un).^2;
for d = 1:ndims(XCEN)
    XCEN = cumsum(XCEN,d);
end
XCEN = XCEN(end)./L2;
% Use l2norm and the ISO 11146 second-moment width
fullwidth=4.*sqrt(l2norm((X-XCEN).*un)./L2);

%% CVS Log
%
% $Log: width.m,v $
% Revision 1.4  2009/10/24 11:08:05  graceej
% * Modified license to make use of the current BSD style open source initiative license.
%
% Revision 1.3  2009/05/06 20:09:46  graceej
% * Updated make package target to include distribution of the license file and README.
%
% Revision 1.2  2009/04/21 12:12:06  graceej
% * Noticed a conceptual bug when dealing with even length meshes.
%
% Revision 1.1  2009/04/17 11:42:32  graceej
% * Brought over auto_checking/BPC-LIB tag=end and checked in.
% * This old repository is now defunct.
%
% Revision 1.2  2009/04/17 11:37:27  graceej
% * Wholescale modification of the entire library.
%
% Revision 1.1.2.1  2008/09/16 14:02:49  graceej
% * Documentation change.
%
% Revision 1.13.2.1  2008/09/09 14:15:29  graceej
% * Modified for N dimensional data.
%
% Revision 1.13  2008/07/14 18:48:31  graceej
% * trivial optimisation
%
% Revision 1.12  2008/06/13 19:29:35  af604
% * Substituted makespace to makefftspace
%
% Revision 1.11  2008/06/10 15:46:34  af604
% * Added makefftspace to see also list
%
% Revision 1.10  2008/06/10 15:37:41  af604
% * Negligible changes
%
% Revision 1.9  2008/06/10 15:36:58  af604
% * Negligible changes
%
% Revision 1.7  2008/06/09 17:44:18  af604
% * Corrected little mistake in the help file.
%
% Revision 1.6  2008/06/09 17:29:30  af604
% * Added author, date, revision record.
%
% Revision 1.5  2008/06/09 17:27:17  af604
% * The function should be to a professional standard now.
%
% Revision 1.4  2008/06/09 17:17:16  af604
% * The function calculates the second moment width correctly
% * The concept of the function has been radically changed with the use of meshgridded inputs.
% 









