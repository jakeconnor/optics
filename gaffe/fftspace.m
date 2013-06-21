%FFTSPACE Create zero first ordered (FFT) coordinate vectors.
%  Given the width of the window (WIDTH) and number of samples (N), this 
%  generates vectors representing the space and k (or omega) space
%  variables. The function also returns the corresponding increments dx 
%  and dkx. 
%
%    [X,DX,KX,DKX] = FFTSPACE(WIDTH,N);
%
%  Alternatively calling it with a previously defined vector, built by
%  FFTSPACE will generate another vector spanning the same range.
%
%    [X,DX,KX,DKX] = FFTSPACE(X,N);
%
%  To reorder the space so that the zero element is in the center of the
%  array use IFFTSHIFT.
% 
% Example: 
%         xwidth=2^5;
%         N=2^8;
%         [x,dx,sx,dsx]=fftspace(xwidth, N);
%         length(x); % Should=2^8
%         length(sx);% Should=2^8
%    
%See also: linspace, logspace, fftshift, ifftshift, FFTNPAD

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

function [x,dx,sx,dsx]=fftspace(X,N)

% Test that N is an integer, throw an exception if it is not.
if abs(floor(N)-ceil(N)) 
  throw(MException('fftspace:noninteger',['Argument N should be' ...
		    ' an integer']));
end
% Test that X is a scalar or a vector, throw an exception if not.
if ~isscalar(X) && ~isvector(X)
    throw(MException('fftspace:nonvector','Argument X should be a scalar or a vector.'));
end


% Catch trivial case of single element
if N==1
    x=0;
    dx=1;
    sx=0;
    dsx=1;
    return
end


% If X is a scalar, giving the mesh width set dx from that. Otherwise set
% the new dx with resepect to the old dx.
if isscalar(X)
    dx=X./N;
else
    dx=(X(2)-X(1)).*length(X)./N; 
end



% Determine spectral sampling interval.
dsx=2*pi/dx/N;

% Return an fftshifted index vector.
nvec = @(M) [0:ceil(M/2)-1 -floor(M/2):1:-1].';

% Determine n vector.
n = nvec(N);

% Generate vectors for coordinates.
x=n*dx;
sx=n*dsx;
    
%% CVS Log
%
% $Log: fftspace.m,v $
% Revision 1.4  2009/10/24 11:08:05  graceej
% * Modified license to make use of the current BSD style open source initiative license.
%
% Revision 1.3  2009/05/06 20:09:45  graceej
% * Updated make package target to include distribution of the license file and README.
%
% Revision 1.2  2009/05/06 17:57:10  graceej
% * Added CVS Revision information.
%
% Revision 1.1  2009/05/06 09:09:55  graceej
% * Renamed makefftspace to fftspace - for consistency with linspace and logspace (and in future fhtspace).
%
% Revision 1.2  2009/04/21 12:50:53  graceej
% * Only allow the use of vectors.
% * Determine old width of space from sample size rather than largest element of space.
%
% Revision 1.1.2.2  2009/04/21 12:47:12  graceej
% * Corrected mesh generation for odd and even sizes with current mesh as argument.
% * Only accepts scalar or vector arguments for the space width.
%
% Revision 1.1.2.1  2009/04/20 17:43:49  graceej
% * First attempt at fixing a potential bug in fftspace.
%
% Revision 1.1  2009/04/17 11:42:32  graceej
% * Brought over auto_checking/BPC-LIB tag=end and checked in.
% * This old repository is now defunct.
%
% Revision 1.2  2009/04/17 11:37:26  graceej
% * Wholescale modification of the entire library.
%
% Revision 1.1.2.5  2008/12/03 17:38:09  graceej
% * Fixed distribution problem + fftspace problem
%
% Revision 1.1.2.4  2008/11/05 13:03:27  graceej
% * Check for N not being an integer.
%
% Revision 1.1.2.3  2008/11/05 13:01:42  graceej
% * Check for N not being an integer.
%
% Revision 1.1.2.2  2008/11/05 12:54:54  graceej
% * Corrected for case of N odd.
%
% Revision 1.1.2.1  2008/09/16 14:34:57  graceej
% * Initial checkin.
%
%
