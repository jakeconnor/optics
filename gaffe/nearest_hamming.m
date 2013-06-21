%NEAREST_HAMMING Given a number N return the nearest hamming numbers.
%
%    Given a number N return the nearest hamming numbers that are greater or
%    less than N. 
%
%      [HPLUS,HMINUS]=NEAREST_HAMMING(N);
%
%    Returns the humble number HPLUS that is equal to or larger than N and
%    HMINUS, the humble number that is equal to or smaller than N.
%
%    If N is a humble number H==H2==N.
%
%See also: HAMMING

% $Author: graceej $ $Revision: 1.2 $
% $Date: 2009/10/26 16:52:03 $

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

function [varargout]=nearest_hamming(n)
global hamming_lut;
if isempty(hamming_lut)
    hamming_lut=load('hammings.mat');
end
if nargout >= 0
    for idx=1:length(n)
        tmp=find(hamming_lut.hamming_lut >= n(idx));
        if ~isempty(tmp)
            varargout{1}(idx)=hamming_lut.hamming_lut(tmp(1));
        else
            varargout{1}(idx)=[];
        end
    end
end
if nargout > 1
    for idx=1:length(n)
        tmp=find(hamming_lut.hamming_lut <= n(idx));
        if ~isempty(tmp)
            varargout{2}(idx)=hamming_lut.hamming_lut(tmp(end));
        else
            varargout{2}(idx)=[];
        end
    end
end

%% CVS Log
%
% $Log: nearest_hamming.m,v $
% Revision 1.2  2009/10/26 16:52:03  graceej
% * Replaced with BSD license.
% * Bumped version number to 0_0_8
%
% Revision 1.1  2009/10/24 11:22:01  graceej
% * Added nearest Hamming number mesh padding lengths.
%
% Revision 1.2  2009/05/06 20:09:45  graceej
% * Updated make package target to include distribution of the license file and README.
%
% Revision 1.1  2009/04/17 11:42:32  graceej
% * Brought over auto_checking/BPC-LIB tag=end and checked in.
% * This old repository is now defunct.
%
% Revision 1.2  2009/04/17 11:37:26  graceej
% * Wholescale modification of the entire library.
%
% Revision 1.1.2.2  2008/11/06 15:32:24  graceej
% * Correctly output number when no output arguments are specified.
%
% Revision 1.1.2.1  2008/09/16 14:34:00  graceej
% * Initial checkin of humble number code.
%
%