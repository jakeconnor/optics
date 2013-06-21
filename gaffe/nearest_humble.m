%NEAREST_HUMBLE Given a number N return the nearest humble numbers.
%
%    Given a number N return the nearest humble numbers that are greater or
%    less than N. 
%
%      [HPLUS,HMINUS]=NEAREST_HUMBLE(N);
%
%    Returns the humble number HPLUS that is equal to or larger than N and
%    HMINUS, the humble number that is equal to or smaller than N.
%
%    If N is a humble number H==H2==N.
%
%See also: HUMBLE

% $Author: graceej $ $Revision: 1.3 $
% $Date: 2009/10/24 11:08:05 $


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

function [varargout]=nearest_humble(n)
global humble_lut;
if isempty(humble_lut)
    humble_lut=load('humbles.mat');
end
if nargout >= 0
    for idx=1:length(n)
        tmp=find(humble_lut.humble_lut >= n(idx));
        if ~isempty(tmp)
            varargout{1}(idx)=humble_lut.humble_lut(tmp(1));
        else
            varargout{1}(idx)=[];
        end
    end
end
if nargout > 1
    for idx=1:length(n)
        tmp=find(humble_lut.humble_lut <= n(idx));
        if ~isempty(tmp)
            varargout{2}(idx)=humble_lut.humble_lut(tmp(end));
        else
            varargout{2}(idx)=[];
        end
    end
end

%% CVS Log
%
% $Log: nearest_humble.m,v $
% Revision 1.3  2009/10/24 11:08:05  graceej
% * Modified license to make use of the current BSD style open source initiative license.
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