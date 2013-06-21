%FIBONACCI Determine the Nth generalised Fibonacci number
%
%    Returns the Fibonacci number F(n) generalised for real or complex n.
%    This includes the special cases of n being an integer >=0, for which the
%    returned value will be the traditional integer Fibonacci number.
%
%      Fn = FIBONACCI(N);
%
%EXAMPLE
%
%    % Plot the generalised Fibonacci numbers for x=-2:0.1:10 and 
%    % include the bona-fide integer Fibonacci numbers at the integer
%    % positions.
%    X = -5:0.1:11;
%    FX = fibonacci(X);
%    XN = [-5:11];
%    FN = fibonacci(int64(XN));
%    plot(X,FX,'b-',XN,FN,'ro');
%    legend('Generalised F(d)','F_n',2);
%    text(XN,double(FN),num2str(FN'));
%    grid on;
%    title('Integer Fibonacci numbers are red.');
%
%BUGS
% 
%    For large values of the input N expect the results to drift with
%    respect to the actual value of F(N) due to internal floating point
%    errors used in the Binet formula.
%
%See also: IFIBONACCI, ISFIBONACCI

% $Author: graceej $ $Date: 2009/10/24 11:08:05 $
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

function fn=fibonacci(d)
% Golden ratio
phi=(1+sqrt(5))/2;

% If called with an integer input, determine the 
% If we are called with an integer input, make sure we produce an integer
% output.
was_int = find(isinteger(d));
int_type = [];
if was_int
    int_type = class(d);
end
fn=zeros(size(d));
if isinteger(d)
    d=double(d);
else
    % Pseudo-integer mask.
    pmax=find(~(d-double(int64(d))));
end
% Generalised Binet formula.
fn = (phi.^d - cos(d.*pi).*phi.^(-d))./sqrt(5);
% Test if the input was strictly an integer.
if was_int
    fn=cast(round(fn),int_type);
else
    fn(pmax) = floor(fn(pmax));
end
end
%% CVS Log
% 
% $Log: fibonacci.m,v $
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
% Revision 1.1.2.1  2008/09/16 14:00:30  graceej
% * Initial checkin.
%
%




