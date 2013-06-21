%STRUCTMERGE Merge together two structures.
%
%   Produce a structure which consists of the values of a given structure,
%   the fields of which exist in the template structure.  Any fields not
%   present in the given structure are set to be equal to the template
%   structure.
%
%   This is recursively applied if the structure contains structures.
%
%   The merging of structures in this manner is useful for managing options
%   lists.
%
%     MERGED_OPTIONS = STRUCTMERGE(DEFAULTS,OPTIONS)
%
%Example
%
%     template = struct('Name','Donald Duck','Age','21 again','Mass',90,'Comment','');
%     new_data = struct('Comment','Diet time!','Mass',99,'New','Ignored');
%     merged = structmerge(template,new_data)
%
%   Will produce the output
%
%     merged = 
%
%       Name: 'Donald Duck'
%        Age: 21
%       Mass: 99
%    Comment: 'Diet time!'
%
%   Note that the field 'New' in the new data is not included as it does
%   not exist in the template structure.
%
%See also: 

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

function merged=structmerge(template,new)
merged=template;
% Get a list of the field names of a structure.
f_names = fieldnames(template);
% Traverse the new structure.  If the field name exists then check to see
% if that object is itself a structure.  If it is then recurse, if not
% override the template with the new value.
for idx=1:length(f_names)
    field_name = char(f_names(idx));
    % Field exists in new.
    if isfield(new,field_name)
       % If the field is a structure
       if isstruct(template.(field_name))
           % Recursively call this function to merge with the one from
           % the new structure.
           merged.(field_name) = structmerge(template.(field_name),new.(field_name));
       else
           % Simply assign the original field with the old data.
           merged.(field_name) = new.(field_name); 
       end
    end
end

% Now form a list of fields in the New list but not in the default s.
s=setdiff(fieldnames(new),fieldnames(merged));
if ~isempty(s)
    str='';
    for n=1:length(s)
        str=strcat(str,sprintf('%s,',s{n}));
    end
    warning('The fields %s were not present in the default structure.',str);
end


end

%% CVS log
%
% $Log: structmerge.m,v $
% Revision 1.4  2009/10/24 11:08:05  graceej
% * Modified license to make use of the current BSD style open source initiative license.
%
% Revision 1.3  2009/05/06 20:09:45  graceej
% * Updated make package target to include distribution of the license file and README.
%
% Revision 1.2  2009/04/20 13:31:06  graceej
% * Warn if a supplied field is not present in the default structure.
%
% Revision 1.1  2009/04/17 11:42:32  graceej
% * Brought over auto_checking/BPC-LIB tag=end and checked in.
% * This old repository is now defunct.
%
% Revision 1.2  2009/04/17 11:37:27  graceej
% * Wholescale modification of the entire library.
%
% Revision 1.1.2.1  2008/09/16 14:03:50  graceej
% * Add CVS log.
%
%