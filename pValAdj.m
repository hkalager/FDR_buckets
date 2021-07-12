function [result] = pValAdj(t,tNull)

% PURPOSE: - compute adjusted p-values for resampling-based stepdown multiple testing
% ------------------------------------------------------------------------
% INPUTS: - t = S x 1 vector of test statistics
%         - tNull = M x S matrix of null resampling statistics
% ------------------------------------------------------------------------
% RETURNS: - result = S x 4 matrix that contains
%            - number of hypothesis (column 1)
%            - test statistics (column 2)
%            - raw p-values (column 3)
%            - adjusted p-values (column 4)
% ------------------------------------------------------------------------
% NOTES: the output will be automatically printed out in 'nice form'   
% ------------------------------------------------------------------------

% written by: Michael Wolf
% CREATED  02/16
% UPDATED by Arman

% Copyright (c) 2016, Michael Wolf
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions 
% are met:
% 
% 1. Redistributions of source code must retain the above copyright
% notice, this list of conditions and the following disclaimer.
% 
% 
% 2. Redistributions in binary form must reproduce the above copyright
% notice, this list of conditions and the following disclaimer in the
% documentation and/or other materials provided with the distribution.
% 
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
% "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
% LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
% A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
% HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
% SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
% LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
% DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
% THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
% (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

S = length(t);
dims = size(tNull);
if (dims(2) ~= S)
    error('length of t and col.number of tNull do not match!');
end
dims = size(t);
if (dims(1) == 1)
    t = t';
end
% raw (or undadjusted) p-values
pvRaw = zeros(S,1);
pvAdj = zeros(S,1);
for (i = 1:S)
    pvRaw(i) = pValRaw(t(i),tNull(:,i));
end
% adjusted p-values in the order r_1, ..., r_S
[tr,r] = sort(t,'descend');
tNullr = tNull(:,r);
maxStat = max(tNullr');
pvAdj(1) = pValRaw(tr(1),maxStat);
%disp('Romano and Wolf, 2016 p-value estimation initiated');
for s = 2:S-1
%     if mod(s,ceil((S-1)./20)) == 0
%         disp([num2str(round(100*(s/(S-1)))) '%  completed'])
%     end
    maxStat = max(tNullr(:,s:S)');
    help = pValRaw(tr(s),maxStat);
    pvAdj(s) = max(help,pvAdj(s-1));
end
help = pValRaw(tr(S),tNullr(:,S));
pvAdj(S) = max(help,pvAdj(S-1));
% re-arrange adjusted p-values into orginal order 1, ..., S
rInv(r) = 1:S;
pvAdj = pvAdj(rInv);
% combined output
result = [(1:S)',t,pvRaw,pvAdj];
% vheader = sprintf('H_%d ', 1:S);
% printmat(result(:,2:4), '', vheader,'t-stat pValRaw pValAdj');
end












