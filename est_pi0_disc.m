% Perform the estimation of pi0 according to Liang, 2015
% Script developed by Arman Hassannia Kalager for completion of 3rd chapter,
% Created on 10 Jul 2017 11:00 BST (estimate),
% Last modified 05 Apr 2018 14:04 BST.
% ****************************************
% Input arguments:
%   p_vals:     The p-values calculated
%   n:          Number of bins to consider between 0 and 1.
%   lambda_max: Maximum level of the tuning parameter lambda (optional)
% 
% ****************************************
% Outputs:
%   pi0:        The estimate of true null 
%   lambda:     the associated tuning parameter

function [pi0,lambda]=est_pi0_disc(p_vals, n, lambda_max)
if ~exist('lambda_max','var')
    lambda_max=0.99;
end
[pvalsorted,~]=sort(p_vals);
m=numel(pvalsorted);
target=1/n*(1:(n-1));
candidate=pvalsorted(and(pvalsorted>0,pvalsorted<1));
lambda_vec=unique(candidate(which_min_left(abs(candidate-target),1e-10)));
cum_count=sum(pvalsorted<=lambda_vec');
cum_count(end+1)=numel(pvalsorted);
p_count=cum_count(1);
p_count=[p_count,diff(cum_count)];

if length(lambda_vec)==1 
    pi0=(m-cum_count(any(lambda_vec==pvalsorted)))/(1-lambda_vec)/m; 
    lambda=lambda_vec;
    return;
end

if lambda_max<1
    lambda_vec=lambda_vec(lambda_vec<=lambda_max);
end

gap = [lambda_vec(1), diff(lambda_vec')];
% index_lambda= find(pvalsorted==lambda_vec');
cum_bin_count=cum_count(1:numel(lambda_vec));
% all but the last bin, length B-1
bin_counts=[cum_bin_count(1), diff(cum_bin_count)];
% number of bins
B=length(lambda_vec)+1;
R = cumsum(bin_counts);

tail_m0 = (m-R)./(1-lambda_vec');
temp = bin_counts./gap - tail_m0;
if sum(temp <= 0)>0
    index =find(temp<=0,1,'first');
else
    index = B-1;
    
end
pi0=tail_m0(index)/m; 
pi0=min(pi0,1);
lambda=lambda_vec(index);

end