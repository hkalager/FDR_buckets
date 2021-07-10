function [LL,LLS]=skewtloglik(x,mu,sigma2,v,lambda)
% Log likelihood for Hansen's Skewed T distribution
%
% USAGE:
%   [LL,LLS]=skewtloglik(X,MU,SIGMA2,V,LAMBDA)
%
% INPUTS:
%   X      - Standardized T random variables, either scalar or column vector
%   MU     - Mean of X, either scalar or size(x)
%   SIGMA2 - Variance of X, either scalar or size(x)
%   V      - Degree of freedom parameters, either scalar or size(x)
%   LAMBDA - Asymmetry parameter, either scalar or size(x)
%
% OUTPUTS:
%   LL    - Log-likelihood evaluated at X
%   LLS   - Vector of log-likelihoods corresponding to X
%
% COMMENTS:
%   V>2, -1<LAMBDA<1
%
% REFERENCES:
%   [1] Hansen (1994), Intl.Econ.Rev. (35)
%
% See also SKEWTCDF, SKEWTINV, SKEWTRND, SKEWTPDF

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 3    Date: 9/1/2005

[T,K]=size(x);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if K~=1
    error('X must be a column vector');
end

if nargin==5
    if length(mu)~=1 && ~all(size(mu)==[T K])
        error('mu must be either a scalar or the same size as X');
    end
    if any(sigma2<=0)
        sigma2(sigma2<=0)=exp(sigma2(sigma2<=0));
        %error('sigma2 must contain only positive elements')
    end
    if length(sigma2)==1
        sigma2=sigma2*ones(T,K);
    elseif size(sigma2,1)~=T || size(sigma2,2)~=1
        error('sigma2 must be a scalar or a vector with the same dimensions as X');
    end
    if length(v) == 1
        
        if v<=2
            v=max(exp(v),3);
        end
        v = ones(T,1) * v;
    elseif size(v,1)~=T || size(v,2)~=1 || any(v < 2)
        error('V must be greater than 2, and the vector must by T by 1');
    end
    
    if length(lambda) == 1 
        if (lambda<=-1 || lambda>=1)
            lambda=sign(lambda)*.99;
        end
        lambda = ones(T,1) * lambda;
    elseif size(lambda,1)~=T || size(lambda,2)~=1 || any(lambda <= -1) || any(lambda >= 1)
        error('All values in LAMBDA must be between -1 and 1, and the vector must by T by 1');
    end
    x=x-mu;
else
    error('Only 5 inputs supported');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
v=real(v);
%Compute the log likelihood
logc = gammaln((v+1)/2) - gammaln(v/2) - 0.5*log(pi*(v-2));
c = exp(logc);

a = 4*lambda.*c.*((v-2)./(v-1));
logb = 0.5*log(1 + 3*lambda.^2 - a.^2);
b = exp(logb);

stdresid=x./sqrt(sigma2);

find1 = (stdresid<(-a./b));
find2 = (stdresid>=(-a./b));
LL1 = 0;
LL2 = 0;
if ~isempty(find1)
    LL1   = logb(find1) + logc(find1) - (v(find1)+1)./2.*log(1+1./(v(find1)-2).*((b(find1).*stdresid(find1)+a(find1))./(1-lambda(find1))).^2);
end
if ~isempty(find2)
    LL2   = logb(find2) + logc(find2) - (v(find2)+1)./2.*log(1+1./(v(find2)-2).*((b(find2).*stdresid(find2)+a(find2))./(1+lambda(find2))).^2);
end
LL    = sum(LL1) + sum(LL2) - 0.5*sum(log(sigma2));

%Compute the individual log likelihoods if needed
if nargout>1
    LLS=zeros(T,1);
    LLS(find1)=LL1 - 0.5*log(sigma2(find1));
    LLS(find2)=LL2 - 0.5*log(sigma2(find2));
end
