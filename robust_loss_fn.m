% Loss function as of Patton (2011)
% Inputs:
%   -input: set of alternative models time series
%   -target: target series to compare with input
%   -b: the b parameter in the Patton Robust Class of Loss Functions
% Outputs:
%   -loss: the mean loss from the model
%   -loss_ser: the series of loss from the model

function [loss,loss_ser]=robust_loss_fn(estimate,target,b)

if size(estimate,1)==size(target,2)
    target=target';
end

if size(estimate,1)~=size(target,1)
    error('Dimensions do not match');
end

missinga=find(target==0);
target(missinga,:)=[];
estimate(missinga,:)=[];

if ischar(b)
    b=lower(b);
    switch b
        case 'mse'
            loss_ser=(estimate-target).^2;
        case 'qlike'
            loss_ser=log10(target)+estimate./target;
        otherwise
            error('invalid loss function');
    end
    
elseif isfloat(b)
    b=round(b);
    if or(b>-1, b<-2)
        loss_ser=(((1+b)*(2+b))^(-1))*(estimate.^(b+2)-target.^(b+2))-...
            ((1+b)^-1)*(target.^(b+1).*(estimate-target));
    elseif b==-1
        loss_ser=target-estimate+log10(estimate./target);

    elseif b==-2
        loss_ser=(estimate./target)-log10(estimate./target)-1;

    end
    
end

loss=mean(loss_ser);
end
