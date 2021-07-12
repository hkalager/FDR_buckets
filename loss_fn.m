% Loss function as of Hansen & Lunde (2005)
% The first 6 are as mentioned in the paper
% The 7th loss function is from Mincer–Zarnowitz (MZ) regression R2
% Inputs:
%   -input: set of alternative models time series
%   -target: target series to compare with input
%   -fn_type: evaluation function
% Outputs:
%   -loss: the allocated loss value

function loss=loss_fn(input,target,fn_type)
switch fn_type
    case 'MSE2'
        loss=sum(((input-target)).^2)/size(input,1);
    case 'QLIKE'
        missinga=find(target==0);
        target(missinga,:)=[];
        input(missinga,:)=[];
        loss=sum(log(input)+input.*(target.^-1))/size(input,1);

    otherwise
        error('The loss function is not valid!');
end
