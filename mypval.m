% Code by Arman Hassanniakalager
% Calculating the efficient p-values based on Bootstap based on the paper
% of Romano and Wolf (2016) paper link:
% http://www.sciencedirect.com/science/article/pii/S0167715216000389

function [pVal] = mypval(Perfs,Perfs_B)
if size(Perfs_B,2)~=numel(Perfs)
    error('Dimensions do not match!');
end
pVal=zeros(numel(Perfs),1);
for i=1:numel(Perfs)
    t=Perfs(i);
    tNull=Perfs_B(:,i);
    pVal(i) = 2*min((length(tNull(tNull >= t))+1)/(length(tNull)+1),(length(tNull(tNull <= t))+1)/(length(tNull))+1);
end


end
