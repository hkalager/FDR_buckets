% Code by Arman Hassanniakalager
% Simple raw calculation of p-value

function [pVal] = mypval(Perfs,Perfs_B)
if size(Perfs_B,2)~=numel(Perfs)
    error('Dimensions do not match!');
end
pVal=zeros(numel(Perfs),1);
for i=1:numel(Perfs)
    t=Perfs(i);
    tNull=Perfs_B(:,i);
    pVal(i) = 2*min((length(tNull(tNull > t)))/(length(tNull)),(length(tNull(tNull < t)))/(length(tNull)));
end


end