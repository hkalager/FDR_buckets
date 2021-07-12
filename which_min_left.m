function outarg=which_min_left(x, tol)
if ~exist('tol','var')
    tol=1e-16;
end
[~,index]=min(x);
for j=1:numel(index)
    if index(j)==1
        outarg(j)=1;
    elseif (x(index(j)-1)-x(index(j)))<tol
        outarg(j)=index(j)-1;
    else
        outarg(j)=index(j);
    end
end
end