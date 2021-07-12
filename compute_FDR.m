function FDR=compute_FDR(gamma,Perfs,pvalues,pi_0)

nbstrats=length(Perfs);

Rplus=sum(Perfs>0 & pvalues<=gamma);

if(Rplus>0)
    Fplus=min(Rplus,.5*nbstrats*pi_0*gamma);
    FDR=Fplus/Rplus;    
else
    FDR=2;    
end


