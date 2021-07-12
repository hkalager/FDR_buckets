function [ADFstat_opt,pval_opt]=opt_lag_adf(res_ser)
minaic=1e10;
popt=1000;
for p=0:20
    Mdl=arima(p,0,0);
    [EstMdl,~,logL,~] = estimate(Mdl,res_ser);
    aic = aicbic(logL,p+1);
    if aic<minaic
        minaic=aic;
        popt=p;
        [~,pval_opt,ADFstat_opt]=adftest(res_ser,'lags',popt);
        
    end
end

end

%[~,~,ppstat]=pptest(USO.Ret,'lags',popt);