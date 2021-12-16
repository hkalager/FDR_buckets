function opt_gamma=gamma_finder(Perfs,pvalues,gamma_range,pi_0)
% Based Barras et al (2010)
opt_gamma=max(gamma_range);
FDR_bias=1;
for gamma=gamma_range
    S_plus_rate=sum((pvalues<=gamma)&(Perfs>0))/numel(pvalues);
    FDR_hat=(pi_0*gamma*.5)/S_plus_rate;
    FDR_bias_temp=abs(FDR_hat-gamma);
    if FDR_bias_temp<FDR_bias
        opt_gamma=gamma;
        FDR_bias=FDR_bias_temp;
    end
end

end