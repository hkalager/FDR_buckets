%% Warning: this code clears the workspace and the command window
% codes developed for the manuscript titled 
% "A False Discovery Rate Approach to Optimal Volatility Forecasting Model Selection"
% codes by Arman Hassanniakalager
% This script generates the summary statistics as presented in Table 2 of
% the manuscript.
% The summary stats are stored in the table variable named "perf_table" and
% written to file named "SummaryStats.xlsx"

clear;clc;
drpadd=pwd;
if ispc()
    addpath([drpadd,'\Dataset_ETFS']);
    % Report to Dropbox
else
    
    addpath([drpadd,'/Dataset_ETFS']);
    % Report to Dropbox
end
tickerlist={'SPY','QQQ','GLD','USO'};
perf_table=table();
iter=0;
for t=1:numel(tickerlist)
    iter=iter+1;
    flname=[tickerlist{t},'_M5_processed.csv'];
    tbl0=readtable(flname);
    perf_table{iter,'Asset'}=tickerlist(t);
    ret_ser=tbl0.OCReturn;
    %% Basic stats
    perf_table{iter,'Observations'}=numel(ret_ser);
    perf_table{iter,'MeanX100'}=mean(ret_ser)*100;
    perf_table{iter,'MedianX100'}=median(ret_ser)*100;
    perf_table{iter,'MaxX100'}=max(ret_ser)*100;
    perf_table{iter,'MinX100'}=min(ret_ser)*100;
    perf_table{iter,'StdX100'}=std(ret_ser)*100;
    perf_table{iter,'Skewness'}=skewness(ret_ser);
    perf_table{iter,'Kurtosis'}=kurtosis(ret_ser)-3;
    %% Jarque-Bera test
    [~,p_val_jb,jbstat]=jbtest(ret_ser);
    if p_val_jb<=1e-3
        sig='***';
    elseif p_val_jb<=1e-2
        sig='**';
    elseif p_val_jb<=5e-2
        sig='**';
    else
        sig='';
    end
    perf_table{iter,'JB'}={[num2str(round(jbstat,2)),sig]};
    %% Ljung-Box Q test
    [~,p_val_q,qstat]=lbqtest(ret_ser);
    if p_val_q<=1e-3
        sig='***';
    elseif p_val_q<=1e-2
        sig='**';
    elseif p_val_q<=5e-2
        sig='**';
    else
        sig='';
    end
    perf_table{iter,'Q20'}={[num2str(round(qstat,2)),sig]};
    %% ADF test - based on min-AIC
    [adfstat,pval_adf]=opt_lag_adf(ret_ser);
    
    if pval_adf<=1e-3
        sig='***';
    elseif pval_adf<=1e-2
        sig='**';
    elseif pval_adf<=5e-2
        sig='**';
    else
        sig='';
    end
    perf_table{iter,'ADF'}={[num2str(round(adfstat,2)),sig]};
    %% Phillips-Perron test - based on min-AIC
    [ppstat,pval_pp]=opt_lag_pp(ret_ser);
    
    if pval_pp<=1e-3
        sig='***';
    elseif pval_pp<=1e-2
        sig='**';
    elseif pval_pp<=5e-2
        sig='**';
    else
        sig='';
    end
    perf_table{iter,'PP'}={[num2str(round(ppstat,2)),sig]};
    
end
writetable(perf_table,'SummaryStats.xlsx')