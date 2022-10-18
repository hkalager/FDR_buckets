% Script developed by Arman Hassanniakalager 
% Rolling forward 1 step ahead (OOS) prediction generator for the
% securities listed in variable "tickerlist"
% Study period defined in "oos_period_range"
% OOS forecast for each model is generated after training each model
% The pool consists of 270 volatility models.
%% Enviroment setting
% Basics
clear;clc;close all;
drpadd=pwd;
if ispc()
    addpath([drpadd,'\univariate']);
    addpath([drpadd,'\multivariate']);
    addpath([drpadd,'\distributions']);
    addpath([drpadd,'\utility']);
    addpath([drpadd,'\timeseries']);
    addpath([drpadd,'\Dataset_ETFS']);
    % Report to Dropbox
    fid=fopen([drpadd,'\VolatilityReport.txt'],'w');
else
    addpath([drpadd,'/univariate']);
    addpath([drpadd,'/multivariate']);
    addpath([drpadd,'/distributions']);
    addpath([drpadd,'/utility']);
    addpath([drpadd,'/timeseries']);
    addpath([drpadd,'/Dataset_ETFS']);
    % Report to Dropbox
    fid=fopen([drpadd,'/VolatilityReport.txt'],'w');
end

txt=['Pool generation started on ',date,' ',num2str(hour(now)),':',num2str(minute(now)),'\n'];
fprintf(fid,txt);
%% Parallel options
delete(gcp('nocreate'))
poolobj=parpool('local',feature('numcores'));
%% Specification of underlying assets
% Series
tickerlist={'SPY','QQQ','GLD','USO'};
oos_period_range=2014:2020;
freq=5; %minutes
%% Loops
for t=1:numel(tickerlist)
    flname=[tickerlist{t},'_M',num2str(freq),'_processed.csv'];
    tbl0=readtable(flname);
    %tbl0{:,'Weekday'}=weekday(tbl0.Date)-1;
    TF1SMP = find(year(tbl0.Date)==oos_period_range(1),1,'first');
    TF2SMP = find(year(tbl0.Date)==oos_period_range(end),1,'last');
    PeriodsPerAnnum=floor((TF2SMP+1-TF1SMP)/numel(oos_period_range));
    num_per=TF2SMP-TF1SMP;
    oos_ser=nan(num_per+1,325);
    parfor s=0:num_per
        tic;
        disp(['Iteration ', num2str(s+1),' ...']);
        startIS=TF1SMP+s-PeriodsPerAnnum;
        finishIS=TF1SMP+s-1;
        tbl1=tbl0(startIS:finishIS,:);
        [predicted,IS_ht]=volpool(tbl1);
        %[predicted,IS_ht]=garchpool(tbl1);
        oos_ser(s+1,:)=predicted;
        %recorded_results{s+1,1}=IS_ht;
        toc;
    end
    oosdate=tbl0.Date(TF1SMP:TF2SMP);
    save([tickerlist{t},'_Pool_M',num2str(freq),...
        '_OOS_',num2str(oos_period_range(1)),'_',...
        num2str(oos_period_range(end))]);
    
    txt=['\nPool created successfully for ',tickerlist{t},' on ', date,' ',num2str(hour(now)),':',num2str(minute(now)),'\n'];
    fprintf(fid,txt);
end
fclose(fid);
record_mdlspec;
delete(poolobj);
