%% This script plots Figures 1 and 2 in the manuscript 
% Script last revised 30 Oct 2022
% @author: Arman Hassanniakalager GitHub: https://github.com/hkalager
% Common disclaimers apply. Subject to change at all time.
clear;
clc;

if ispc
    act_fld=[pwd,'\'];
    addpath([act_fld,'\Dataset_ETFS']);
    addpath([act_fld,'\kfwe']);
elseif ismac
    act_fld=[pwd,'/'];
    addpath([act_fld,'/Dataset_ETFS']);
    addpath([act_fld,'/kfwe']);
end

main_ticker={'SPY','QQQ','GLD','USO'};
tickernames=main_ticker;
%tickernames={'S&P500','NASDAQ','Gold','Light Oil'};
loss_b_range=-2;
switch loss_b_range
    case -2
        loss_label='QLIKE';
    case 0
        loss_label='MSE';
    otherwise
        loss_label=['beta=',num2str(loss_b_range)];
end
Benchmark={'GARCH','GJR-GARCH','HAR'};

try
    Spec_data=load('RV_Pool_325_Spec_tbl.mat');
catch
    record_mdlspec;
    Spec_data=load('RV_Pool_325_Spec_tbl.mat');
end

%% The benchmark indexes
bench_ind=zeros(size(Benchmark));
family_class=Spec_data.Mdl_Class;
for s=1:numel(bench_ind)
    sel_bench=Benchmark{s};
    idx_bench=find(strcmp(family_class,sel_bench));
    bench_ind(s)=idx_bench(1);
    
end
freq=5; %mins
IS_per=252; % the number of days in the training period
i_l=5;
nbin=10;
rng(0);
Perf_bar_matrix=[];
for ta=1:numel(main_ticker)
    flname=[act_fld,main_ticker{ta},'_Pool_M',num2str(freq),...
            '_OOS_2014_2020_',num2str(IS_per),'.mat'];
    load(flname,'oosdate','oos_ser','TF1SMP','TF2SMP','tbl0');
    oos_ser(oos_ser>.01)=.01;
    oos_ser(oos_ser<1e-8)=1e-8;
    oos_ser(isnan(oos_ser))=.01;
    oos_ser_tested=oos_ser;
    poolset_ser=oos_ser_tested;
    target=tbl0{TF1SMP:end,'RVDaily'};
    modelscount=size(oos_ser,2);

    Perf(1,:)=robust_loss_fn(poolset_ser,target,loss_b_range);
        
    Perf=abs(Perf);
    Perf(Perf>5*mean(Perf(:,1)))=5*mean(Perf(:,1));
    Perf_bar=real(mean(Perf,1,'omitnan'));
    Perf_bar_matrix(:,ta)=Perf_bar';
   
end

boxplot(Perf_bar_matrix,'PlotStyle','traditional','Symbol','bo','Colors','k','Widths',.75,'Jitter',.15,'Labels',main_ticker);title(loss_label)

%set(gcf,'WindowState','maximized');