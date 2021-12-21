% Created revised 21 Jun 2021, 20:03 BST.
% Last revised 20 Dec 2021.

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
%tickerlistaa={'SPY','QQQ','SHV','LQD','GLD','USO'};
main_ticker={'SPY','QQQ','GLD','USO'};
loss_b_range=[-5,-2,0,1];
Benchmark={'GARCH','GJR-GARCH','HAR'};
freq=5; %mins
%oos_period_range_end=oos_period_range_test+oos_per-1;
try
    Spec_data=load('RV_Pool_270_Spec_tbl');
catch
    record_mdlspec;
    Spec_data=load('RV_Pool_270_Spec_tbl');
end
Bsize=1000;
Bwindow=10;
%% Study period
oos_period_range_test=2014:2020;

%% Liang specification
Max_lambda=0.95;
N_bins=20;
gamma_range=.05:.05:.95;
%% FDR setting
fdrtarget=0.1;
rng(0);
%% Setting for RSW
gamma_rsw=.1;

%% The benchmark indexes
bench_ind=zeros(size(Benchmark));
family_class=Spec_data.Mdl_Class;
for s=1:numel(bench_ind)
    sel_bench=Benchmark{s};
    idx_bench=find(strcmp(family_class,sel_bench));
    bench_ind(s)=idx_bench(1);
    
end
perf_table=table();
iter=0;
for t=1:numel(main_ticker)
    flname=[act_fld,main_ticker{t},'_Pool_M',num2str(freq),'_OOS_2014_2020.mat'];
    load(flname,'oosdate','oos_ser','TF1SMP','TF2SMP','tbl0');
    oosdate_red=oosdate;
    oosdate_red(weekday(oosdate_red)==1)=[];
    disp(['Calculations started for ',main_ticker{t}]);
    
    oos_ser(oos_ser>.01)=.01;
    oos_ser(oos_ser<1e-8)=1e-8;
    oos_ser(isnan(oos_ser))=.01;
    
    oos_ser_tested=oos_ser;
    modelscount=size(oos_ser_tested,2);
    tic;
    poolset_ser=oos_ser_tested;
    target=tbl0{TF1SMP:end,'RVDaily'};
    
    %sigma2=tbl0{TF1SMP+poolsetind(1)-1:TF1SMP+poolsetind(end)-1,4+voltyp};
    indices=stationary_bootstrap((1:size(poolset_ser,1))',Bsize,Bwindow);
    Perf=zeros(1,modelscount);
    Perf_B=zeros(Bsize,modelscount);
    for l=1:numel(loss_b_range)
        iter=iter+1;
        perf_table{iter,'Asset'}=main_ticker(t);
        perf_table{iter,'Year'}={'2014-2020'};
        perf_table{iter,'Robust_b'}=loss_b_range(l);
        [Perf,loss_ser]=robust_loss_fn(poolset_ser,target,loss_b_range(l));
        for b=1:Bsize
            bsdata=poolset_ser(indices(:,b),:);
            bstarget=target(indices(:,b),:);
            Perf_B(b,:)=robust_loss_fn(bsdata,bstarget,loss_b_range(l));
        end
 
        for bi=1:numel(Benchmark)
            
            Bench_Perf=Perf(bench_ind(bi));
            [~,maxind]=max(Perf);
            Bench_Perf_B=Perf_B(:,bench_ind(bi));
            
            % SPA set
            n=size(poolset_ser,1);
            d_bar=Bench_Perf-Perf';
            d_bar_B=(Bench_Perf_B-Perf_B);
            omega_hat=var((n^.5)*d_bar_B)';
            I_val=((n^.5)*d_bar./omega_hat<=-(2*log10(log10(n)))^.5);
            mu_k=I_val.*d_bar;
            reject_set=find(mu_k>0);
            %SPA_PVal=bsds(loss_ser(:,bench_ind(bi)),loss_ser,Bsize,Bwindow,'STANDARD','STATIONARY');
            lbl_column_spa=['SPA_',Benchmark{bi}];
            perf_table{iter,lbl_column_spa}=numel(reject_set);
            
            % RSW test
            reject_set_rsw1=kfwe(Bench_Perf-Perf,(Perf_B-Perf),1,fdrtarget,modelscount);
            
            lbl_column_rsw1=['KStepM1_',Benchmark{bi}];
            perf_table{iter,lbl_column_rsw1}=numel(reject_set_rsw1);
            
            
            reject_set_rsw5=kfwe(Bench_Perf-Perf,(Perf_B-Perf),5,fdrtarget,modelscount);
            
            lbl_column_rsw5=['KStepM5_',Benchmark{bi}];
            perf_table{iter,lbl_column_rsw5}=numel(reject_set_rsw5);
            
            reject_set_rsw10=kfwe(Bench_Perf-Perf,(Perf_B-Perf),10,fdrtarget,modelscount);
            
            lbl_column_rsw10=['KStepM10_',Benchmark{bi}];
            perf_table{iter,lbl_column_rsw10}=numel(reject_set_rsw10);
            
            k_rsw=1;
            reject_set_rsw_fdp=kfwe(Bench_Perf-Perf,(Perf_B-Perf),k_rsw,fdrtarget,modelscount);
            while numel(reject_set_rsw_fdp)>=(k_rsw/gamma_rsw-1)
                k_rsw=k_rsw+1;
                reject_set_rsw_fdp=kfwe(Bench_Perf-Perf,(Perf_B-Perf),k_rsw,fdrtarget,modelscount);
            end
            
            lbl_column_rsw_fdp=['kStepM_FDP_',Benchmark{bi}];
            perf_table{iter,lbl_column_rsw_fdp}=numel(reject_set_rsw_fdp);
            
            disp(['Number of significant models with SPA and benchmark ',...
                Benchmark{bi},' is ', num2str(numel(reject_set))]);
            disp(['Number of significant models with kStepM-FDP and benchmark ',...
                Benchmark{bi},' is ', num2str(numel(reject_set_rsw_fdp))]);

        end
    end
    toc;
    
    
end

fl_lbl=['AllAssets_',num2str(oos_period_range_test(1)),...
        '_',num2str(oos_period_range_test(end)),'_M',num2str(freq),'_OnePiece_FWE.csv'];
    writetable(perf_table,fl_lbl);
