%% This script generates the results in Table 4 in the manuscript
% Created 08 Jul 2021, 09:35 BST.
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
%tickerlistaa={'SPY','QQQ','SHV','LQD','GLD','USO'};
main_ticker={'SPY','QQQ','GLD','USO'};
loss_b_range=[-5,-2,0,1];
Benchmark={'GARCH','GJR-GARCH','HAR'};
freq=5; %mins
%oos_period_range_end=oos_period_range_test+oos_per-1;
try
    Spec_data=load('RV_Pool_325_Spec_tbl.mat');
catch
    record_mdlspec;
    Spec_data=load('RV_Pool_325_Spec_tbl.mat');
end
Bsize=1000;
Bwindow=10;
IS_per=252;
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

    flname=[act_fld,main_ticker{t},'_Pool_M',num2str(freq),...
        '_OOS_2014_2020_',num2str(IS_per),'.mat'];
    load(flname,'oosdate','oos_ser','TF1SMP','TF2SMP','tbl0');
    oosdate_red=oosdate;
    oosdate_red(weekday(oosdate_red)==1)=[];
    disp(['Calculations started for ',main_ticker{t},' ...']);

    oos_ser(oos_ser>.01)=.01;
    oos_ser(oos_ser<1e-8)=1e-8;
    oos_ser(isnan(oos_ser))=.01;

    oos_ser_tested=oos_ser;
    modelscount=size(oos_ser_tested,2);
    oos_idx_range=1:size(oos_ser_tested,1);
    target_idx=oos_idx_range+TF1SMP-1;
        
    poolset_ser=oos_ser_tested(oos_idx_range,:);
    target=tbl0{target_idx,'RVDaily'};

    %sigma2=tbl0{TF1SMP+poolsetind(1)-1:TF1SMP+poolsetind(end)-1,4+voltyp};
    indices=stationary_bootstrap((1:size(poolset_ser,1))',Bsize,Bwindow);
    Perf=zeros(1,modelscount);
    Perf_B=zeros(Bsize,modelscount);
    for l=1:numel(loss_b_range)
        tic;
        iter=iter+1;
        disp(['Case robust loss b = ',num2str(loss_b_range(l)),' ...']);
        perf_table{iter,'Asset'}=main_ticker(t);
        perf_table{iter,'Robust_b'}=loss_b_range(l);
        [Perf,loss_ser]=robust_loss_fn(poolset_ser,target,loss_b_range(l));
        for b=1:Bsize
            bsdata=poolset_ser(indices(:,b),:);
            bstarget=target(indices(:,b),:);
            Perf_B(b,:)=robust_loss_fn(bsdata,bstarget,loss_b_range(l));
        end

        % MCS test first
        [INCLUDEDR] = mcs(loss_ser,fdrtarget,Bsize,Bwindow) ;
        perf_table{iter,'MCS_included'}=numel(INCLUDEDR);
        bucket_mcs=poolset_ser(:,INCLUDEDR);
        bucket_mcs_mse=mean(abs(robust_loss_fn(bucket_mcs,target,0)),'omitnan');
        perf_table{iter,'MCS_bucket_MSE'}=bucket_mcs_mse;

        bucket_mcs_qlike=mean(abs(robust_loss_fn(bucket_mcs,target,-2)),'omitnan');
        perf_table{iter,'MCS_bucket_QLIKE'}=bucket_mcs_qlike;
        disp(['Number of included models by MCS is ', num2str(numel(INCLUDEDR))]);
        for bi=1:numel(Benchmark)
            Bench_Perf=Perf(bench_ind(bi));
            Bench_Ser=poolset_ser(:,bench_ind(bi));
            [~,maxind]=max(Perf);
            Bench_Perf_B=Perf_B(:,bench_ind(bi));
            pvalues=mypval(Bench_Perf-Perf',(Perf_B-Perf));

            try
                [pi_0hat,lambda]=est_pi0_disc(pvalues, N_bins,Max_lambda);
            catch
                pi_0hat=1;
            end
            %pi_0hat=max(pi_0hat,.5);
            opt_gamma=gamma_finder(Bench_Perf-Perf',pvalues,gamma_range,pi_0hat);
            [pi_aplushat, pi_aminushat] = compute_pi_ahat(pvalues, Bench_Perf-Perf', pi_0hat, opt_gamma);
            [PORTFDR, FDRhat] = my_portfolio_FDR_mod(fdrtarget, Bench_Perf-Perf', pvalues, pi_0hat);

            % The unlilely case where a benchmark has the best
            % performance
            PORTFDR=(FDRhat~=2).*PORTFDR;
            lbl_column_fdr=['FDR_',Benchmark{bi}];
            perf_table{iter,lbl_column_fdr}=sum(PORTFDR);
            disp(['Number of significant models with FDR and benchmark ',...
                Benchmark{bi},' is ', num2str(sum(PORTFDR))]);
            if sum(PORTFDR)==0
                PORTFDR(Benchmark{bi})=1;
            end
            bucket_fdr=poolset_ser(:,PORTFDR==1);
            bench_MSE_lbl=[Benchmark{bi},'_MSE'];
            bench_mse_val=abs(robust_loss_fn(Bench_Ser,target,0));
            perf_table{iter,bench_MSE_lbl}=bench_mse_val;

            bench_qlike_lbl=[Benchmark{bi},'_QLIKE'];
            bench_qlike_val=abs(robust_loss_fn(Bench_Ser,target,-2));
            perf_table{iter,bench_qlike_lbl}=bench_qlike_val;

            lbl_column_fdr_bucket_MSE=['FDR_',Benchmark{bi},'_bucket_MSE'];
            bucket_fdr_mse=mean(abs(robust_loss_fn(bucket_fdr,target,0)),'omitnan');
            perf_table{iter,lbl_column_fdr_bucket_MSE}=bucket_fdr_mse;

            lbl_column_fdr_bucket_QLIKE=['FDR_',Benchmark{bi},'_bucket_QLIKE'];
            bucket_fdr_qlike=mean(abs(robust_loss_fn(bucket_fdr,target,-2)),'omitnan');
            perf_table{iter,lbl_column_fdr_bucket_QLIKE}=bucket_fdr_qlike;


            % RSW-FDP Test

            k_rsw=1;
            reject_set_rsw=kfwe(Bench_Perf-Perf,(Perf_B-Perf),k_rsw,fdrtarget,modelscount);
            while numel(reject_set_rsw)>=(k_rsw/gamma_rsw-1)
                k_rsw=k_rsw+1;
                reject_set_rsw=kfwe(Bench_Perf-Perf,(Perf_B-Perf),k_rsw,fdrtarget,modelscount);
            end

            lbl_column_rsw=['KStepM_',Benchmark{bi}];
            perf_table{iter,lbl_column_rsw}=numel(reject_set_rsw);

            disp(['Number of significant models with kStepM and benchmark ',...
                Benchmark{bi},' is ', num2str(numel(reject_set_rsw))]);

            if numel(reject_set_rsw)==0
                reject_set_rsw=bench_ind(bi);
            end

            bucket_rsw=poolset_ser(:,reject_set_rsw);

            lbl_column_rsw_bucket_MSE=['KStepM_',Benchmark{bi},'_bucket_MSE'];
            bucket_rsw_mse=mean(abs(robust_loss_fn(bucket_rsw,target,0)),'omitnan');
            perf_table{iter,lbl_column_rsw_bucket_MSE}=bucket_rsw_mse;

            lbl_column_rsw_bucket_QLIKE=['KStepM_',Benchmark{bi},'_bucket_QLIKE'];
            bucket_rsw_qlike=mean(abs(robust_loss_fn(bucket_rsw,target,-2)),'omitnan');
            perf_table{iter,lbl_column_rsw_bucket_QLIKE}=bucket_rsw_qlike;

        end
        toc;
    end

end
    
fl_lbl=['Quantify_',num2str(oos_period_range_test(1)),...
        '_',num2str(oos_period_range_test(end)),'_M',num2str(freq),...
                '_',num2str(IS_per),'.csv'];
writetable(perf_table,fl_lbl);
