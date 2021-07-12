% Last revised 21 Jun 2021, 20:03 BST.
clear;
clc;
%act_fld='/Volumes/My Passport/Empirics4/';
%act_fld='C:\Users\Hassannia\Documents\MATLAB\Empirics4\';
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
loss_b_range=[-5,-2,-1,0,1];
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
%% The benchmark indexes
bench_ind=zeros(size(Benchmark));
family_class=Spec_data.Mdl_Class;
for s=1:numel(bench_ind)
    sel_bench=Benchmark{s};
    idx_bench=find(strcmp(family_class,sel_bench));
    bench_ind(s)=idx_bench(1);
    
end

for t=1:numel(main_ticker)
    perf_table=table();
    iter=0;
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
        
        % MCS test first
        [INCLUDEDR] = mcs(loss_ser,fdrtarget,Bsize,Bwindow) ;
        dropped_mcs=size(loss_ser,2)-numel(INCLUDEDR);
        perf_table{iter,'MCS_included'}=numel(INCLUDEDR);
        disp(['Number of dropped models by MCS is ', num2str(size(loss_ser,2)-numel(INCLUDEDR))]);
        for bi=1:numel(Benchmark)
            
            Bench_Perf=Perf(bench_ind(bi));
            [~,maxind]=max(Perf);
            Bench_Perf_B=Perf_B(:,bench_ind(bi));
            pvalues=mypval(Bench_Perf-Perf',(Perf_B-Perf));
            
            [pi_0hat,lambda]=est_pi0_disc(pvalues, N_bins,Max_lambda);
            %pi_0hat=max(pi_0hat,.5);
            opt_gamma=gamma_finder(Bench_Perf-Perf',pvalues,gamma_range,pi_0hat);
            [pi_aplushat, pi_aminushat] = compute_pi_ahat(pvalues, Bench_Perf-Perf', pi_0hat, opt_gamma);
            [PORTFDR, FDRhat] = my_portfolio_FDR_mod(fdrtarget, Bench_Perf-Perf', pvalues, pi_0hat);
            
            % The unlilely case where a benchmark has the best
            % performance
            PORTFDR=(FDRhat~=2).*PORTFDR;
            lbl_column_fdr=['FDR_',Benchmark{bi}];
            perf_table{iter,lbl_column_fdr}=sum(PORTFDR);
            
            % RSW Set
            reject_set_rsw=kfwe(Bench_Perf-Perf,(Perf_B-Perf),5,fdrtarget,20);
            
            lbl_column_rsw=['RSW_',Benchmark{bi}];
            perf_table{iter,lbl_column_rsw}=numel(reject_set_rsw);

            % Show what you got!
            disp(['Number of significant models with RSW and benchmark ',...
                Benchmark{bi},' is ', num2str(numel(reject_set_rsw))]);
            disp(['Number of significant models with FDR and benchmark ',...
                Benchmark{bi},' is ', num2str(sum(PORTFDR))]);
        end
    end
    toc;
    
    lbl=main_ticker{t};
    fl_lbl=[lbl,'_',num2str(oos_period_range_test(1)),...
        '_',num2str(oos_period_range_test(end)),'_OnePiece.csv'];
    writetable(perf_table,fl_lbl);
end

