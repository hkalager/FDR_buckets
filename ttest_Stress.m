% Import tables as qlike_base for main data
% qlike_high for periods of high stress and 
% qlike_low for periods with low stress

%% QLIKE

[h,p]=ttest(qlike_base{:,4},qlike_high{:,4},'Tail','both','Alpha',.05)
[h,p]=ttest(qlike_base{:,5},qlike_high{:,5},'Tail','both','Alpha',.05)
[h,p]=ttest(qlike_base{:,6},qlike_high{:,6},'Tail','both','Alpha',.05)
[h,p]=ttest(qlike_base{:,7},qlike_high{:,7},'Tail','both','Alpha',.05)
[h,p]=ttest(qlike_base{:,8},qlike_high{:,8},'Tail','both','Alpha',.05)
[h,p]=ttest(qlike_base{:,9},qlike_high{:,9},'Tail','both','Alpha',.05)
[h,p]=ttest(qlike_base{:,10},qlike_high{:,10},'Tail','both','Alpha',.05)
[h,p]=ttest(qlike_base{:,11},qlike_high{:,11},'Tail','both','Alpha',.05)
[h,p]=ttest(qlike_base{:,12},qlike_high{:,12},'Tail','both','Alpha',.05)
[h,p]=ttest(qlike_base{:,13},qlike_high{:,13},'Tail','both','Alpha',.05)
[h,p]=ttest(qlike_base{:,10},qlike_high{:,10},'Tail','both','Alpha',.05)
[h,p]=ttest(qlike_base{:,13},qlike_high{:,13},'Tail','both','Alpha',.05)

[h,p]=ttest(qlike_base{:,4},qlike_low{:,4},'Tail','both','Alpha',.05)
[h,p]=ttest(qlike_base{:,5},qlike_low{:,5},'Tail','both','Alpha',.05)
[h,p]=ttest(qlike_base{:,6},qlike_low{:,6},'Tail','both','Alpha',.05)
[h,p]=ttest(qlike_base{:,7},qlike_low{:,7},'Tail','both','Alpha',.05)
[h,p]=ttest(qlike_base{:,8},qlike_low{:,8},'Tail','both','Alpha',.05)
[h,p]=ttest(qlike_base{:,9},qlike_low{:,9},'Tail','both','Alpha',.05)
[h,p]=ttest(qlike_base{:,10},qlike_low{:,10},'Tail','both','Alpha',.05)
[h,p]=ttest(qlike_base{:,11},qlike_low{:,11},'Tail','both','Alpha',.05)
[h,p]=ttest(qlike_base{:,12},qlike_low{:,12},'Tail','both','Alpha',.05)
[h,p]=ttest(qlike_base{:,13},qlike_low{:,13},'Tail','both','Alpha',.05)
%% MSE
[~,p]=ttest(mse_base{:,4},mse_high{:,4},'Tail','both','Alpha',.05)
[~,p]=ttest(mse_base{:,5},mse_high{:,5},'Tail','both','Alpha',.05)
[~,p]=ttest(mse_base{:,6},mse_high{:,6},'Tail','both','Alpha',.05)
[~,p]=ttest(mse_base{:,7},mse_high{:,7},'Tail','both','Alpha',.05)
[~,p]=ttest(mse_base{:,8},mse_high{:,8},'Tail','both','Alpha',.05)
[~,p]=ttest(mse_base{:,9},mse_high{:,9},'Tail','both','Alpha',.05)
[~,p]=ttest(mse_base{:,10},mse_high{:,10},'Tail','both','Alpha',.05)
[~,p]=ttest(mse_base{:,11},mse_high{:,11},'Tail','both','Alpha',.05)
[~,p]=ttest(mse_base{:,12},mse_high{:,12},'Tail','both','Alpha',.05)
[~,p]=ttest(mse_base{:,13},mse_high{:,13},'Tail','both','Alpha',.05)


[~,p]=ttest(mse_base{:,4},mse_low{:,4},'Tail','both','Alpha',.05)
[~,p]=ttest(mse_base{:,5},mse_low{:,5},'Tail','both','Alpha',.05)
[~,p]=ttest(mse_base{:,6},mse_low{:,6},'Tail','both','Alpha',.05)
[~,p]=ttest(mse_base{:,7},mse_low{:,7},'Tail','both','Alpha',.05)
[~,p]=ttest(mse_base{:,8},mse_low{:,8},'Tail','both','Alpha',.05)
[~,p]=ttest(mse_base{:,9},mse_low{:,9},'Tail','both','Alpha',.05)
[~,p]=ttest(mse_base{:,10},mse_low{:,10},'Tail','both','Alpha',.05)
[~,p]=ttest(mse_base{:,11},mse_low{:,11},'Tail','both','Alpha',.05)
[~,p]=ttest(mse_base{:,12},mse_low{:,12},'Tail','both','Alpha',.05)
[~,p]=ttest(mse_base{:,13},mse_low{:,13},'Tail','both','Alpha',.05)