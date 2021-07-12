clear;clc;
drpadd=pwd;
if ispc()
    addpath([drpadd,'\Dataset_ETFS']);
    % Report to Dropbox
else
    
    addpath([drpadd,'/Dataset_ETFS']);
    % Report to Dropbox
end

flname='SPY_M5_processed.csv';
tbl0=readtable(flname);
tbl0(1:252,:)=[];
date_Ser=tbl0.Date;
FSI_Vol=tbl0.FSIVol;

figure,
plot(date_Ser,FSI_Vol,'ob','MarkerSize',3,'DatetimeTickFormat','MMM-yy');
hold on;
plot(date_Ser,zeros(size(date_Ser,1)),'-.k','LineWidth',1.1);hold off
xlim([datetime(2013,12,01),datetime(2021,02,01)]);
title('OFR FSI Volatility') 
xtickangle(45)
