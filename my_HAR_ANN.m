% Function to estimate the conditional volatility based on HAR model of Corsi (2004)
% and MLP as the non-linear estimator
% Script developed by Arman Hassanniakalager 
% Created on 18 Oct 2022.
% Last modified 18 Oct 2022.
%% Inputs:
%   -date_ser: date series to help seperate weeks & months
%   -rv      : realized volatility series
%   -mdl_type: model type HAR-RV (1) or HAR-log(RV) (2)
function [ht,oos_predict,b]=my_HAR_ANN(date_ser,rv,mdl_type,layer_size,my_activation)
if nargin<2
    error('Not enough input arguments')
elseif nargin==2
    mdl_type=1;
end
YYYY = year(date_ser);
MM = month(date_ser);
DD = day(date_ser);
Date_NO = datenum(YYYY,MM,DD);
YYYYMM = 1000*YYYY + MM;
RV_Day = rv;
% empirical length of week
w_length=floor(numel(RV_Day)/52);
% empirical length of month
m_length=floor(numel(RV_Day)/12);
iter=0;
for s=m_length:numel(RV_Day)
    iter=iter+1;
    RV_w(iter)=mean(RV_Day(s-w_length+1:s));
    RV_m(iter)=mean(RV_Day(s-m_length+1:s));
end

%% HAR regression
if mdl_type==1
    X=[RV_Day(m_length:end-1),RV_w(1:end-1)',RV_m(1:end-1)'];
    Y=RV_Day(m_length+1:end);
    ann_mdl=fitrnet(X,Y,'LayerSizes',layer_size,'Activations',my_activation);
    
    Forecasted_Y = predict(ann_mdl,X);
    oos_predict=predict(ann_mdl,[RV_Day(end),RV_w(end),RV_m(end)]);
    ht=[nan(m_length,1);Forecasted_Y];
else
    
    X=[log(RV_Day(m_length:end-1)),log(RV_w(1:end-1)'),log(RV_m(1:end-1)')];
    Y=log(RV_Day(m_length+1:end));
    ann_mdl=fitrnet(X,Y,'LayerSizes',layer_size,'Activations',my_activation);
    Forecasted_logY = predict(ann_mdl,X);
    Forecasted_Y=exp(Forecasted_logY);
    logoos_predict=predict(ann_mdl,[log(RV_Day(end)),log(RV_w(end)),log(RV_m(end))]);
    oos_predict=exp(logoos_predict);
    ht=[nan(m_length,1);Forecasted_Y];

end
% % Plot the sequences
% plot(Date_NO(2:end),Y,'b-',Date_NO(2:end),Forecasted_Y,'r-')
% datetick('x','yyyymm','keeplimits')
% title('Observed and forecasted RV based on HAR model: HARRV')
% xlabel('Time')
% ylabel('Realized Volatility')
% legend('Observed RV','Forecasted RV')

end