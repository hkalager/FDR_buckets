% Function to estimate the conditional volatility based on HAR model of Corsi (2004)
% Script developed by Arman Hassanniakalager 
% The script is based on some anonymous script found online.
% Created on 23 Jan. 2018
% Last modified 10 Jun. 2019 19:39 BST
%% Inputs:
%   -date_ser: date series to help seperate weeks & months
%   -rv      : realized volatility series
%   -mdl_type: model type HAR-RV (1) or HAR-log(RV) (2)
function [ht,oos_predict,b]=my_HAR_TV(date_ser,rv,mdl_type)
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
    X=[ones(numel(RV_w)-1,1),...
        RV_Day(m_length:end-1),RV_Day(m_length:end-1).*(abs(RV_Day(m_length:end-1)-RV_m(1:end-1)')),...
        RV_w(1:end-1)',RV_w(1:end-1)'.*(abs(RV_w(1:end-1)'-RV_m(1:end-1)')),...
        RV_m(1:end-1)'];
    Y=RV_Day(m_length+1:end);
    b = (X'*X)\(X'*Y);
    Forecasted_Y = X*b;
    oos_predict=[1,RV_Day(end),RV_Day(end)*abs(RV_Day(end)-RV_m(end)),RV_w(end),RV_w(end)*abs(RV_w(end)-RV_m(end)),RV_m(end)]*b;
    ht=[nan(m_length,1);Forecasted_Y];
else
    X=[ones(numel(RV_w)-1,1),...
        log(RV_Day(m_length:end-1)),log(RV_Day(m_length:end-1).*(abs(RV_Day(m_length:end-1)-RV_m(1:end-1)'))),...
        log(RV_w(1:end-1)'),log(RV_w(1:end-1)'.*(abs(RV_w(1:end-1)'-RV_m(1:end-1)'))),...
        log(RV_m(1:end-1)')];
    
    Y=log(RV_Day(m_length+1:end));
    b = (X'*X)\(X'*Y);
    Forecasted_logY = X*b;
    Forecasted_Y=exp(Forecasted_logY);
    logoos_predict=[1,log(RV_Day(end)),log(RV_Day(end)*abs(RV_Day(end)-RV_m(end))),log(RV_w(end)),log(RV_w(end)*abs(RV_w(end)-RV_m(end))),log(RV_m(end))]*b;
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