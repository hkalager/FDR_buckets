% Function to estimate the conditional volatility based on HAR-RV-CJ model of 
% Andersen, Tim Bollerslev, and Francis X. Diebold 2006
% 
% Script developed by Arman Hassanniakalager 
% Created on 18 Oct. 2022
% Last modified 18 Oct. 2022
%% Inputs:
%   -date_ser: date series to help seperate weeks & months
%   -rv      : realized volatility series
%   -J_Comp      : Jump component as introduced in ATF(2006)
%   -mdl_type: model type HAR-RV-CJ (1), HAR-RV-NONLINEAR-CJ (2) or HAR-log(RV)-CJ (3)
function [ht,oos_predict,b]=my_HAR_CJ(date_ser,rv,J_Comp,mdl_type)
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
J_Day=J_Comp;
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
    X=[ones(numel(RV_w)-1,1),RV_Day(m_length:end-1),...
        RV_w(1:end-1)',RV_m(1:end-1)',J_Day(m_length:end-1)];
    Y=RV_Day(m_length+1:end);
    b = (X'*X)\(X'*Y);
    Forecasted_Y = X*b;
    oos_predict=[1,RV_Day(end),RV_w(end),RV_m(end),J_Day(end)]*b;
    ht=[nan(m_length,1);Forecasted_Y];
elseif mdl_type==2
    X=[ones(numel(RV_w)-1,1),RV_Day(m_length:end-1).^(.5),...
        RV_w(1:end-1)'.^(.5),RV_m(1:end-1)'.^(.5),J_Day(m_length:end-1).^(.5)];
    Y=RV_Day(m_length+1:end).^.5;
    b = (X'*X)\(X'*Y);
    Forecasted_sqrtY = X*b;
    Forecasted_Y=(Forecasted_sqrtY).^2;
    sqrt_oos_predict=[1,RV_Day(end)^.5,RV_w(end)^.5,RV_m(end).^5,J_Day(end)^.5]*b;
    oos_predict=sqrt_oos_predict^2;
    ht=[nan(m_length,1);Forecasted_Y];

elseif mdl_type==3
    
    X=[ones(numel(RV_w)-1,1),log(RV_Day(m_length:end-1)),...
        log(RV_w(1:end-1)'),log(RV_m(1:end-1)'),log(1+J_Day(m_length:end-1))];
    Y=log(RV_Day(m_length+1:end));
    b = (X'*X)\(X'*Y);
    Forecasted_logY = X*b;
    Forecasted_Y=exp(Forecasted_logY);
    logoos_predict=[1,log(RV_Day(end)),log(RV_w(end)),log(RV_m(end)),log(1+J_Day(end))]*b;
    oos_predict=exp(logoos_predict);
    ht=[nan(m_length,1);Forecasted_Y];
end


end