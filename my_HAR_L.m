% Function to estimate the conditional volatility based on leveraged HAR (LHAR) model of Corsi and Reno (2012)
% In this adaption we drop the jump component for simplicity.
% Script developed by Arman Hassanniakalager 
% The script is based on some anonymous script found online.
% Created on 23 Jan. 2018
% Last modified 10 Jun. 2019 19:03 BST
%% Inputs:
%   -date_ser: date series to help seperate weeks & months
%   -rv      : realized volatility series
function [ht,oos_predict,b]=my_HAR_L(date_ser,yser,rv)
if nargin<3
    error('Not enough input arguments')
end
YYYY = year(date_ser);
MM = month(date_ser);
DD = day(date_ser);
Date_NO = datenum(YYYY,MM,DD);
YYYYMM = 1000*YYYY + MM;
RV_Day = rv;
% empirical length of week
w_length=floor(numel(RV_Day)/52);
iter=0;
for s=w_length:numel(RV_Day)
    iter=iter+1;
    RV_w(iter)=mean(RV_Day(s-w_length+1:s));
    ret_w(iter)=mean(yser(s-w_length+1:s));
end
ret_d=yser(w_length:end);
ret_d_p=ret_d.*(ret_d>=0);
ret_d_n=ret_d.*(ret_d<0);
ret_w_p=ret_w.*(ret_w>=0);
ret_w_n=ret_w.*(ret_w<0);

    
X=[ones(numel(RV_w)-1,1),log(RV_Day(w_length:end-1)),log(RV_w(1:end-1)'),ret_d_p(1:end-1),ret_d_n(1:end-1),ret_w_p(1:end-1)',ret_w_n(1:end-1)'];
Y=log(RV_Day(w_length+1:end));
b = (X'*X)\(X'*Y);
Forecasted_logY = X*b;
Forecasted_Y=exp(Forecasted_logY);
logoos_predict=[1,log(RV_Day(end)),log(RV_w(end)),ret_d_p(end),ret_d_n(end),ret_w_p(end),ret_w_n(end)]*b;
oos_predict=exp(logoos_predict);
ht=[nan(w_length,1);Forecasted_Y];

% % Plot the sequences
% plot(Date_NO(2:end),Y,'b-',Date_NO(2:end),Forecasted_Y,'r-')
% datetick('x','yyyymm','keeplimits')
% title('Observed and forecasted RV based on HAR model: HARRV')
% xlabel('Time')
% ylabel('Realized Volatility')
% legend('Observed RV','Forecasted RV')

end