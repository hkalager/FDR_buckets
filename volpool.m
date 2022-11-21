% The volatility models pool generator based on the following papers:
% Hansen & Lunde (2005), Baillie et al (1996), Corsi (2009), Chan & Grant
% (2016), and Bollerslev, Patton, and Quaedvlieg (2016), Andersen, 
% Bollerslev, and Diebold (2007), and Rapach, Strauss, and Zhou (2010)
% Created on 16 Jan. 2018
% Script last revised 2 Nov 2022
% @author: Arman Hassanniakalager GitHub: https://github.com/hkalager
% Common disclaimers apply. Subject to change at all time.
function [predicted,IS_ht,Mdl_Class,Mdl_Spec]=volpool(tbl1)
%% Add necessary paths
warning('off')
%% Read Import
try
    retser=tbl1.Ret;
catch
    retser=tbl1.OCReturn;
    retSqser=tbl1.OCReturnSq;
    RVser=tbl1.RVDaily;
    RQser=tbl1.RQDaily;
    J_Comp=tbl1.JDaily;
end
yser=retser-mean(retser);
dateser=tbl1.Date;

MuCondVar=mean(RVser);
% Distribution types
distname={'NORMAL','STUDENTST','GED','SKEWT'};

% The Gaussian distribution is estimated with a t-dist with Infinite dof
mdlpool=struct('Prediction',[],'Params',[],'ErrDist',[],'MdlClass',[],'Spec',[],'Forecast',[]);
%% Optimization options
%opts=optimoptions(@fmincon,'Display','off','Diagnostics','off');
opts = optimset('fmincon');
opts.Display='off';
opts.Diagnostics = 'off';
%opts.Algorithm = 'sqp';
opts.Algorithm = 'interior-point';
iter=0;
step=1;
%% Creating the pool
%% ARCH Models (1) -- 8 models
for dist=1:numel(distname) % All distributions
    for q=0 % All GARCH Lags
        for p=1:2 % All ARCH lags
            iter=iter+1;
            %mdlARCH = garch('Distribution',struct('Name','t','DoF',DoFF(dist)),'ARCHLags',q);
            try
                    opts  =  optimset(opts , 'MaxFunEvals' , 200*(2+p+q));
                [parameters, ~, ht] = tarch(yser, p, 0, q, distname{dist}, [], [], opts);
                predicted_ht = parameters(1);
                for j=1:p
                    predicted_ht = predicted_ht + parameters(j+1)*(yser(end-j+1)).^2;
                end
                for j=1:q
                    predicted_ht = predicted_ht + parameters(j+p+0+1)*ht(end-j+1) ;
                end
            catch
                ht=NaN;
                parameters=NaN;
                predicted_ht=0;
            end
            %EstMdlARCH = estimate(mdlARCH,yser,'options',opts);
            %predicted=forecast(EstMdlARCH,5,'Y0',yser)';
            mdlpool.Prediction{iter}=ht;
            mdlpool.Forecast(iter,step)=predicted_ht;
            mdlpool.Params{iter}={parameters};
            mdlpool.ErrDist{iter}=distname(dist);
            mdlpool.MdlClass{iter}='ARCH';
            mdlpool.Spec{iter}=['q=',num2str(q),',p=',num2str(p)];
        end
    end
end
    

%% GARCH models (2) -- 16 models

for dist=1:numel(distname) % All distributions
    for q=1:2 % All GARCH Lags
        for p=1:2 % All ARCH lags
            iter=iter+1;
            try
                opts  =  optimset(opts , 'MaxFunEvals' , 200*(2+p+q));
                [parameters, ~, ht] = tarch(yser, p, 0, q, distname{dist}, [], [], opts);
                predicted_ht = parameters(1);
                for j=1:p
                    predicted_ht = predicted_ht + parameters(j+1)*(yser(end-j+1)).^2;
                end
                o=0;
                for j=1:q
                    predicted_ht = predicted_ht + parameters(j+p+0+1)*ht(end-j+1) ;
                end
            catch
                ht=NaN;
                parameters=NaN;
                predicted_ht=0;
            end
            %mdlGARCH = garch('Distribution',struct('Name','t','DoF',DoFF(dist)),'GARCHLags',1:p,'ARCHLags',1:q);
            %EstMdlGARCH = estimate(mdlGARCH,yser,'options',opts);
            %predicted=forecast(EstMdlGARCH,5,'Y0',yser)';
            mdlpool.Prediction{iter}=ht;
            mdlpool.Forecast(iter,step)=predicted_ht;
            mdlpool.Params{iter}={parameters};
            mdlpool.ErrDist{iter}=distname(dist);
            mdlpool.MdlClass{iter}='GARCH';
            mdlpool.Spec{iter}=['q=',num2str(q),',p=',num2str(p)];
        end
    end
end


%% IGARCH models (3) -- 16 models

for dist=1:numel(distname) % All distributions
    for q=1:2 % All GARCH Lags
        for p=1:2 % All ARCH lags
            iter=iter+1;
            try
                opts  =  optimset(opts , 'MaxFunEvals' , 200*(2+p+q));
                [parameters, ~, ht] = igarch(yser, p, q, distname{dist}, [], [],[], opts);

                archParameters = parameters(2:p+q);
                finalParameter = 1 - sum(archParameters);
                predicted_ht = parameters(1);
                for j=1:p
                    predicted_ht = predicted_ht + parameters(j+1)*(yser(end+1-j)).^2;
                end
                for j=1:(q-1)
                    predicted_ht = predicted_ht + parameters(j+p+1)*ht(end+1-j) ;
                end
                % This line is needed since the number of ARCH + GARCH coefficients is
                % 1 less than in a usual GARCH
                predicted_ht = predicted_ht + finalParameter*ht(end+1-q) ;
                %mdlGARCH = garch('Distribution',struct('Name','t','DoF',DoFF(dist)),'GARCHLags',1:p,'ARCHLags',1:q);
                %EstMdlGARCH = estimate(mdlGARCH,yser,'options',opts);
                %predicted=forecast(EstMdlGARCH,5,'Y0',yser)';
            catch
                ht=NaN;
                parameters=NaN;
                predicted_ht=0;
            end               
                
            mdlpool.Prediction{iter}=ht;
            mdlpool.Forecast(iter,step)=predicted_ht;
            mdlpool.Params{iter}={parameters};
            mdlpool.ErrDist{iter}=distname(dist);
            mdlpool.MdlClass{iter}='IGARCH';
            mdlpool.Spec{iter}=['q=',num2str(q),',p=',num2str(p)];
        end
    end
end

%% Taylor/Schwert models (4) -- 16 models

for dist=1:numel(distname) % All distributions
    for q=1:2 % All GARCH Lags
        for p=1:2 % All ARCH lags
            iter=iter+1;
            try
                opts  =  optimset(opts , 'MaxFunEvals' , 200*(2+p+q));
                [parameters, ~, ht] = tarch(yser, p, 0, q, distname{dist}, 1, [], opts);
                sqrtpredicted_ht = parameters(1);
                for j=1:p
                    sqrtpredicted_ht = sqrtpredicted_ht + parameters(j+1)*abs((yser(end-j+1)));
                end
                o=0;
                for j=1:q
                    sqrtpredicted_ht = sqrtpredicted_ht + parameters(j+p+0+1)*ht(end-j+1).^.5 ;
                end
                predicted_ht=sqrtpredicted_ht^2;
            catch
                ht=NaN;
                parameters=NaN;
                predicted_ht=0;
            end
            %mdlGARCH = garch('Distribution',struct('Name','t','DoF',DoFF(dist)),'GARCHLags',1:p,'ARCHLags',1:q);
            %EstMdlGARCH = estimate(mdlGARCH,yser,'options',opts);
            %predicted=forecast(EstMdlGARCH,5,'Y0',yser)';
            mdlpool.Prediction{iter}=ht;
            mdlpool.Forecast(iter,step)=predicted_ht;
            mdlpool.Params{iter}={parameters};
            mdlpool.ErrDist{iter}=distname(dist);
            mdlpool.MdlClass{iter}='Taylor/Schwert';
            mdlpool.Spec{iter}=['q=',num2str(q),',p=',num2str(p)];
        end
    end
end

%% A-GARCH models (5) -- 16 models

for dist=1:numel(distname) % All distributions
    for q=1:2 % All GARCH Lags
        for p=1:2 % All ARCH lags
            iter=iter+1;
            try
                opts  =  optimset(opts , 'MaxFunEvals' , 200*(2+p+q));
                [parameters, ~, ht] = agarch(yser, p, q, 'AGARCH', distname{dist},[], opts);

                gamma = parameters(2+p);
                shock1 = (yser(end-(p+q):end)-gamma).^2;
                shock2 = (yser(end-(p+q):end)-gamma*sqrt(ht(end-(p+q):end))).^2;
                model_type=1;
                if model_type == 1
                    predicted_ht = parameters(1);
                    for j=1:p
                        predicted_ht = predicted_ht + parameters(j+1)*shock1(end+1-j);
                    end
                    for j=1:q
                        predicted_ht = predicted_ht+ parameters(j+2+p)*ht(end+1-j);
                    end
                else
                    predicted_ht = parameters(1);
                    for j=1:p
                        predicted_ht = predicted_ht + parameters(j+1)*shock2(end+1-j);
                    end
                    for j=1:q
                        predicted_ht = predicted_ht+ parameters(j+2+p)*ht(end+1-j);
                    end

                end
            catch
                ht=NaN;
                parameters=NaN;
                predicted_ht=0;
            end
            %mdlGARCH = garch('Distribution',struct('Name','t','DoF',DoFF(dist)),'GARCHLags',1:p,'ARCHLags',1:q);
            %EstMdlGARCH = estimate(mdlGARCH,yser,'options',opts);
            %predicted=forecast(EstMdlGARCH,5,'Y0',yser)';
            mdlpool.Prediction{iter}=ht;
            mdlpool.Forecast(iter,step)=predicted_ht;
            mdlpool.Params{iter}={parameters};
            mdlpool.ErrDist{iter}=distname(dist);
            mdlpool.MdlClass{iter}='A-GARCH';
            mdlpool.Spec{iter}=['q=',num2str(q),',p=',num2str(p)];
        end
    end
end


%% NA-GARCH models (6) -- 16 models

for dist=1:numel(distname) % All distributions
    for q=1:2 % All GARCH Lags
        for p=1:2 % All ARCH lags
            iter=iter+1;
            try
                opts  =  optimset(opts , 'MaxFunEvals' , 200*(2+p+q));
                [parameters, ~, ht] = agarch(yser, p, q, 'NAGARCH', distname{dist},[], opts);

                gamma = parameters(2+p);
                shock1 = (yser(end-(p+q):end)-gamma).^2;
                shock2 = (yser(end-(p+q):end)-gamma*sqrt(ht(end-(p+q):end))).^2;
                model_type=21;
                if model_type == 1
                    predicted_ht = parameters(1);
                    for j=1:p
                        predicted_ht = predicted_ht + parameters(j+1)*shock1(end+1-j);
                    end
                    for j=1:q
                        predicted_ht = predicted_ht+ parameters(j+2+p)*ht(end+1-j);
                    end
                else
                    predicted_ht = parameters(1);
                    for j=1:p
                        predicted_ht = predicted_ht + parameters(j+1)*shock2(end+1-j);
                    end
                    for j=1:q
                        predicted_ht = predicted_ht+ parameters(j+2+p)*ht(end+1-j);
                    end

                end
            catch
                ht=NaN;
                parameters=NaN;
                predicted_ht=0;
            end
            %mdlGARCH = garch('Distribution',struct('Name','t','DoF',DoFF(dist)),'GARCHLags',1:p,'ARCHLags',1:q);
            %EstMdlGARCH = estimate(mdlGARCH,yser,'options',opts);
            %predicted=forecast(EstMdlGARCH,5,'Y0',yser)';
            mdlpool.Prediction{iter}=ht;
            mdlpool.Forecast(iter,step)=predicted_ht;
            mdlpool.Params{iter}={parameters};
            mdlpool.ErrDist{iter}=distname(dist);
            mdlpool.MdlClass{iter}='NA-GARCH';
            mdlpool.Spec{iter}=['q=',num2str(q),',p=',num2str(p)];
        end
    end
end

%% TGARCH (7) -- 16 models
for dist=1:numel(distname) % All distributions
    for q=1:2 % All GARCH Lags
        for p=1:2 % All ARCH lags
            iter=iter+1;
            try
                opts  =  optimset(opts , 'MaxFunEvals' , 200*(2+p+q));
                [parameters, ~, ht] = tarch(yser, p, p, q, distname{dist}, 1, [], opts);
                sqrtpredicted_ht = parameters(1);
                for j=1:p
                    sqrtpredicted_ht = sqrtpredicted_ht + parameters(j+1)*abs(yser(end-j+1));
                end
                o=p;
                for j=1:o
                    sqrtpredicted_ht = sqrtpredicted_ht + parameters(j+p+1)*(abs(yser(end-j+1)).*(yser(end-j+1)<0));
                end
                for j=1:q
                    sqrtpredicted_ht = sqrtpredicted_ht + parameters(j+p+o+1)*(ht(end-j+1).^.5) ;
                end
                predicted_ht=sqrtpredicted_ht^2;
            catch
                ht=NaN;
                parameters=NaN;
                predicted_ht=0;
            end
            %mdlGARCH = garch('Distribution',struct('Name','t','DoF',DoFF(dist)),'GARCHLags',1:p,'ARCHLags',1:q);
            %EstMdlGARCH = estimate(mdlGARCH,yser,'options',opts);
            %predicted=forecast(EstMdlGARCH,5,'Y0',yser)';
            mdlpool.Prediction{iter}=ht;
            mdlpool.Forecast(iter,step)=predicted_ht;
            mdlpool.Params{iter}={parameters};
            mdlpool.ErrDist{iter}=distname(dist);
            mdlpool.MdlClass{iter}='TGARCH';
            mdlpool.Spec{iter}=['q=',num2str(q),',p=',num2str(p)];
        end
    end
end

%% GJR-GARCH (8) -- 16 models

for dist=1:numel(distname) % All distributions
    for q=1:2 % All GARCH Lags
        for p=1:2 % All ARCH lags
            iter=iter+1;
            try
                opts  =  optimset(opts , 'MaxFunEvals' , 200*(2+p+q));
                [parameters, ~, ht] = tarch(yser, p, p, q, distname{dist}, [], [], opts);
                predicted_ht = parameters(1);
                for j=1:p
                    predicted_ht = predicted_ht + parameters(j+1)*(yser(end-j+1)).^2;
                end
                o=p;
                for j=1:o
                    predicted_ht = predicted_ht + parameters(j+p+1)*((yser(end-j+1).^2).*(yser(end-j+1)<0));
                end
                for j=1:q
                    predicted_ht = predicted_ht + parameters(j+p+o+1)*(ht(end-j+1)) ;
                end
            catch
                ht=NaN;
                parameters=NaN;
                predicted_ht=0;
            end
            %mdlGARCH = garch('Distribution',struct('Name','t','DoF',DoFF(dist)),'GARCHLags',1:p,'ARCHLags',1:q);
            %EstMdlGARCH = estimate(mdlGARCH,yser,'options',opts);
            %predicted=forecast(EstMdlGARCH,5,'Y0',yser)';
            mdlpool.Prediction{iter}=ht;
            mdlpool.Forecast(iter,step)=predicted_ht;
            mdlpool.Params{iter}={parameters};
            mdlpool.ErrDist{iter}=distname(dist);
            mdlpool.MdlClass{iter}='GJR-GARCH';
            mdlpool.Spec{iter}=['q=',num2str(q),',p=',num2str(p)];
        end
    end
end


%% log-GARCH (9) -- 16 models

for dist=1:numel(distname) % All distributions
    for q=1:2 % All GARCH Lags
        for p=1:2 % All ARCH lags
            iter=iter+1;
            
            try
                opts  =  optimset(opts , 'MaxFunEvals' , 200*(2+p+q));
                [parameters, ~, ht] = egarchmfe(yser, p, 0, q, distname{dist}, [], opts);
                subconst=sqrt(2/pi);
                data2=yser(end-(p+q):end);
                absdata2=abs(data2)-subconst;
                lnpredicted_ht = parameters(1);
                for j=1:p
                    lnpredicted_ht = lnpredicted_ht + parameters(j+1)*absdata2(end+1-j);
                end
                o=0;
                for j=1:o
                    lnpredicted_ht = lnpredicted_ht + parameters(j+p+1)*data2(end+1-j) ;
                end
                for j=1:q
                    lnpredicted_ht = lnpredicted_ht + parameters(j+p+o+1)*log(ht(end+1-j)) ;
                end
                predicted_ht=exp(lnpredicted_ht);
            catch
                ht=NaN;
                parameters=NaN;
                predicted_ht=0;
            end
            %mdlGARCH = garch('Distribution',struct('Name','t','DoF',DoFF(dist)),'GARCHLags',1:p,'ARCHLags',1:q);
            %EstMdlGARCH = estimate(mdlGARCH,yser,'options',opts);
            %predicted=forecast(EstMdlGARCH,5,'Y0',yser)';
            mdlpool.Prediction{iter}=ht;
            mdlpool.Forecast(iter,step)=predicted_ht;
            mdlpool.Params{iter}={parameters};
            mdlpool.ErrDist{iter}=distname(dist);
            mdlpool.MdlClass{iter}='log-GARCH';
            mdlpool.Spec{iter}=['q=',num2str(q),',p=',num2str(p)];
        end
    end
end

%% E-GARCH (10) -- 16 models

for dist=1:numel(distname) % All distributions
    for q=1:2 % All GARCH Lags
        for p=1:2 % All ARCH lags
            iter=iter+1;
            
            % The following lines try to pass the starting values to EGARCH
            % This is to avoid getting bad starting points
%             start_val_cte=1/3*dist^3-2.5*dist^2+37/6*dist-3;
%             strt_val=zeros(p+p+q+start_val_cte,1);
            try 
                opts  =  optimset(opts , 'MaxFunEvals' , 200*(2+p+p+q));
                [parameters, ~, ht] = egarchmfe(yser, p, p, q, distname{dist},[], opts);
                
                subconst=sqrt(2/pi);
                data2=yser(end-(p+q):end);
                absdata2=abs(data2)-subconst;
                lnpredicted_ht = parameters(1);
                for j=1:p
                    lnpredicted_ht = lnpredicted_ht + parameters(j+1)*absdata2(end+1-j);
                end
                o=p;
                for j=1:o
                    lnpredicted_ht = lnpredicted_ht + parameters(j+p+1)*data2(end+1-j) ;
                end
                for j=1:q
                    lnpredicted_ht = lnpredicted_ht + parameters(j+p+o+1)*log(ht(end+1-j)) ;
                end
                predicted_ht=exp(lnpredicted_ht);
            catch
                ht=NaN;
                parameters=NaN;
                predicted_ht=0;
            end
            %mdlGARCH = garch('Distribution',struct('Name','t','DoF',DoFF(dist)),'GARCHLags',1:p,'ARCHLags',1:q);
            %EstMdlGARCH = estimate(mdlGARCH,yser,'options',opts);
            %predicted=forecast(EstMdlGARCH,5,'Y0',yser)';
            mdlpool.Prediction{iter}=ht;
            mdlpool.Forecast(iter,step)=predicted_ht;
            mdlpool.Params{iter}={parameters};
            mdlpool.ErrDist{iter}=distname(dist);
            mdlpool.MdlClass{iter}='EGARCH';
            mdlpool.Spec{iter}=['q=',num2str(q),',p=',num2str(p)];
        end
    end
end

%% N-GARCH (11) -- 16 models

for dist=1:numel(distname) % All distributions
    for q=1:2 % All GARCH Lags
        for p=1:2 % All ARCH lags
            iter=iter+1;
            try
                opts  =  optimset(opts , 'MaxFunEvals' , 200*(2+p+q));
                [parameters, ~, ht] = aparch(yser, p, 0, q, distname{dist}, [], [], opts);

                o=0;
                delta = parameters(1+p+o+q+1);
                deltainv=2/delta;
                omega = parameters(1);
                alpha = parameters(2:p+1);
                gamma = parameters(p+2:p+o+1);
                beta = parameters(p+o+2:p+o+q+1);
                htdeltaser=ht.^(1/deltainv);
                htdelta = omega;
                for j=1:p
                    if o>=j
                        htdelta = htdelta + alpha(j)*(abs(yser(end+1-j))+gamma(j)*yser(end+1-j))^delta;
                    else
                        htdelta = htdelta + alpha(j)*(abs(yser(end+1-j)))^delta;
                    end
                end
                for j=1:q
                    htdelta = htdelta + beta(j)*htdeltaser(end+1-j);
                end
                htTemp = htdelta^deltainv;
                predicted_ht=htTemp;
            catch
                ht=NaN;
                parameters=NaN;
                predicted_ht=0;
            end
            %mdlGARCH = garch('Distribution',struct('Name','t','DoF',DoFF(dist)),'GARCHLags',1:p,'ARCHLags',1:q);
            %EstMdlGARCH = estimate(mdlGARCH,yser,'options',opts);
            %predicted=forecast(EstMdlGARCH,5,'Y0',yser)';
            mdlpool.Prediction{iter}=ht;
            mdlpool.Forecast(iter,step)=predicted_ht;
            mdlpool.Params{iter}={parameters};
            mdlpool.ErrDist{iter}=distname(dist);
            mdlpool.MdlClass{iter}='NGARCH';
            mdlpool.Spec{iter}=['q=',num2str(q),',p=',num2str(p)];
        end
    end
end


%% A-PARCH (12) -- 16 models
for dist=1:numel(distname) % All distributions
    for q=1:2 % All GARCH Lags
        for p=1:2 % All ARCH lags
            iter=iter+1;
            try
                opts  =  optimset(opts , 'MaxFunEvals' , 200*(2+p+p+q));
                [parameters, ~, ht] = aparch(yser, p, p, q, distname{dist}, [], [], opts);

                o=p;
                delta = parameters(1+p+o+q+1);
                deltainv=2/delta;
                omega = parameters(1);
                alpha = parameters(2:p+1);
                gamma = parameters(p+2:p+o+1);
                beta = parameters(p+o+2:p+o+q+1);
                htdeltaser=ht.^(1/deltainv);
                htdelta = omega;
                for j=1:p
                    if o>=j
                        htdelta = htdelta + alpha(j)*(abs(yser(end+1-j))+gamma(j)*yser(end+1-j))^delta;
                    else
                        htdelta = htdelta + alpha(j)*(abs(yser(end+1-j)))^delta;
                    end
                end
                for j=1:q
                    htdelta = htdelta + beta(j)*htdeltaser(end+1-j);
                end
                htTemp = htdelta^deltainv;
                predicted_ht=htTemp;
            catch
                ht=NaN;
                parameters=NaN;
                predicted_ht=0;
            end
            %mdlGARCH = garch('Distribution',struct('Name','t','DoF',DoFF(dist)),'GARCHLags',1:p,'ARCHLags',1:q);
            %EstMdlGARCH = estimate(mdlGARCH,yser,'options',opts);
            %predicted=forecast(EstMdlGARCH,5,'Y0',yser)';
            mdlpool.Prediction{iter}=ht;
            mdlpool.Forecast(iter,step)=predicted_ht;
            mdlpool.Params{iter}={parameters};
            mdlpool.ErrDist{iter}=distname(dist);
            mdlpool.MdlClass{iter}='APARCH';
            mdlpool.Spec{iter}=['q=',num2str(q),',p=',num2str(p)];
        end
    end
end

%% FI-GARCH (13) -- 16 models

for dist=1:numel(distname) % All distributions
    for q=0:1 % W/WO MA term (beta)
        for p=0:1 % W/WO AR term (phi)
            iter=iter+1;
            try
                opts  =  optimset(opts , 'MaxFunEvals' , 400);
                [parameters, ~, ht] = figarch(yser, p, q, distname{dist}, [], [], opts);

                omega = parameters(1);
                figarchWeightParameters = parameters(2:2+p+q);
                T = size(yser,1);
                truncLag=1000;
                archWeights = figarch_weights(figarchWeightParameters,p,q,truncLag);
                tau = truncLag+1:truncLag+T;
                backCastLength = max(floor(length(yser)^(1/2)),1);
                backCastWeights = .05*(.9.^(0:backCastLength ));
                backCastWeights = backCastWeights/sum(backCastWeights);
                backCast = backCastWeights*((yser(1:backCastLength+1)).^2);
                if backCast==0
                    backCast=cov(epsilon);
                end
                epsilon2 = [zeros(truncLag,1);yser.^2];
                epsilon2(1:truncLag) = backCast;
                predicted_ht=omega + archWeights'*epsilon2(end:-1:end+1-truncLag);
            catch
                ht=NaN;
                parameters=NaN;
                predicted_ht=0;
            end
            mdlpool.Prediction{iter}=ht;
            mdlpool.Forecast(iter,step)=predicted_ht;
            mdlpool.Params{iter}={parameters};
            mdlpool.ErrDist{iter}=distname(dist);
            mdlpool.MdlClass{iter}='FI-GARCH';
            mdlpool.Spec{iter}=['q=',num2str(q),',p=',num2str(p)];
        end
    end
end


%% GARCH-MA (14) -- 16 models

[~, ~, yserMA]=armaxfilter(yser,0,1);
for dist=1:numel(distname) % All distributions
    for q=1:2 % All GARCH Lags
        for p=1:2 % All ARCH lags
            iter=iter+1;
            try
                opts  =  optimset(opts , 'MaxFunEvals' , 200*(2+p+q));
                [parameters, ~, ht] = tarch(yserMA, p, 0, q, distname{dist}, [], [], opts);
                predicted_ht = parameters(1);
                for j=1:p
                    predicted_ht = predicted_ht + parameters(j+1)*(yser(end-j+1)).^2;
                end
                o=0;
                for j=1:o
                    predicted_ht = predicted_ht + parameters(j+p+1)*((yser(end-j+1).^2).*(yser(end-j+1)<0));
                end
                for j=1:q
                    predicted_ht = predicted_ht + parameters(j+p+o+1)*(ht(end-j+1)) ;
                end
            catch
                ht=NaN;
                parameters=NaN;
                predicted_ht=0;
            end
            %mdlGARCH = garch('Distribution',struct('Name','t','DoF',DoFF(dist)),'GARCHLags',1:p,'ARCHLags',1:q);
            %EstMdlGARCH = estimate(mdlGARCH,yser,'options',opts);
            %predicted=forecast(EstMdlGARCH,5,'Y0',yser)';
            mdlpool.Prediction{iter}=ht;
            mdlpool.Forecast(iter,step)=predicted_ht;
            mdlpool.Params{iter}={parameters};
            mdlpool.ErrDist{iter}=distname(dist);
            mdlpool.MdlClass{iter}='GARCH-MA';
            mdlpool.Spec{iter}=['q=',num2str(q),',p=',num2str(p)];
        end
    end
end


%% SV (15) -- 8 models

for dist=1:4 % All distributions
    for q=1:2
        for p=0 % All ARCH lags
            iter=iter+1;
            try
                opts  =  optimset(opts , 'MaxFunEvals' , 200*(2+p+q));
                [parameters, ~, ht] = svegarch(yser, p, 0, q, distname{dist},MuCondVar, [], opts);

                subconst=sqrt(2/pi);
                data2=yser(end-(p+q):end);
                absdata2=abs(data2)-subconst;

                lnpredicted_ht = log(MuCondVar);
                for j=1:p
                    lnpredicted_ht = lnpredicted_ht + parameters(j+1)*absdata2(end+1-j);
                end
                o=0;
                for j=1:o
                    lnpredicted_ht = lnpredicted_ht + parameters(j+p+1)*data2(end+1-j) ;
                end
                for j=1:q
                    lnpredicted_ht = lnpredicted_ht + parameters(j+p+o+1)*(log(ht(end+1-j))-log(MuCondVar)) ;
                end
                predicted_ht=exp(lnpredicted_ht);
            catch
                ht=NaN;
                parameters=NaN;
                predicted_ht=0;
            end
            %mdlGARCH = garch('Distribution',struct('Name','t','DoF',DoFF(dist)),'GARCHLags',1:p,'ARCHLags',1:q);
            %EstMdlGARCH = estimate(mdlGARCH,yser,'options',opts);
            %predicted=forecast(EstMdlGARCH,5,'Y0',yser)';
            mdlpool.Prediction{iter}=ht;
            mdlpool.Forecast(iter,step)=predicted_ht;
            mdlpool.Params{iter}={parameters};
            mdlpool.ErrDist{iter}=distname(dist);
            mdlpool.MdlClass{iter}='SV';
            mdlpool.Spec{iter}=['q=',num2str(q)];
        end
    end
end


%% SV-L (16) -- 16 models

for dist=1:4 % All distributions
    for q=1:2 % All GARCH lags
        for p=1:2 % All ARCH lags
            iter=iter+1;
            try
                opts  =  optimset(opts , 'MaxFunEvals' , 200*(2+p+q));
                [parameters, ~, ht] = svegarch(yser, p, p, q, distname{dist},MuCondVar, [], opts);

                subconst=sqrt(2/pi);
                data2=yser(end-(p+q):end);
                absdata2=abs(data2)-subconst;

                lnpredicted_ht = log(MuCondVar);
                for j=1:p
                    lnpredicted_ht = lnpredicted_ht + parameters(j+1)*absdata2(end+1-j);
                end
                o=p;
                for j=1:o
                    lnpredicted_ht = lnpredicted_ht + parameters(j+p+1)*data2(end+1-j) ;
                end
                for j=1:q
                    lnpredicted_ht = lnpredicted_ht + parameters(j+p+o+1)*(log(ht(end+1-j))-log(MuCondVar)) ;
                end
                predicted_ht=exp(lnpredicted_ht);
            catch
                ht=NaN;
                parameters=NaN;
                predicted_ht=0;
            end
            
            %mdlGARCH = garch('Distribution',struct('Name','t','DoF',DoFF(dist)),'GARCHLags',1:p,'ARCHLags',1:q);
            %EstMdlGARCH = estimate(mdlGARCH,yser,'options',opts);
            %predicted=forecast(EstMdlGARCH,5,'Y0',yser)';
            mdlpool.Prediction{iter}=ht;
            mdlpool.Forecast(iter,step)=predicted_ht;
            mdlpool.Params{iter}={parameters};
            mdlpool.ErrDist{iter}=distname(dist);
            mdlpool.MdlClass{iter}='SV-L';
            mdlpool.Spec{iter}=['q=',num2str(q),',p=',num2str(p)];
        end
    end
end

%% SV-MA (17) -- 8 models

[~, ~, yserMA]=armaxfilter(yser,0,1);
for dist=1:4 % All distributions
    for q=1:2
        for p=0 % All ARCH lags
            iter=iter+1;
            try
                opts  =  optimset(opts , 'MaxFunEvals' , 200*(2+p+q));
                MuCondVar=mean(yser.^2);
                [parameters, ~, ht] = svegarch(yserMA, p, 0, q, distname{dist},MuCondVar, [], opts);

                subconst=sqrt(2/pi);
                data2=yser(end-(p+q):end);
                absdata2=abs(data2)-subconst;

                lnpredicted_ht = log(MuCondVar);
                for j=1:p
                    lnpredicted_ht = lnpredicted_ht + parameters(j+1)*absdata2(end+1-j);
                end
                o=0;
                for j=1:o
                    lnpredicted_ht = lnpredicted_ht + parameters(j+p+1)*data2(end+1-j) ;
                end
                for j=1:q
                    lnpredicted_ht = lnpredicted_ht + parameters(j+p+o+1)*(log(ht(end+1-j))-log(MuCondVar)) ;
                end
                predicted_ht=exp(lnpredicted_ht);
            catch
                ht=NaN;
                parameters=NaN;
                predicted_ht=0;
            end
            %mdlGARCH = garch('Distribution',struct('Name','t','DoF',DoFF(dist)),'GARCHLags',1:p,'ARCHLags',1:q);
            %EstMdlGARCH = estimate(mdlGARCH,yser,'options',opts);
            %predicted=forecast(EstMdlGARCH,5,'Y0',yser)';
            mdlpool.Prediction{iter}=ht;
            mdlpool.Forecast(iter,step)=predicted_ht;
            mdlpool.Params{iter}={parameters};
            mdlpool.ErrDist{iter}=distname(dist);
            mdlpool.MdlClass{iter}='SV-MA';
            mdlpool.Spec{iter}=['q=',num2str(q)];
        end
    end
end

%% EWMA (18) -- 10 models

for lambda=.8:0.02:.98
    iter=iter+1;
    ht = riskmetrics(retser,lambda,[]);
    ht=reshape(ht,[numel(dateser),1]);
    predicted_ht=(1-lambda)*retser(end)^2 + lambda * ht(end);
    mdlpool.Prediction{iter}=ht;
    mdlpool.Forecast(iter,step)=predicted_ht;
    mdlpool.Params{iter}='No Parameters';
    mdlpool.ErrDist{iter}='No Dist';
    mdlpool.MdlClass{iter}='EWMA';
    mdlpool.Spec{iter}=['lambda=',num2str(lambda)];
end

%% HAR (19) -- 10 Models
% HAR-RV -- 2 models

mdlset={'HAR-RV','HAR-log'};
for mdltype=1:2
    iter=iter+1;
    [ht,predicted_ht] = my_HAR(dateser,RVser,mdltype);
    mdlpool.Prediction{iter}=ht;
    mdlpool.Forecast(iter,step)=predicted_ht;
    mdlpool.Params{iter}='No Parameters';
    mdlpool.ErrDist{iter}='No Dist';
    mdlpool.MdlClass{iter}='HAR';
    mdlpool.Spec{iter}=['Type=',mdlset{mdltype}];
end
% HARQ -- 2 models
mdlset={'HAR-Q','HAR-Q-log'};
for mdltype=1:2
    iter=iter+1;
    [ht,predicted_ht] = my_HARQ(dateser,RVser,RQser,mdltype);
    mdlpool.Prediction{iter}=ht;
    mdlpool.Forecast(iter,step)=predicted_ht;
    mdlpool.Params{iter}='No Parameters';
    mdlpool.ErrDist{iter}='No Dist';
    mdlpool.MdlClass{iter}='HAR';
    mdlpool.Spec{iter}=['Type=',mdlset{mdltype}];
end

% HAR-CJ -- 3 models
mdlset={'HAR-CJ','HAR-CJ-sqrt','HAR-CJ-log'};
for mdltype=1:3
    iter=iter+1;
    [ht,predicted_ht] = my_HAR_CJ(dateser,RVser,J_Comp,mdltype);
    mdlpool.Prediction{iter}=ht;
    mdlpool.Forecast(iter,step)=predicted_ht;
    mdlpool.Params{iter}='No Parameters';
    mdlpool.ErrDist{iter}='No Dist';
    mdlpool.MdlClass{iter}='HAR';
    mdlpool.Spec{iter}=['Type=',mdlset{mdltype}];
end

% TV-HAR -- 2 models
mdlset={'HAR-TV','HAR-TV-log'};
for mdltype=1:2
    iter=iter+1;
    [ht,predicted_ht] = my_HAR_TV(dateser,RVser,mdltype);
    mdlpool.Prediction{iter}=ht;
    mdlpool.Forecast(iter,step)=predicted_ht;
    mdlpool.Params{iter}='No Parameters';
    mdlpool.ErrDist{iter}='No Dist';
    mdlpool.MdlClass{iter}='HAR';
    mdlpool.Spec{iter}=['Type=',mdlset{mdltype}];
end

% HAR-L -- 1 models
mdlset={'HAR-L'};
for mdltype=1
    iter=iter+1;
    [ht,predicted_ht] = my_HAR_L(dateser,yser,RVser);
    mdlpool.Prediction{iter}=ht;
    mdlpool.Forecast(iter,step)=predicted_ht;
    mdlpool.Params{iter}='No Parameters';
    mdlpool.ErrDist{iter}='No Dist';
    mdlpool.MdlClass{iter}='HAR';
    mdlpool.Spec{iter}=['Type=',mdlset{mdltype}];
end
%% SVR (20) -- 20 Models
% HAR-SVR -- 20 models
mdlset={'HAR-SVR','HAR-SVR-log'};

for kernel=["rbf","linear"]
    for epsilon_svr=[1e-6,1e-5,1e-4,1e-3,1e-2]
        for mdltype=1:2
            iter=iter+1;
            [ht,predicted_ht] = my_HAR_SVR(dateser,RVser,mdltype,kernel,epsilon_svr);
            mdlpool.Prediction{iter}=ht;
            mdlpool.Forecast(iter,step)=predicted_ht;
            mdlpool.Params{iter}='No Parameters';
            mdlpool.ErrDist{iter}='No Dist';
            mdlpool.MdlClass{iter}='SVR';
            mdlpool.Spec{iter}=['Type=',mdlset{mdltype},';kernel=',kernel,'; epsilon=',num2str(epsilon_svr)];
        end
    end
end
%% ANN (21) -- 30 Models

% HAR-ANN -- 30 models
mdlset={'HAR-ANN','HAR-ANN-log'};

for my_activation=["none","tanh","sigmoid"]
    for layer_size=[1,5,10,50,100]
        for mdltype=1:2
            iter=iter+1;
            [ht,predicted_ht] = my_HAR_ANN(dateser,RVser,mdltype,layer_size,my_activation);
            mdlpool.Prediction{iter}=ht;
            mdlpool.Forecast(iter,step)=predicted_ht;
            mdlpool.Params{iter}='No Parameters';
            mdlpool.ErrDist{iter}='No Dist';
            mdlpool.MdlClass{iter}='ANN';
            mdlpool.Spec{iter}=['Type=',mdlset{mdltype},';neurons=',num2str(layer_size),'; activation=',my_activation];
        end
    end
end


%% SMA (22) -- 5 Models

for num_per=[5,10,22,63,126]
    iter=iter+1;
    ht=zeros(numel(dateser),1);
    for idx=2:numel(dateser)
        num_per_adj=min(idx,num_per);
        ht(idx)=mean(retSqser(idx-num_per_adj+1:idx));
    end
    
    predicted_ht=mean(retSqser(max(end-num_per+1,1):end));
    mdlpool.Prediction{iter}=ht;
    mdlpool.Forecast(iter,step)=predicted_ht;
    mdlpool.Params{iter}='No Parameters';
    mdlpool.ErrDist{iter}='No Dist';
    mdlpool.MdlClass{iter}='SMA';
    mdlpool.Spec{iter}=['p=',num2str(num_per)];
end

%% FC (23) Forecast combination -- 2 Models
predicted=mdlpool.Forecast';
IS_ht=zeros(numel(dateser),length(predicted));

for idx=1:numel(predicted)
    IS_ht(:,idx)=mdlpool.Prediction{1,idx};
end
missing_idx=isnan(IS_ht);
sel_row=(sum(missing_idx,2)==0);
IS_ht=IS_ht(sel_row,:);
target_IS=RVser(sel_row);

for theta=[.9,1]
    iter=iter+1;
    weight_time=ones(numel(target_IS),1);
    for idx=1:numel(target_IS)
        weight_time(idx)=theta^(numel(target_IS)-idx);
    end
    mae=abs(IS_ht-target_IS);
    weighted_mae=weight_time.*mae;
    inv_w=sum(weighted_mae,1);
    weights_fc=(1./inv_w)./sum(1./inv_w);
    predicted_ht=weights_fc*predicted';
    mdlpool.Prediction{iter}=ht;
    mdlpool.Forecast(iter,step)=predicted_ht;
    mdlpool.Params{iter}='No Parameters';
    mdlpool.ErrDist{iter}='No Dist';
    mdlpool.MdlClass{iter}='Comb';
    mdlpool.Spec{iter}=['Theta=',num2str(theta)];

end

%% Outputs

predicted=mdlpool.Forecast';
IS_ht=mdlpool.Prediction';
Mdl_Class=mdlpool.MdlClass;
Mdl_Spec=mdlpool.Spec;
end

