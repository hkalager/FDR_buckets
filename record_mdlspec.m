% The volatility models pool generator based on the following papers:
% Hansen & Lunde (2005), Baillie et al (1996), Corsi (2008) & Chan & Grant (2016)
% The Gaussian distribution is estimated with a t-dist with Infinite dof
% Created on 16 Jan. 2018
% Script last revised 30 Oct 2022
% @author: Arman Hassanniakalager GitHub: https://github.com/hkalager
% Common disclaimers apply. Subject to change at all time.

mdlpool=struct('ErrDist',[],'Family',[],'MdlClass',[],'Spec',[],'Forecast',[]);
distname={'NORMAL','STUDENTST','GED','SKEWT'};
iter=0;
step=1;
%% Creating the pool
%% ARCH Models (1) -- 8 models
for dist=1:numel(distname) % All distributions
    for q=0 % All GARCH Lags
        for p=1:2 % All ARCH lags
            iter=iter+1;
            %mdlARCH = garch('Distribution',struct('Name','t','DoF',DoFF(dist)),'ARCHLags',q);
            mdlpool.Family{iter}='GARCH';
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
            mdlpool.Family{iter}='GARCH';
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
            mdlpool.Family{iter}='GARCH';
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
            mdlpool.Family{iter}='GARCH';
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
            mdlpool.Family{iter}='GARCH';
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
            mdlpool.Family{iter}='GARCH';
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
            mdlpool.Family{iter}='GARCH';
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
            mdlpool.Family{iter}='GARCH';
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
            mdlpool.Family{iter}='GARCH';
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
            mdlpool.Family{iter}='GARCH';
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
            mdlpool.Family{iter}='GARCH';
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
            mdlpool.Family{iter}='GARCH';
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
            mdlpool.Family{iter}='GARCH';
            mdlpool.ErrDist{iter}=distname(dist);
            mdlpool.MdlClass{iter}='FI-GARCH';
            mdlpool.Spec{iter}=['q=',num2str(q),',p=',num2str(p)];
        end
    end
end


%% GARCH-MA (14) -- 16 models

for dist=1:numel(distname) % All distributions
    for q=1:2 % All GARCH Lags
        for p=1:2 % All ARCH lags
            iter=iter+1;
            mdlpool.Family{iter}='GARCH';
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
            mdlpool.Family{iter}='SV';
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
            mdlpool.Family{iter}='SV';
            mdlpool.ErrDist{iter}=distname(dist);
            mdlpool.MdlClass{iter}='SV-L';
            mdlpool.Spec{iter}=['q=',num2str(q),',p=',num2str(p)];
        end
    end
end

%% SV-MA (17) -- 8 models

for dist=1:4 % All distributions
    for q=1:2
        for p=0 % All ARCH lags
            iter=iter+1;
            mdlpool.Family{iter}='SV';
            mdlpool.ErrDist{iter}=distname(dist);
            mdlpool.MdlClass{iter}='SV-MA';
            mdlpool.Spec{iter}=['q=',num2str(q)];
        end
    end
end

%% EWMA (18) -- 10 models

for lambda=.8:0.02:.98
    iter=iter+1;
    mdlpool.Family{iter}='EWMA';
    mdlpool.ErrDist{iter}='No Dist';
    mdlpool.MdlClass{iter}='EWMA';
    mdlpool.Spec{iter}=['lambda=',num2str(lambda)];
end

%% HAR (19) -- 10 Models
% HAR-RV -- 2 models
mdlset={'HAR-RV','HAR-log'};
for mdltype=1:2
    iter=iter+1;   
    mdlpool.Family{iter}='HAR';
    mdlpool.ErrDist{iter}='No Dist';
    mdlpool.MdlClass{iter}='HAR';
    mdlpool.Spec{iter}=['Type=',mdlset{mdltype}];
end
% HARQ -- 2 models
mdlset={'HAR-Q','HAR-Q-log'};
for mdltype=1:2
    iter=iter+1;
    mdlpool.Family{iter}='HAR';
    mdlpool.ErrDist{iter}='No Dist';
    mdlpool.MdlClass{iter}='HAR';
    mdlpool.Spec{iter}=['Type=',mdlset{mdltype}];
end

% HAR-CJ -- 3 models
mdlset={'HAR-CJ','HAR-CJ-sqrt','HAR-CJ-log'};
for mdltype=1:3
    iter=iter+1;
    mdlpool.Family{iter}='HAR';
    mdlpool.Params{iter}='No Parameters';
    mdlpool.ErrDist{iter}='No Dist';
    mdlpool.MdlClass{iter}='HAR';
    mdlpool.Spec{iter}=['Type=',mdlset{mdltype}];
end

% TV-HAR -- 2 models
mdlset={'HAR-TV','HAR-TV-log'};
for mdltype=1:2
    iter=iter+1;
    mdlpool.Family{iter}='HAR';
    mdlpool.ErrDist{iter}='No Dist';
    mdlpool.MdlClass{iter}='HAR';
    mdlpool.Spec{iter}=['Type=',mdlset{mdltype}];
end

% HAR-L -- 1 models
mdlset={'HAR-L'};
for mdltype=1
    iter=iter+1;
    mdlpool.Family{iter}='HAR';
    mdlpool.ErrDist{iter}='No Dist';
    mdlpool.MdlClass{iter}='HAR';
    mdlpool.Spec{iter}=['Type=',mdlset{mdltype}];
end

%% SVR (20) -- 20 Models
mdlset={'HAR-SVR','HAR-SVR-log'};

for kernel=["rbf","linear"]
    for epsilon_svr=[1e-6,1e-5,1e-4,1e-3,1e-2]
        for mdltype=1:2
            iter=iter+1;
            mdlpool.Params{iter}='No Parameters';
            mdlpool.ErrDist{iter}='No Dist';
            mdlpool.Family{iter}='SVR';
            mdlpool.MdlClass{iter}='SVR';
            mdlpool.Spec{iter}=['Type=',mdlset{mdltype},';kernel=',kernel,'; epsilon=',num2str(epsilon_svr)];
        end
    end
end

%% ANN (21) -- 30 Models
mdlset={'HAR-ANN','HAR-ANN-log'};

for my_activation=["none","tanh","sigmoid"]
    for layer_size=[1,5,10,50,100]
        for mdltype=1:2
            iter=iter+1;
            mdlpool.Family{iter}='ANN';
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
    mdlpool.ErrDist{iter}='No Dist';
    mdlpool.Family{iter}='SMA';
    mdlpool.MdlClass{iter}='SMA';
    mdlpool.Spec{iter}=['p=',num2str(num_per)];
end

%% FC (23) Forecast combination -- 2 Models

for theta=[.9,1]
    iter=iter+1;
    mdlpool.Params{iter}='No Parameters';
    mdlpool.ErrDist{iter}='No Dist';
    mdlpool.MdlClass{iter}='Comb';
    mdlpool.Family{iter}='Comb';
    mdlpool.Spec{iter}=['Theta=',num2str(theta)];

end

%% Outputs
Mdl_Dist=mdlpool.ErrDist;
Mdl_Class=mdlpool.MdlClass;
Mdl_Spec=mdlpool.Spec;
Mdl_Family=mdlpool.Family;
spec_tbl=table((1:iter)',Mdl_Family',Mdl_Dist',Mdl_Class',Mdl_Spec');
spec_tbl.Properties.VariableNames={'ID','Family','Error','Class','Spec'};
writetable(spec_tbl,'ModelSpec_Revised.xlsx');
save(['RV_Pool_',num2str(iter),'_Spec_tbl']);