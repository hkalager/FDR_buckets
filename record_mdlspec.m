% The volatility models pool generator based on the following papers:
% Hansen & Lunde (2005), Baillie et al (1996), Corsi (2008) & Chan & Grant (2016)
% Script developed by Arman Hassanniakalager 
% Created on 16 Jan. 2018
% Last modified 21 Jun. 2021 09:07 BST
% The Gaussian distribution is estimated with a t-dist with Infinite dof
mdlpool=struct('ErrDist',[],'MdlClass',[],'Spec',[],'Forecast',[]);
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
            mdlpool.ErrDist{iter}=distname(dist);
            mdlpool.MdlClass{iter}='SV-MA';
            mdlpool.Spec{iter}=['q=',num2str(q)];
        end
    end
end

%% EWMA (18) -- 10 models

for lambda=.8:0.02:.98
    iter=iter+1;
    mdlpool.ErrDist{iter}='No Dist';
    mdlpool.MdlClass{iter}='EWMA';
    mdlpool.Spec{iter}=['lambda=',num2str(lambda)];
end

%% HAR (19) -- 7 Models
% HAR-RV -- 2 models
mdlset={'HAR-RV','HAR-log'};
for mdltype=1:2
    iter=iter+1;    
    mdlpool.ErrDist{iter}='No Dist';
    mdlpool.MdlClass{iter}='HAR';
    mdlpool.Spec{iter}=['Type=',mdlset{mdltype}];
end
% HARQ -- 2 models
mdlset={'HAR-Q','HAR-Q-log'};
for mdltype=1:2
    iter=iter+1;
    mdlpool.ErrDist{iter}='No Dist';
    mdlpool.MdlClass{iter}='HAR';
    mdlpool.Spec{iter}=['Type=',mdlset{mdltype}];
end

% TV-HAR -- 2 models
mdlset={'HAR-TV','HAR-TV-log'};
for mdltype=1:2
    iter=iter+1;
    mdlpool.ErrDist{iter}='No Dist';
    mdlpool.MdlClass{iter}='HAR';
    mdlpool.Spec{iter}=['Type=',mdlset{mdltype}];
end

% HAR-L -- 1 models
mdlset={'HAR-L'};
for mdltype=1
    iter=iter+1;
    mdlpool.ErrDist{iter}='No Dist';
    mdlpool.MdlClass{iter}='HAR';
    mdlpool.Spec{iter}=['Type=',mdlset{mdltype}];
end


%% SMA (20) -- 5 Models

for num_per=[5,10,22,63,126]
    iter=iter+1;
    mdlpool.ErrDist{iter}='No Dist';
    mdlpool.MdlClass{iter}='SMA';
    mdlpool.Spec{iter}=['p=',num2str(num_per)];
end



%% Outputs
Mdl_Dist=mdlpool.ErrDist;
Mdl_Class=mdlpool.MdlClass;
Mdl_Spec=mdlpool.Spec;
spec_tbl=table((1:iter)',Mdl_Dist',Mdl_Class',Mdl_Spec');
spec_tbl.Properties.VariableNames={'ID','Error','Class','Spec'};
writetable(spec_tbl,'ModelSpec_Revised.xlsx');
save(['RV_Pool_',num2str(iter),'_Spec_tbl']);