function [c_K,Rj] = critvaluesureduced(w_Tboot,S,w_Tsort,k,alph,Nmax)
% iteratively calculates the critical value d_K from the bootstrap matrix of studentized test
% statistics called bsstatistics (rows are simulation runs, colums variables). Does iterative testing. S is size(data,2), sigi are estimated
% s.d. in original word, theta_0 are hypothesized values, k, alph, Nmax as in kStepM.

j=0;
%theta_0 = ones(1,S)*-eps;       %vector of hypothesized values w.r.t. test statistic sorted IX!!!
Rjmin1 = 0; Rj = Inf;

%determine set over which to find kmax-values, compute kmax-values
%theta = repmat(w_Tsort,[M 1]);

while Rj>=k && not(Rjmin1==Rj) && not(Rj==S) %convergence criterion implemented as while loop criterion
    j=j+1;
    
    if j==1 %if Rj==Inf, then we are in the first step and Kset={1,...,S} follows
        Kset = 1:S;
        wma = w_Tboot;
        wmat = sort(wma,2);                 %sort elements of each row of matrix w_Tboot,r_s-thetahat_T,r_s in ascending order
        if length(Kset)>=k                  %casewise computation of kmax=|K|-k+1 to avoid |K|-k+1=0
            kmaxaux = length(Kset)-k+1;
        else kmaxaux = 1;
        end
        kmax = wmat(:,kmaxaux);
        
        %invert generalized confidence regions, i.e. reject hypothesis s, 1<=s<=S
        rej = zeros(1,S);
        c_K = quantileR(kmax,1-alph,1);
        for n = 1:S
            if c_K < w_Tsort(n)
                rej(n) = 1;
            end
        end
        Rjmin1 = Rj;
        Rj = sum(rej);
        
    else %if j>1, we need to compute all possible sets K, which is denoted by Kset here
        %Computational shortcut is implemented if number of combinations exceeds tol.
        if nchoosek(Rj,k-1)<=Nmax
            Nstar=Rj;
        else                                        %operative method
            Nstar=k-1;
            while nchoosek(Nstar,k-1)<=Nmax
                Nstar = Nstar+1;
            end
            Nstar=Nstar-1;
        end
        
        Kset = nchoosek((Rj-Nstar+1):Rj,k-1);       %all possible combinations of drawing k-1 indizes without replacement out of {1,...,R_(j-1)}.
        e = size(Kset,1);
        c_Kvec = zeros(1,e);
        
        for d=1:e                                   %compute c_K(1-alpha,k,Phat_T) for each row of Kset
            Ktemp = horzcat(Kset(d,:),(Rj+1):S);    %union of Kset and {R_(j-1)+1,...,S}
            X_K = sort(wma(:,Ktemp),2);             %sort elements of each row (bootstrap run) of matrix X_K in ascending order
            if length(Ktemp)>=k                     %casewise kmax=|K|-k+1 to avoid |K|-k+1=0
                kmaxaux = length(Ktemp)-k+1;
            else kmaxaux = 1;
            end
            kmax = X_K(:,kmaxaux);
            c_Kvec(d) = quantileR(kmax,1-alph,1);   %this is c_K(1-alpha,k,Phat_T) as defined in the paper
        end
        c_K(:,j) = max(c_Kvec);
        c_K(j);
        
        %invert generalized confidence regions, i.e. test hypothesis H_r_s for R_(j-1)+1<=s<=S
        rej = zeros(1,S);                               %vector of w_T-ordered hedge funds. 1 if fund s was significantly rejected.
        for n = (Rj+1):S
            if c_K(j) < w_Tsort(n)
                rej(n) = 1;
            end
        end
        Rjmin1 = Rj;
        Rj = Rj + sum(rej);
    end
end     %end of while loop
