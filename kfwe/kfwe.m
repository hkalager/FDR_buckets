function [toplist, c_K] = kfwe(teststatvec,bootteststatmat,k,alph,Nmax)
%INPUTS:
%teststatvec:    1xS vector of (centered) test statistics, corresponding to the S null hypothesis (called z in the R routine)
%                e.g. if H0: mu=3, muhat = 1, then unstudentized test statistic is (1-3) or studentized (1-3)/stdhat(muhat) 
%bootteststatmat:MxS vector of (centered) bootstrap test statistics (called z.null in the R routine)
%                Note: Here, the observed value muhat (instead of the null-hypothesized value) 
%                needs to be used to center the bootstrap test statistic 
%k:              control of k-Familywise Error Rate (k>=1)
%alph:           alpha as in k-FWE (e.g. 0.1)
%Nmax:           as in operative method of k-StepM (e.g. 20)
%
%OUTPUTS:
%toplist:        indices of null hypothesis that were rejected, according to the columns of teststatvec
%c_K             critical values in each step of the kStepM-routine


S=size(teststatvec,2); M=size(bootteststatmat,1);

%block-size, control of k-FWE, alpha as in GCR, R for number of rejections
Rjmin1 = 0; Rj = Inf;           %number of rejections so that while condition sum(rej)>=k is initially satisfied

        z_T=0;
        %global j theta_0 Rjmin1 Rj X Xbootind theta bl L r t w_Tsort w_Tboot S mu k alph Nmax
        %Computation of test statistic vector and reordering of X according to size of test statistic
		S;
		[w_Tsort,IX] = sort(teststatvec);       %IX returns vector of original indizes ordered in increasing order acc.t. w_T
        w_Tsort = w_Tsort*flipdim(eye(S),1);    %But we need decreasing order or w_T, hence flip w_Tsort and IX by flipped identity matrix
        IX_sort = IX*flipdim(eye(S),1);

		w_Tboot = bootteststatmat(:,IX_sort);

        [c_K,Rj] = critvaluesureduced(w_Tboot,S,w_Tsort,k,alph,Nmax);

toplist = IX_sort(1:Rj);
