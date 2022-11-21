function [pVal] = pValRaw(t,tNull)

% PURPOSE: compute raw (or unadjusted) p-value
% ------------------------------------------------------------------------
% INPUTS: - t = univariate test statistic
%         - tNull = vector of null resampling statistics
% ------------------------------------------------------------------------
% RETURNS: - pVal = univariate p-value
% ------------------------------------------------------------------------
% NOTES:   
% ------------------------------------------------------------------------

% written by: Michael Wolf
% CREATED  02/16
% UPDATED 

pVal = (length(tNull(tNull >= t))+1)/(length(tNull)+1);












