function [pi_aplushat,pi_aminushat]=compute_pi_ahat(pvalues, Perfs, pi_0hat, gamma)

nbstrats=length(pvalues);

T_plus=sum((pvalues<=gamma)&(Perfs>0))-.5*pi_0hat*nbstrats*gamma;
T_minus=sum((pvalues<=gamma)&(Perfs<0))-.5*pi_0hat*nbstrats*gamma;

T_plus=max(T_plus,0);
T_minus=max(T_minus,0);

pi_aplushat=max(0,T_plus/nbstrats);
pi_aminushat=max(0,T_minus/nbstrats);

