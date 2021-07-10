University of Zurich, 22.2.2010, Dan Wunderli

These routines are designed to control the generalized familywise
error rate, k-FWE, as described in the paper "Formalized Data Snooping
Based on Generalized Error Rates".

The main routine is kfwe(). A help routine is critvaluesureduced(). 
quantileR() is used for comparability with the R code k.fwe(). Make sure 
all three routines are in your working directory (type pwd to see it), 
or in a place specified in the matlabpath (can be altered under the menu 
File->Set Path…, or by adding a directory dir to the matlabpath by typing
addpath dir)).

Type help kfwe to see the syntax and appropriate inputs of kfwe(). There 
are no default values, so you need to feed in all inputs. Typical 
default values can be seen by typing help kfwe.

The files ztwo.csv, znulltwo.csv, ztworand.csv, and znulltworand.csv are 
examples to get started quickly. Read them in with 
ztwo = csvread('ztwo.csv') and so on.
For example, you can get started by reading in ztworand and znulltworand.
Then, execute
 rejections = kfwe(znulltworand,ztworand,1,0.1,20)
which gives you the indices of the null hypothesis that can be rejected
while ensuring that the multiple error type I 'k=1'-FWE is controlled.

kfwe() is a generic routine to control the k-FWE based on the Operative
Method of Remark 4.1. The user has to supply a 1 x S vector of test
statistics (computed from the original data), teststatvec, and a 
M x S matrix of `bootstrap' statistics, bootteststatmat. 
It is up to the user to supply:
- basic or studentized statistics
- statistics designed for the one-sided setup or for the two-sided setup
  (in the latter case, absolute values should be used everywhere)
Obviously, how the 'bootstrap' statistics have to be computed depends on
the context (e.g. whether the data are independent or constitute a
time series). Furthermore, one could use a parametric bootstrap, a
nonparametric bootstrap, or even subsampling; various possibilities
are discussed in the paper. However, again, this task is completely up
to the user.

Note that the teststatistic vector and the bootstrap test statitic 
vectors need to be BOTH centered. The teststatistic vector is centered 
with the null-hypothesized value. Each bootstrap test statistic vector 
needs to be centered with the observed value (instead of with the 
null-hypothesized value). The kfwe routine does not center for 
itself (how could it possibly, the null-hypothesized value is not fed 
into kfwe). Look at help kfwe for further explanations.
