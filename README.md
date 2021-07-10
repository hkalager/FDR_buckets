# MHT_buckets
The codes and files in this repository replicate the results in a manuscript entitled 
"Profiling Volatility Forecasting Models with Multiple Hypothesis Testing"
You can access a free copy of the manuscript through SSRN at:
https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3737477
For any enquiries please get in touch via hassannia@outlook.com

Reproduction:

In order to replicate the results you need to run the m files in MATLAB. The steps are identified in file names as "step0_XXX.m", "step1_XXX.m" et cetra. Each file is accompanied by the necessary guidance within the script.

Data license:

The datasets provided under folder "Dataset_ETFS" are samples from TAQ at WRDS with no column headings. The provided csv files are for evaulation purposes only and do not bear any official information. It is at the viewer's discretion to interpret the columns. We are not in the capacity to publish the data publicly. You can obtain the data from https://wrds-web.wharton.upenn.edu/wrds/

Third-party scripts: 

I) The majority of the codes for the GARCH and SV models are from Kevin Sheppard's MATLAB repository at https://www.kevinsheppard.com/code/matlab/mfe-toolbox/
All script in folders "distributions", "univariate", "multivariate", "utility", and "timeseries" are from Kevin Sheppard along with the scripts "bsds.m" and "stationary_bootstrap.m"; these codes are provided as is subject to change at all time without previous notice.

II) The scripts performing the Romano, Wolf, and Sheikh (2008)'s k-FWE test is from Michael Wolf's repository with the University of Zurich at https://www.econ.uzh.ch/en/people/faculty/wolf/publications.html#Programming_Code
All contents under the folder "kfwe" is from this source ; these codes are provided as is subject to change at all time without previous notice.

III) Parts of the FDR+/- codes are from the accompanying codes for Bajgrowicz and Scaillet (2012).

III) The script "est_pi0_disc.m" is a function estimating the tuning parameter lambda for the estimating pi0 as in Liang (2016)
