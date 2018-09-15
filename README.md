# Estimation of extreme depth‐based quantile regions: dataset and computer codes
Dataset and computer codes used in He and Einmahl (2017), `Estimation of extreme depth‐based quantile regions', JRSS-b

Two sets of programs are available: A. Main programs (Figures 1-5, Table 1); B. Data file and auxiliary programs. All the programs should be put under the same folder and run in Matlab. The auxiliary programs should not be run independently.

---
A.	Main programs

A.1 Simulation Study (Figure 1-3, Table 1) : “simstudymain.m”

A.2 Application to real data (Figure 4 and 5): “applicationmain.m” (require data file “stockpricedata.csv”)

---
B.	Data file and auxiliary programs

B.1 The data file “stockpricedata.csv” contains the daily market indices.

B.2 The list of auxiliary programs is as follows.

“cloverrnd.m”: Generates random samples of bivariate clover distribution

“CmpErr.m”: Calculation of the Estimation Relative Errors

“drawMulCon”: Draw the bivariate quantile contours

“elliprnd.m”: Generates random samples of the bivariate elliptical distribution (Section 3)

“genrest.m”: Generates radius estimate for different values of p or beta

“hallinQ.m”: Generates an estimate of the regression coefficients in Hallin, Paindaveine and Siman (2010) for a single scenario

“hddepth.m”: Calculation of the halfspace depth

“MrvHDOHPS.m”: Generates estimates of quantile regions using HPS approach (Hallin, Paindaveine and Siman, 2010) for all scenarios

“MrvHDOKM.m”: Generates estimates of quantile regions using KM approach (Kong and Mizera, 2012) for all scenarios

“MrvHDQEVT.m”: Generates estimates of quantile regions using our extreme value approach for all scenarios

“MrvHDQNpar.m”: Generates empirical estimates of quantile regions for all scenarios

“Q_hat.m”: Generate an extreme estimate of quantile regions for a single scenario

“rq.m”: Quantile regression program provided by Roger Koenker
