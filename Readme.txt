Requirements: 
MATLAB R2020b, MATLAB Global Optimization Toolbox, MATLAB Statistics and Machine Learning Toolbox, and EnergyPlus 9.4.0;

%%%%%%%%%%%%%%%%%%%%%% Descriptions of folders %%%%%%%%%%%%%%%%%%%%%%%%%%
Example1 gives MATLAB figure files of all figures for Example 1 shown in the paper.
Example2 gives MATLAB figure files of all figures for Example 2 shown in the paper.
Example3 gives MATLAB figure files of all figures for Example 3 shown in the paper.

%%%%%%%%%%%%%%%%%%%%%% Descriptions of files %%%%%%%%%%%%%%%%%%%%%%%%%%%
1.	MainExample1.m — Replicates the simulations and reproduces all results for Example 1 except Figure H.4.
2.	MainExample1Size2.m — Replicates the simulations that provide the data used to plot Figure H.4 of Example 1 and reproduces Figure H.4.
3.	MainExample2.m — Replicates the simulations and reproduces all results for Example 2 except Figure I.3.
4.	MainExample2Size2.m — Replicates the simulations that provide the data used to plot Figure I.3 of Example 2 and reproduces Figure I.3.
5.	MainExample3.m — Replicates the simulations and reproduces all results for Example 3 except Figure F.6.
6.	MainExample3Size2.m — Replicates the simulations that provide the data used to plot Figure F.6 of Example 3 and reproduces Figure F.6.

Note that the above six scripts obtained their important outputs for each method from the six scripts (scripts no. 7 to no. 12) described below. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
7.	CalibrationAGP.m — Implements the bi-fidelity MBC-AGP, BC-AGP, MID-AGP or SR-AGP method. The table output of this script, i.e., RecordTable=table(D,Level,SVec,MaxAFVals,xhatstarMLs,Shminhats,ShxhatstarMLs,sigma2ls,thetals,rhos,gamma2s,sigma2hs,thetahs,phis,mus,condVprimes,MinM2LogLikelihoods), contains the following important variables:
(a) xhatstarMLs (n0-th row to last row, where n0 is the total number of initial design points): Gives the estimates of the MLE of the calibration parameter vector after data from the initial design and each follow up design point are obtained, where the n-th row of the matrix gives the estimate obtained with n design points.
(b) ShxhatstarMLs (n0-th row to last row): Gives the value of the true HF SSE evaluated at each row of xhatstarMLs, where the n-th element of the vector gives the true HF SSE evaluated at the n-th row of xhatstarMLs.
(c) D: Gives all design points, i.e., the initial design points, and all follow-up design points. All design points are represented by row vectors; follow-up design points are ordered chronologically by the time they were run, with earlier points listed as rows above later points; the rows representing the initial design points are listed above those representing the follow-up design points.
(d) Level: Gives the fidelity level of each design point (2 indicates a HF run and 1 indicates a LF run).
(e) MaxAFVals (n0-th row to last row): Gives the AF value at each follow-up design point (which maximizes the AF) at the row corresponding to the follow-up design point.
(f) Svec: Gives, at each HF design point in D, the true value of the HF SSE at that point, and at each LF design point in D, the true value of the LF LSSE (for the BC-AGP and SR-AGP methods) or modified LF SSE (for the MBC-AGP and MID-AGP methods) at that point.
(g) Shminhats (n0-th row to last row): The n-th element of this vector gives the 0.9 quantile of the HF SSE at the n-th row of xhatstarMLs.
Note: This code can only work if the HF SSE and LF SSE functions are both bounded from below by a positive number over the design region.

8.	CalibrationBCGP.m — Implements the single fidelity BC-GP method. The table output of this script, i.e.,
RecordTable=table(D,Level,ShVec,MaxAFVals,xhatstarMLs,Shminhats,ShxhatstarMLs,sigma2s,phis,thetas,mus,condRs,MinM2LogLikelihoods), contains the following important variables:
(a) xhatstarMLs (n0-th row to last row, where n0 is the total number of initial design points): Gives the estimates of the MLE of the calibration parameter vector after data from the initial design and each follow up design point are obtained, where the n-th row of the matrix gives the estimate obtained with n design points.
(b) ShxhatstarMLs (n0-th row to last row): Gives the value of the true HF SSE evaluated at each row of xhatstarMLs, where the n-th element of the vector gives the true HF SSE evaluated at the n-th row of xhatstarMLs.
(c) D: Gives all design points, i.e., the initial design points, and all follow-up design points. All design points are represented by row vectors; follow-up design points are ordered chronologically by the time they were run, with earlier points listed as rows above later points; the rows representing the initial design points are listed above those representing the follow-up design points.
(d) Level: Gives the fidelity level of each design point (2 indicates a HF run and 1 indicates a LF run).
(e) MaxAFVals (n0-th row to last row): Gives the AF value at each follow-up design point (which maximizes the AF) at the row corresponding to the follow-up design point.
(f) ShVec: Gives, at each design point in D, the true value of the HF SSE at that point.
(g) Shminhats (n0-th row to last row): The n-th element of this vector gives the posterior 0.9 quantile of the HF SSE at the n-th row of xhatstarMLs.
Note: This code can only work if the HF SSE function is bounded from below by a positive number over the design region.

9.	CalibrationSRGP.m — Implements the single fidelity SR-GP method. The table output of this script, i.e.,
RecordTable=table(D,Level,ZhVec,MaxAFVals,xhatstarMLs,Zhminhats,ShxhatstarMLs,sigma2s,thetas,mus,condRs,MinM2LogLikelihoods), contains the following important variables:
(a) xhatstarMLs (n0-th row to last row, where n0 is the total number of initial design points): Gives the estimates of the MLE of the calibration parameter vector after data from the initial design and each follow up design point are obtained, where the n-th row of the matrix gives the estimate obtained with n design points.
(b) ShxhatstarMLs (n0-th row to last row): Gives the value of the true HF SSE evaluated at each row of xhatstarMLs, where the n-th element of the vector gives the true HF SSE evaluated at the n-th row of xhatstarMLs.
(c) D: Gives all design points, i.e., the initial design points, and all follow-up design points. All design points are represented by row vectors; follow-up design points are ordered chronologically by the time they were run, with earlier points listed as rows above later points; the rows representing the initial design points are listed above those representing the follow-up design points.
(d) Level: Gives the fidelity level of each design point (2 indicates a HF run and 1 indicates a LF run).
(e) MaxAFVals (n0-th row to last row): Gives the AF value at each follow-up design point (which maximizes the AF) at the row corresponding to the follow-up design point.
(f) ZhVec: Gives, at each point in D, the square root of the true value of the HF MSE (HF SSE divided by the number of field observations N) at that point.
(g) Zhminhats (n0-th row to last row): The n-th element of this vector gives the posterior 0.9 quantile of the square root of the HF MSE evaluated at the n-th row of xhatstarMLs.

10.	CalibrationNested.m — Implements the bi-fidelity Nested method. The table output of this script, i.e.,
RecordTable=table(D,Level,ZVec,MaxAFVals,xhatstarMLs,Zhminhats,ShxhatstarMLs,sigma2ls,thetals,rhos,sigma2hs,thetahs,mus,condRlRhs,MinM2LogLikelihoodhs,MinM2LogLikelihoodls,MinM2LogLikelihoods), contains the following important variables:
(a) xhatstarMLs (n0-th row to last row, where n0 is the total number of initial design points): Gives the estimates of the MLE of the calibration parameter vector after data from the initial design and each follow up design point are obtained, where the n-th row of the matrix gives the estimate obtained with n design points.
(b) ShxhatstarMLs (n0-th row to last row): Gives the value of the true HF SSE evaluated at each row of xhatstarMLs, where the n-th element of the vector gives the true HF SSE evaluated at the n-th row of xhatstarMLs.
(c) D: Gives all design points, i.e., the initial design points, and all follow-up design points. All design points are represented by row vectors; follow-up design points are ordered chronologically by the time they were run, with earlier points listed as rows above later points; the rows representing the initial design points are listed above those representing the follow-up design points.
(d) Level: Gives the fidelity level of each design point (2 indicates a HF run and 1 indicates a LF run).
(e) MaxAFVals (n0-th row to last row): Gives the AF value at each follow-up design point (which maximizes the AF) at the row corresponding to the follow-up design point.
(f) ZVec: Gives, at each design point in D, the square root of the true value of the HF SSE at that point if that point is a HF design point, or the square root of the true value of the LF SSE at that point if that point is a LF design point.
(g) Zhminhats (n0-th row to last row): The n-th element of this vector gives the smallest square root of the HF SSE value observed at a HF design point after n-th experiment runs have been made.
Note: Nested method only works if the vector of square roots of LF SSE values at points in the intersection of the initial HF and LF designs does not have identical elements.

11.	CalibrationSVD.m — Implements the single fidelity SVD method. The table output of this script, i.e., RecordTable=table(D,Level,SSEs,MaxAFVals,xhatstarMLs,Shminhats,ShxhatstarMLs,ps,sigma2s,thetas,CondRs,MinM2LogLikelihoods), contains the following important variables:
(a) xhatstarMLs (n0-th row to last row, where n0 is the total number of initial design points): Gives the estimates of the MLE of the calibration parameter vector after data from the initial design and each follow up design point are obtained, where the n-th row of the matrix gives the estimate obtained with n design points.
(b) ShxhatstarMLs (n0-th row to last row): Gives the value of the true HF SSE evaluated at each row of xhatstarMLs, where the n-th element of the vector gives the true HF SSE evaluated at the n-th row of xhatstarMLs.
(c) D: Gives all design points, i.e., the initial design points, and all follow-up design points. All design points are represented by row vectors; follow-up design points are ordered chronologically by the time they were run, with earlier points listed as rows above later points; the rows representing the initial design points are listed above those representing the follow-up design points.
(d) Level: Gives the fidelity level of each design point (2 indicates a HF run and 1 indicates a LF run).
(e) MaxAFVals (n0-th row to last row): Gives the AF value at each follow-up design point (which maximizes the AF) at the row corresponding to the follow-up design point.
(f) SSEs: Gives, at each point in D, the HF SSE value at that point.
(g) Shminhats (n0-th row to last row): The n-th element of this vector gives the posterior mean of the HF SSE evaluated at the n-th row of xhatstarMLs.

12.	CalibrationSVDAGP.m — Implements the bi-fidelity SVD-AGP method. The table output of this script, i.e.,
RecordTable=table(D,Level,SSEs,MaxAFVals,xhatstarMLs,Shminhats,ShxhatstarMLs,pls,phs,thetals,sigma2ls,thetahs,sigma2hs,condRlRhs,MinM2LogLikelihoodhs,MinM2LogLikelihoodls), contains the following important variables:
(a) xhatstarMLs (n0-th row to last row, where n0 is the total number of initial design points): Gives the estimates of the MLE of the calibration parameter vector after data from the initial design and each follow up design point are obtained, where the n-th row of the matrix gives the estimate obtained with n design points.
(b) ShxhatstarMLs (n0-th row to last row): Gives the value of the true HF SSE evaluated at each row of xhatstarMLs, where the n-th element of the vector gives the true HF SSE evaluated at the n-th row of xhatstarMLs.
(c) D: Gives all design points, i.e., the initial design points, and all follow-up design points. All design points are represented by row vectors; follow-up design points are ordered chronologically by the time they were run, with earlier points listed as rows above later points; the rows representing the initial design points are listed above those representing the follow-up design points.
(d) Level: Gives the fidelity level of each design point (2 indicates a HF run and 1 indicates a LF run).
(e) MaxAFVals (n0-th row to last row): Gives the AF value at each follow-up design point (which maximizes the AF) at the row corresponding to the follow-up design point.
(f) SSEs: Gives, at each LF design point in D, the LF SSE value at that point, and at each HF design point in D, the HF SSE value at that point.
(g) Shminhats (n0-th row to last row): The n-th element of this vector gives the posterior mean of the HF SSE evaluated at the n-th row of xhatstarMLs.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
13.	ComputeRmatrix.m — Computes a matrix of the prior correlations between two sets of points given by the product Matern correlation function with smoothness parameter 3/2. 
14.	ComputeRmatrix2.m — Computes the prior correlation matrix for a set of points given by the product Matern correlation function with smoothness parameter 3/2 with a nugget added to each diagonal element of the matrix. 
15.	EPBasicFile.idf — Gives a basic EnergyPlus input editor file for running EnergyPlus simulations.
16.	EPWeather.epw — Includes Hong Kong weather data for use in the EnergyPlus simulations.
17.	Example1FindTrueMLE — Finds the true MLE of the calibration parameter vector for Example 1.
18.	Example1GenerateGridData5forAppendixH — Computes the HF and LF outputs of the simulator for Example 1 at the points in the 26^3 grid {0,1/25,...,1}^3, and records the run times for all those simulation runs (the data are used to compute various integrals whose values are needed for the results reported in Appendix H, and the ratio of the average run times for the HF simulations to that for the LF simulations can be used to justify setting c_h/c_l=4).
19.	Example2FindTrueMLE — Finds the true MLE of the calibration parameter vector for Example 2.
20.	GenerateNestedLHD.m — Generates one pair of nested LHDs using the method described in Qian (2009). Reference: Qian, P.Z.G. (2009) Nested Latin hypercube designs. Biometrika, 96(4), 957-970.
21.	invandlogdet.m — Gives the inverse, log determinant, and condition number of a positive definite matrix.
22.	Simulator.m — Gives the HF and LF simulator outputs for Examples 1-3.
23.	TransformData.m — Gives the Box-Cox, square root, or identity transformation, and the log absolute value of the Jacobian for the transformation.
24.	TransformData_inv.m — Gives the inverse of the Box-Cox, square root, or identity transformation.
25.	ApproximateSh_SVD.m — Computes, for the 100 trials in Example 2, the true minimizer, minimum, and range over {0,0.05,…,1}^3 of S_{h,1}, the S_{h,1} value at the estimate of the calibration parameter vector given by the SVD method at termination, and other information needed to produce some numerical results stated in Appendix J. 
26.	ApproximateSh_SVDAGP.m — Computes, for the 100 trials in Example 2, the true minimizer, minimum, and range over {0,0.05,…,1}^3 of S_{h,2}, the S_{h,2} value at the estimate of the calibration parameter vector given by the SVD-AGP method at termination, and other information needed to produce some numerical results stated in Appendix J.
27.	AnalyzeApproximateShData.m — Produces the numbers reported in Appendix J using the data obtained by running ApproximateSh_SVD.m and ApproximateSh_SVDAGP.m.

28.	Example1.mat — Contains the numerical results for Example 1.
29.	Example1InputData.mat — Contains the initial designs and initial data for Example 1.
30.	Example1Size2.mat — Contains the numerical results used to create Figure H.4 of Example 1. 
31.	Example1Size2InputData.mat — Contains the initial designs and initial data used to create Figure H.4 of Example 1.
32.	Example1GridData26.mat — Gives the LF and HF simulator output data on a 26^3-point grid for Example 1.
33.	Example2.mat — Contains the numerical results for Example 2.
34.	Example2InputData.mat — Contains the initial designs and initial data for Example 2.
35.	Example2Size2.mat — Contains the numerical results used to create Figure I.3 of Example 2.
36.	Example2Size2InputData.mat — Contains the initial designs and initial data used to create Figure I.3 of Example 2.
37.	Example3.mat — Contains the numerical results for Example 3.
38.	Example3InputData.mat — Contains the initial designs and initial data for Example 3.
39.	Example3Size2.mat — Contains the numerical results used to create Figure F.6 of Example 3. 
40.	Example3Size2InputData.mat — Contains the initial designs and initial data used to create Figure F.6 of Example 3.
41.     DataforAppendixJ_SVD.mat — Contains data used to obtain numerical results for the SVD method reported in Appendix J.
42.     DataforAppendixJ_SVDAGP.mat — Contains data used to obtain numerical results for the SVD-AGP method reported in Appendix J.
43.     Example1TrueMLE.mat — Contains the true MLE of the calibration parameter vector and value of the HF SSE at the true MLE of the calibration parameter vector for Example 1.
44.     Example2TrueMLE.mat — Contains the true MLE of the calibration parameter vector and value of the HF SSE at the true MLE of the calibration parameter vector for Example 2. 

%%%%%%%%%%%%%%%%%% Steps for replicating the results for Example 1 %%%%%%%%%%%%%%%%%%%
Step 1: Install EnergyPlus 9.4.0 in the folder C:\EnergyPlusV9-4-0.
Step 2: Put all the MATLAB codes, the EPBasicFile.idf file, and the EPWeather.epw file in a single folder.
Step 3: Run MainExample1.m and MainExample1Size2.m.

%%%%%%%%%%%%%%%%%% Step for replicating the results for Example 2 %%%%%%%%%%%%%%%%%%%
Step 1: Put all the MATLAB codes in a single folder.
Step 2: Run MainExample2.m and MainExample2Size2.m.

%%%%%%%%%%%%%%%%%% Step for replicating the results for Example 3 %%%%%%%%%%%%%%%%%%%
Step 1: Put all the MATLAB codes in a single folder.
Step 2: Run MainExample3.m and MainExample3Size2.m.

%%%%%%%%%%%%%% Step for obtaining the results at the end of Appendix J %%%%%%%%%%%%%%%%%%%
Steps: Run ApproximateSh_SVD.m and then ApproximateSh_SVDAGP.m. Finally, run AnalyzeApproximateShData.m.

Note: 
To check for differences in numerical computation results between computers, you may want to check if the simulator output values in the following data files remain the same on your computer: Example1InputData.mat, Example1Size2InputData.mat, Example1GridData26, Example2InputData.mat, Example2Size2InputData.mat, Example2TrueMLE.mat, Example3InputData.mat, Example3Size2InputData.mat.
For example, the values in MultiDataInput(id).Yl, MultiDataInput(id).Yh, and SingleDataInput(id).Yh contain simulator outputs at the points in MultiDataInput(id).Dl, MultiDataInput(id).Dh, and SingleDataInput(id).Dh for id=1 to id=100. If the output values are different, you may want to recompute all the output values.
