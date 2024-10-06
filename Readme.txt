Requirements: 
Matlab R2020b, Matlab Global Optimization Toolbox, Matlab Statistics and Machine Learning Toolbox, and EnergyPlus 9.4.0;

%%%%%%%%%%%%%%%%%%%%%% Descriptions of folders%%%%%%%%%%%%%%%%%%%%%
Example1 gives MATLAB figure files of all figures for Example 1 shown in the paper.
Example2 gives MATLAB figure files of all figures for Example 2 shown in the paper.
Example3 gives MATLAB figure files of all figures for Example 3 shown in the paper.
FindMLEforExample1 gives all MATLAB codes and input data files used to find the true MLE of the calibration parameter vector for Example 1. Data files resulting from running the codes are also given.
FindMLEforExample2 gives all MATLAB codes and input data files used to find the true MLE of the calibration parameter vector for Example 2. Data files resulting from running the codes are also given.

%%%%%%%%%%%%%%%%%%%%%% Descriptions of files %%%%%%%%%%%%%%%%%%%%%
CalibrationAGP.m----------Implements the bi-fidelity MBC-AGP, BC-AGP, MID-AGP or SR-AGP method.
CalibrationBCGP.m----------Implements the single fidelity BC-GP method.
CalibrationNested.m----------Implements the bi-fidelity Nested method.
CalibrationSRGP.m----------Implements the single fidelity SR-GP method.
CalibrationSVD.m ----------Implements the single fidelity SVD method.
CalibrationSVDAGP.m ----------Implements the bi-fidelity SVD-AGP method.
ComputeRmatrix.m----------Computes a matrix of the prior correlations between two sets of points given by the product Matern correlation function with smoothness parameter 3/2. 
ComputeRmatrix2.m----------Computes the prior correlation matrix for a set of points given by the product Matern correlation function with smoothness parameter 3/2 with a nugget added to each diagonal element of the matrix. 
EPBasicFile.idf----------Gives a basic EnergyPlus input editor file for running EnergyPlus simulations.
EPWeather.epw---------- Includes Hong Kong weather data for use in the EnergyPlus simulations.
GenerateNestedLHD.m----------Generates one pair of nested LHDs using the method described in Qian (2009). Reference: Qian, P. Z. G. (2009). Nested Latin hypercube designs. Biometrika, 96(4), 957–970.
MainExample1.m----------Replicates the simulations and gives all results for Example 1.
MainExample1Size2.m ----------Replicates the simulations and gives all results for Example 1, larger initial design size.
MainExample2.m----------Replicates the simulations and gives all results for Example 2.
MainExample2Size2.m----------Replicates the simulations and gives all results for Example 2, larger initial design size.
MainExample3.m----------Replicates the simulations and gives all results for Example 3.
MainExample3Size2.m ----------Replicates the simulations and gives all results for Example 3, smaller initial design size.
Simulator.m----------Gives the HF and LF simulator outputs for Examples 1-3.
TransformData.m----------Gives the Box-Cox, square root, or identity transformation, and the log absolute value of the Jacobian for the transformation.
TransformData_inv.m----------Gives the inverse of the Box-Cox, square root, or identity transformation.
invandlogdet.m----------Gives the inverse and log determinant of a positive definite matrix.

Example1.mat----------Initial designs and numerical results for Example 1.
Example1Size2.mat ----------Initial designs and numerical results for Example 1, larger initial design size 
Example1GridData5.mat----------Gives the LF and HF simulator output data on a 26^3-point grid for Example 1.
Example2.mat----------Initial designs and numerical results for Example 2.
Example2Size2.mat ----------Initial designs and numerical results for Example 2, larger initial design size 
Example3.mat----------Initial designs and numerical results for Example 3.
Example3Size2.mat ----------Initial designs and numerical results for Example 3, smaller initial design size 

%%%%%%%%%%%%%%%%%%%%% Steps for replicating results for Example 1 %%%%%%%%%%%%%%%%%%%%%
Step 1: Install EnergyPlus 9.4.0 in the folder C:\EnergyPlusV9-4-0
Step 2: Put all the Matlab codes, the EPBasicFile.idf file, and the EPWeather.epw file in a single folder.
Step 3: Run MainExample1.m and MainExample1Size2.m

%%%%%%%%%%%%%%%%%%%%% Step for replicating results for Example 2 %%%%%%%%%%%%%%%%%%%%%
Step: Run MainExample1.m and MainExample2Size2.m to obtain the results for Example 2

%%%%%%%%%%%%%%%%%%%%% Step for replicating results for Example 3 %%%%%%%%%%%%%%%%%%%%%
Step: Run MainExample3.m and MainExample3Size2.m to obtain the results for Example 3
