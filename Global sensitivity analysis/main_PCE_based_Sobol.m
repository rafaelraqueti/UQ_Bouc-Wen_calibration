%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rafael da Silva Raqueti
% Bruna Silveira Pavlack
% João Pedro Canisso Valese Norenberg
% Luccas Pereira Miguel
% Uncertainty Quantification
% Global sensitivity analysis
%  - main_PCE_based_Sobol.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This spcrip is the main file for the prior global sensitivity analysis.

% Before beginning:
clearvars % Clear all variables from Workspace;
rng(100,'twister') % Control the random number generator;
uqlab % And initialize UQLab.


%% 1 - DATASET

%% 2 - COMPUTATIONAL MODEL

% More information about the model module is available at
% https://www.uqlab.com/model-user-manual

% Define the computational model as an UQLab model object:
ModelOpts.mFile = 'myMASE';
ModelOpts.isVectorized = true;

% Set the static configuration parameters:
ModelOpts.Parameters.t = t;
ModelOpts.Parameters.dataA = dataA;
ModelOpts.Parameters.dataC = dataC;

% Create UQLab model object:
myModel = uq_createModel(ModelOpts);


%% 3 - INPUT PARAMETERS

% More information about the input module is available at
% https://www.uqlab.com/input-user-manual

InputOpts.Marginals(1).Name = 'c';
InputOpts.Marginals(1).Type = 'Uniform';
InputOpts.Marginals(1).Parameters = [0.0e+00, 1.0e-02];

InputOpts.Marginals(2).Name = 'k';
InputOpts.Marginals(2).Type = 'Uniform';
InputOpts.Marginals(2).Parameters = [9.99e-01, 1.001e+00];

InputOpts.Marginals(3).Name = 'alpha';
InputOpts.Marginals(3).Type = 'Uniform';
InputOpts.Marginals(3).Parameters = [0.0e+00, 1.0e+00];

InputOpts.Marginals(4).Name = 'gamma';
InputOpts.Marginals(4).Type = 'Uniform';
InputOpts.Marginals(4).Parameters = [1.0e+03, 1.0e+04];

InputOpts.Marginals(5).Name = 'delta';
InputOpts.Marginals(5).Type = 'Uniform';
InputOpts.Marginals(5).Parameters = [-1.0e+03, 1.0e+03];

InputOpts.Marginals(6).Name = 'nu';
InputOpts.Marginals(6).Type = 'Uniform';
InputOpts.Marginals(6).Parameters = [1.0e+00, 3.0e+00];

% Create UQLab input object:
myInput = uq_createInput(InputOpts);


%% 4 - PCE-BASED METAMODEL

% Select the metamodeling tool in UQLab and the polynomial chaos expansion
% (PCE) type:
PCEOpts.Type = 'Metamodel';
PCEOpts.MetaType = 'PCE';

% Assign the full computational model:
PCEOpts.FullModel = myModel;

% Specify the maximum polynomial degree:
PCEOpts.Method = 'LARS';
PCEOpts.Degree = 1:15;
PCEOpts.TruncOptions.qNorm = 0.75;

% Specify the size of the experimental design (i.e., the total
% computational cost of constructing the metamodel):
PCEOpts.ExpDesign.NSamples = 2.0e+03;

% Calculate the PCE:
myPCE = uq_createModel(PCEOpts);


%% 5 - SENSITIVITY ANALYSIS

% PCE-based Sobol' indices

SobolOpts.Type = 'Sensitivity';
SobolOpts.Method = 'Sobol';
SobolOpts.Input = myInput;
SobolOpts.Model = myPCE;
SobolOpts.Sobol.SampleSize = 1.0e+02;

% Create UQLab analysis object
mySobolAnalysis = uq_createAnalysis(SobolOpts);

% Print and display analysis
uq_print(mySobolAnalysis)
uq_display(mySobolAnalysis)

return
%% 5.1 - PCE-BASED METAMODEL VALIDATION

Xval = uq_getSample(1e+02);
Yval = uq_evalModel(myModel,Xval);
YPCE = uq_evalModel(myPCE,Xval);
