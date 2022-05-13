%%=======================================================================%%
%  Author: Rafael da Silva Raqueti
%  Advisor: Samuel da Silva
%  On the Calibration of Reduced-Order Models to Describe the 
%  Viscoelasticity in Steady-State Rolling Tires 
%  Methodology >> 1_Sensitivity_analysis
%    main_SA_PCE_based.m
%%=======================================================================%%

%  main_SA_PCE_based.m
%    This is the main file to perform a prior global sensitivity analysis
%    by calculating Sobol' indices from the coefficients of a Polynomial
%    Chaos Expansion (PCE) based metamodel using UQLab framework.

%__________________________________________________________________________
%% 0 - INITIALIZE UQLAB

clearvars
rng(100,'twister')
uqlab

%__________________________________________________________________________
%% 1 - MICHELIN DATASET PRE-PROCESSING


%__________________________________________________________________________
%% 2 - COMPUTATIONAL MODEL

% More information about the Model module is available at
% https://www.uqlab.com/model-user-manual

% Define the computational model as an UQLab model object:
ModelOpts.mFile = 'myMASE';
ModelOpts.isVectorized = true;

% Set computational model static configuration parameters:
ModelOpts.Parameters.t = t;
ModelOpts.Parameters.dataA = tempA;
ModelOpts.Parameters.dataC = tempC;

% Create UQLab model object:
myModel = uq_createModel(ModelOpts);

%__________________________________________________________________________
%% 3 - PROBABILISTIC INPUT MODEL

% More information about the Input module is available at
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

InputOpts.Marginals(4).Name = 'nu';
InputOpts.Marginals(4).Type = 'Uniform';
InputOpts.Marginals(4).Parameters = [1.0e+00, 3.0e+00];

% Create UQLab input object:
myInput = uq_createInput(InputOpts);

%__________________________________________________________________________
%% 4 - POLYNOMIAL CHAOS EXPANSION (PCE) METAMODEL

% More information about PCE metamodeling tool is available at
% https://www.uqlab.com/pce-user-manual

PCEOpts.Type = 'Metamodel';           % UQLab metamodeling tool;
PCEOpts.MetaType = 'PCE';             % Polynomial Chaos Expansion (PCE);
PCEOpts.Input = myInput;              % Probabilistic input model;
PCEOpts.Degree = 1:15;                % Polynomial degrees;
PCEOpts.TruncOptions.qNorm = 0.75;    % Hyperbolic truncation scheme;
PCEOpts.FullModel = myModel;          % Experimental design (model);
PCEOpts.ExpDesign.NSamples = 4.0e+03; % Experimental design (size).

% Create a PCE metamodel:
tic
myPCE = uq_createModel(PCEOpts);
toc

return
%% 

% Load Cross-Entropy method results:
load('...')

myData.y = epsilonCE;
myData.MOMap = [1; 1];


%%

SigmaOpts.Marginals.Type = 'Uniform';
SigmaOpts.Marginals.Parameters = [0, 1e-01];

SigmaDist = uq_createInput(SigmaOpts);

DiscrepancyOpts.Type = 'Gaussian';
DiscrepancyOpts.Prior = SigmaDist;


%%

Solver.Type = 'MCMC';
Solver.MCMC.Sampler = 'AIES';
Solver.MCMC.Steps = 2e+04;
Solver.MCMC.NChains = 10;

BayesOpts.Type = 'Inversion';
BayesOpts.Prior = myInput;
BayesOpts.ForwardModel = myPCE;
BayesOpts.Data = myData;
BayesOpts.Discrepancy = DiscrepancyOpts;
BayesOpts.Solver = Solver;

myBayesianAnalysis = uq_createAnalysis(BayesOpts);
uq_print(myBayesianAnalysis)
uq_display(myBayesianAnalysis)
% uq_display(myBayesianAnalysis,'trace','all')
% uq_display(myBayesianAnalysis,'meanConvergence','all')
% uq_display(myBayesianAnalysis,'acceptance',true)

return
%%

PostSample3D = myBayesianAnalysis.Results.Sample;
PostSample2D = reshape(permute(PostSample3D,[2 1 3]),size(PostSample3D,2),[]);

XHat = PostSample2D(1:6,:);

InputOptsHat.Inference.Data = XHat';
InputOptsHat.Copula.Type = 'Independent';

myInputHat = uq_createInput(InputOptsHat);
InputOpts.Marginals(6).Bounds = [1, inf];

uq_print(myInputHat)
uq_display(myInputHat)

%%

SobolOpts.Type = 'Sensitivity';
SobolOpts.Method = 'Sobol';

SobolOpts.Sobol.Order = 1;

SobolOpts.Sobol.SampleSize = 1e+05;

mySobolAnalysis = uq_createAnalysis(SobolOpts);

uq_print(mySobolAnalysis)
uq_display(mySobolAnalysis)
