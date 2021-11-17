%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rafael da Silva Raqueti
% Bruna Silveira Pavlack
% João Pedro Canisso Valese Norenberg
% Luccas Pereira Miguel
% Uncertainty Quantification
% Model calibration >> UQLab_Bayesian
%  - main_UQLab_Bayesian.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This spcrip is the main file for the bayesian inference.

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
ModelOpts.Parameters.dataC = tempC;
ModelOpts.Parameters.dataA = tempA;

% Create UQLab model object:
myForwardModel = uq_createModel(ModelOpts);


%% 3 - PRIOR DISTRIBUTION OF MODEL PARAMETERS

% More information about the input module is available at
% https://www.uqlab.com/input-user-manual

PriorOpts.Marginals(1).Name = 'c';
PriorOpts.Marginals(1).Type = 'Uniform';
PriorOpts.Marginals(1).Parameters = [+0.0e+00, +5.0e-03];

PriorOpts.Marginals(2).Name = 'k';
PriorOpts.Marginals(2).Type = 'Uniform';
PriorOpts.Marginals(2).Parameters = [+9.995e-01, 1.0005e+00];

PriorOpts.Marginals(3).Name = 'alpha';
PriorOpts.Marginals(3).Type = 'Uniform';
PriorOpts.Marginals(3).Parameters = [+0.0e+00, 1.0e+00];

PriorOpts.Marginals(4).Name = 'nu';
PriorOpts.Marginals(4).Type = 'Uniform';
PriorOpts.Marginals(4).Parameters = [+1.0e+00, 3.0e+00];

% Create UQLab input object:
myPriorDist = uq_createInput(PriorOpts);


%% 4 - SELECTED DATA

myData.y = epsilonCE;
myData.MOMap = [1;... % Model ID
                1];   % Output ID


%% 5 - DISCREPANCY MODEL

% Specify discrepancy distribution as an UQLab input object:
SigmaOpts.Marginals.Name = 'sigma2';
SigmaOpts.Marginals.Type = 'Uniform';
SigmaOpts.Marginals.Parameters = [+0.0e+00, +1.0e+00];

mySigmaDist = uq_createInput(SigmaOpts);

% Assign the distribution of $\sigma^2$ to the discrepancy options:
DiscrepancyOpts.Type = 'Gaussian';
DiscrepancyOpts.Prior = mySigmaDist;


%% 6 - BAYESIAN ANALYSIS

% More information about Bayesian inference is available at
% https://www.uqlab.com/inversion-user-manual

% Solver options:
Solver.Type = 'MCMC';
Solver.MCMC.Sampler = 'MH';
myProposal.PriorScale = 0.0225;     % For MH algorithm;
Solver.MCMC.Proposal = myProposal;  % For MH algorithm;
Solver.MCMC.NChains = 10; % Parallel chains
Solver.MCMC.Steps = +1.0e+04; % Iterations

% Posterior sample generation:
BayesOpts.Type = 'Inversion';
BayesOpts.Name = 'Bayesian model';
BayesOpts.Prior = myPriorDist;
BayesOpts.Data = myData;
BayesOpts.Discrepancy = DiscrepancyOpts;
BayesOpts.Solver = Solver;

% Run the Bayesian inversion analysis:
myBayesianAnalysis = uq_createAnalysis(BayesOpts);

% Print out a report of the results and create a graphical representation
% of the results:
uq_print(myBayesianAnalysis)
uq_display(myBayesianAnalysis)
