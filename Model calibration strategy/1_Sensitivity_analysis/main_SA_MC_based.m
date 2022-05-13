%%=======================================================================%%
%  Author: Rafael da Silva Raqueti 
%  Advisor: Samuel da Silva 
%  On the Calibration of Reduced-Order Models to Describe the 
%  Viscoelasticity in Steady-State Rolling Tires 
%  Methodology >> 1_Sensitivity_analysis
%    main_SA_MC_based.m
%%=======================================================================%%

%  main_SA_MC_based.m
%    This is the main file to perform a prior global sensitivity analysis
%    by calculating Sobol' indices via Monte Carlo simulation using UQLab
%    framework.

%__________________________________________________________________________
%% 0 - INITIALIZE UQLAB

clearvars
rng(100,'twister')
uqlab

%__________________________________________________________________________
%% 1 - MICHELIN DATASET PRE-PROCESSING


%__________________________________________________________________________
%% 2 - COMPUTATIONAL MODEL

% More information about the model module is available at
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

%__________________________________________________________________________
%% 4 - SENSITIVITY ANALYSIS

% More information about the Sensitivity Analysis module is available at
% https://www.uqlab.com/pce-user-manual

SobolOpts.Type = 'Sensitivity';         % UQLab sensitivity module;
SobolOpts.Method = 'Sobol';             % Sobol' indices;
SobolOpts.Input = myInput;              % Specification of the used input object;
SobolOpts.Model = myModel;              % Specification of the used model;
SobolOpts.Sobol.SampleSize = 1.0e+02;   % # of model evaluations.

% Create a sensitivity analysis:
tic
mySobolAnalysis = uq_createAnalysis(SobolOpts);
toc

uq_print(mySobolAnalysis)   % Printing;
uq_display(mySobolAnalysis) % Graphically display the results.
return

%__________________________________________________________________________
%% 5 - POST-PROCESSING RESULTS

% Total Sobol' indices:

figure1 = figure;
axes1 = axes('Parent',figure1);
hold(axes1,'on');
bar(mySobolAnalysis.Results.Total);
ylabel('$S^{\textnormal{Tot}}_{i}$','Interpreter','latex');
xlabel({'Parameters'},'Interpreter','latex');
title({'Total Sobol'' indices'},'Interpreter','latex');
box(axes1,'on');
hold(axes1,'off');
set(axes1,'FontSize',12,'TickLabelInterpreter','latex','XTick',...
    [1 2 3 4 5 6],'XTickLabel',...
    {'$c$','$k$','$\alpha$','$\gamma$','$\delta$','$\nu$'},'YGrid','on');


% Sobol' indices Order 1

figure2 = figure;
axes2 = axes('Parent',figure2);
hold(axes2,'on');
bar(mySobolAnalysis.Results.FirstOrder);
ylabel('$S^{(1)}_{u}$','Interpreter','latex');
xlabel({'Parameters'},'Interpreter','latex');
title({'Sobol'' indices Order 1'},'Interpreter','latex');
box(axes2,'on');
hold(axes2,'off');
set(axes2,'FontSize',12,'TickLabelInterpreter','latex','XTick',...
    [1 2 3 4 5 6],'XTickLabel',...
    {'$c$','$k$','$\alpha$','$\gamma$','$\delta$','$\nu$'},'YGrid','on');
