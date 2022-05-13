%%=======================================================================%%
%  Author: Rafael da Silva Raqueti
%  Advisor: Samuel da Silva
%  On the Calibration of Reduced-Order Models to Describe the 
%  Viscoelasticity in Steady-State Rolling Tires 
%  Methodology >> 3_Statistical_inference
%    myLikelihood.m
%%=======================================================================%%

function [Y,L] = myLikelihood(X,Parameters)
%
%  myLikelihood
%    This function evaluates the likelihood function.
%
%  USAGE: [Y,L] = myLikelihood(X,Parameters)
%__________________________________________________________________________
%  OUTPUTS
%    Y : State vector [A, Z] where:
%      - A is the viscoelastic internal variable;
%      - Z is the histeretic output.
%    L : Likelihood function value.
%__________________________________________________________________________
%  INPUTS
%    X : Vector of unknown parameters:
%      - X(1) = c;
%      - X(2) = k;
%      - X(3) = \alpha;
%      - X(4) = \nu.
%    Parameters : Computational model static configuration parameters:
%      - Parameters.t : Time vector;
%      - Parameters.dataA : Viscoelastic internal variables data;
%      - Parameters.dataC : Components of right Cauchy-Green tensor data;
%      - Parameters.epsilonCE : Cross-Entropy method MASE value.
%__________________________________________________________________________
% 1. Retrieve the static configuration parameters:
t = Parameters.t;
dataA = Parameters.dataA;
dataC = Parameters.dataC;
epsilonCE = Parameters.epsilonCE;

Ni = size(dataA,1); % # of time samples;
Nj = size(dataA,2); % # of selected viscoelastic internal variables data.

% 2. Evaluate the Likelihood function:
Y = zeros(Ni,Nj);      % Preallocating Reduced-Order Model responses;
epsilon = zeros(1,Nj); % Preallocating MASE values;
C = cell(1,Nj);        % Preallocate C(t);
IC = zeros(2,Nj);      % Preallocate Initial Condition;
den = zeros(1,Nj);     % Preallocate MASE denomiator.
parfor j = 1:Nj
    % Assign input C(t) and Initial Condition:
    C{1,j} = griddedInterpolant(t,dataC(:,j));
    IC(1,j) = dataA(1,j);
    % Compute MASE denominator:
    den(1,j) = sum(abs(dataA(2:end,j)-dataA(1:end-1,j)),1)/(Ni-1);
end
% Define odeset options:
% options = odeset('RelTol',1e-06,'AbsTol',1e-08);
parfor j = 1:Nj
    % Solve the Reduced-Order Model:
    % [~,y] = ode45(@(t,y) myODE2(t,y,C{j},X(1,1:end-1)),t,IC(:,j),options);
    [~,y] = ode45(@(t,y) myODE2(t,y,C{j},X(1,1:end-1)),t,IC(:,j));
    Y(:,j) = y(:,1);
    % Evaluate MASE:
    epsilon(j) = sum(abs((dataA(:,j)-y(:,1))/den(j)),1);
end
% Compute the mean of MASE values:
epsilon = sum(epsilon/Ni,2)/Nj;
% Evaluate and return the Likelihood function:
L = exp(-0.5*(epsilon-epsilonCE)^2/X(end))/sqrt(X(end));
% where X(end) is the discrepancy parameter.
end