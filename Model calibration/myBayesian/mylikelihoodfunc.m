%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rafael da Silva Raqueti
% Bruna Silveira Pavlack
% João Pedro Canisso Valese Norenberg
% Luccas Pereira Miguel
% Uncertainty Quantification
% Model calibration >> myBayesian
%  - mylikelihoodfunc.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This spcrip evaluates the likelihood function.

function [Y,L] = mylikelihoodfunc(X,Parameters)

% 1. Retrieve the static configuration parameters:
t = Parameters.t;                   % Time vector;
dataA = Parameters.dataA;           % Data of viscoelastic internal variables;
dataC = Parameters.dataC;           % Data of right Cauchy-Green deformation;
epsilonCE = Parameters.epsilonCE;   % Cross-Entropy Method MASE measure;
Ni = size(dataA,1);                 % Number of time points;
Nj = size(dataA,2);                 % Number of selected data.

% 2. Evaluate the likelihood function:
% This script cannot be parallelized! Candidate depends on previous value.
epsilon = zeros(1,Nj);      % Pre-allocating for MASE values;
Y = zeros(Ni,Nj);           % Pre-allocating for reduced-order model responses;
options = odeset('RelTol',1e-06,'AbsTol',1e-08);
for j = 1:Nj
    % Solve reduced-order model:
    C = griddedInterpolant(t,dataC(:,j));   % Input C(t);
    IC = [dataA(1,j);0];                    % Initial conditions.
    [~,y] = ode45(@(t,Y) myode(t,Y,C,X(1,1:end-1)),t,IC,options);
    Y(:,j) = y(:,1);
    % Compute MASE denominator:
    % Atention! It can be possible to divide by zero.
    den = sum(abs(dataA(2:end,j)-dataA(1:end-1,j)),1)/(Ni-1);
    % Evaluate MASE:
    epsilon(j) = sum(abs((dataA(:,j)-y(:,1))/den),1);
end
epsilon = sum(epsilon/Ni,2);

% Evaluate and return the likelihood function:
L = exp(-0.5*(epsilon-epsilonCE)^2/X(end))/sqrt(X(end));
% In which X(end) is the discrepancy parameter.

end
