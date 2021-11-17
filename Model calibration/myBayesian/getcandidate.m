%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rafael da Silva Raqueti
% Bruna Silveira Pavlack
% João Pedro Canisso Valese Norenberg
% Luccas Pereira Miguel
% Uncertainty Quantification
% Model calibration >> myBayesian
%  - getcandidate.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This script gets a candidade for Metropolis-Hastings algorithm.

function Xcd = getcandidate(X,MHParameters)

% 1. Retrieve the static configuration parameters:
Xmin = MHParameters.Xmin;           % Lower and;
Xmax = MHParameters.Xmax;           % Upper bounds;
sigma = MHParameters.sigma;         % Random walk step size;
Nl = size(Xmin,2);                  % Number of parameters.

% 2. Generate a candidade:
x = (X - Xmin)./(Xmax - Xmin);      % Dimensionless precedent value;
xcd = x + sigma*randn(1,Nl);        % Get a candidate;
% Check if all candidate values are in bounds:
if any(xcd < 0) || any(xcd > 1)
    cond = true;
    while cond
        xcd = x + sigma*randn(1,length(Xmin));
        if any(xcd < 0) || any(xcd > 1)
            cond = true;
        else
            cond = false;
            break
        end
    end
end

% Return the candidate:
Xcd = (Xmax - Xmin).*xcd + Xmin;

end

