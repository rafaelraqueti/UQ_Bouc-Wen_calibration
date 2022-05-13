%%=======================================================================%%
%  Author: Rafael da Silva Raqueti
%  Advisor: Samuel da Silva
%  On the Calibration of Reduced-Order Models to Describe the 
%  Viscoelasticity in Steady-State Rolling Tires 
%  Methodology >> 3_Statistical_inference
%    getcandidate.m
%%=======================================================================%%

function Xcd = getcandidate(X,MHParameters)
%
%  getcandidate
%    This function generates a candidade for Metropolis-Hastings algorithm.
%
%  USAGE: Xcd = getcandidate(X,MHParameters)
%__________________________________________________________________________
%  OUTPUTS
%    Xcd : Vector of unknown parameters candidate.
%__________________________________________________________________________
%  INPUTS
%    X : Vector of parameters:
%      - X(1) = c;
%      - X(2) = k;
%      - X(3) = \alpha;
%      - X(4) = \nu.
%    MHParameters : Metropolis-Hastings static configuration parameters:
%      - MHParameters.Xmin : Lower bound;
%      - MHParameters.Xmax : Upper bound;
%      - MHParameters.sigma : Random walk step size.

%__________________________________________________________________________
% 1. Retrieve the static configuration parameters:
Xmin = MHParameters.Xmin;   % Lower and;
Xmax = MHParameters.Xmax;   % Upper bounds;
sigma = MHParameters.sigma; % Random walk step size;

Nl = size(Xmin,2); % Number of parameters.

% 2. Generate a candidade:
x = (X-Xmin)./(Xmax-Xmin);  % Dimensionless precedent value;
xcd = x+sigma*randn(1,Nl);  % Get a candidate;
% Check if all candidate values are in bounds:
if any(xcd < 0) || any(xcd > 1)
    cond = true;
    while cond
        xcd = x + sigma*randn(1,length(Xmin));
        if any(xcd < 0) || any(xcd > 1)
            cond = true;
        else
            cond = false;
            % break
        end
    end
end
% Return the candidate:
Xcd = (Xmax-Xmin).*xcd+Xmin;
end