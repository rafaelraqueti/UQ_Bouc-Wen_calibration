%%=======================================================================%%
%  Author: Rafael da Silva Raqueti
%  Advisor: Samuel da Silva
%  On the Calibration of Reduced-Order Models to Describe the 
%  Viscoelasticity in Steady-State Rolling Tires 
%  Methodology >> 2_Optimization
%    myMASE.m
%%=======================================================================%%

function epsilon = myMASE(X,Parameters)
%
%  myMASE
%    This function returns the Mean Absolute Scaled Error (MASE) mean
%    between viscoelastic internal variables and Reduced-Order Model
%    responses. The objective is to determine the optimal set of parameters
%    that minimizes this value.
%
%  USAGE: epsilon = myMASE(X,Parameters)
%__________________________________________________________________________
%  OUTPUTS
%    epsilon : Mean of MASE values;
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
%      - Parameters.dataC : Components of right Cauchy-Green tensor data.
%__________________________________________________________________________
% 1. Retrieve the static configuration parameters:
t = Parameters.t;
dataA = Parameters.dataA;
dataC = Parameters.dataC;

Ni = size(dataA,1); % # of time samples;
Nj = size(dataA,2); % # of selected viscoelastic internal variables data;
Nk = size(X,1);     % # of unknown parameters samples.

% 2. Calculate the model responses for each unknown parameters sample:
epsilon = zeros(Nk,Nj); % Preallocate MASE values;
C = cell(1,Nj);         % Preallocate C(t);
IC = zeros(2,Nj);       % Preallocate Initial Condition;
den = zeros(1,Nj);      % Preallocate MASE denomiator.
parfor j = 1:Nj
    % Assign input C(t) and Initial Condition:
    C{1,j} = griddedInterpolant(t,dataC(:,j));
    IC(1,j) = dataA(1,j);
    % Compute MASE denominator:
    den(1,j) = sum(abs(dataA(2:end,j)-dataA(1:end-1,j)),1)/(Ni-1);
end
% Define odeset options:
% options = odeset('RelTol',1e-06,'AbsTol',1e-08);
parfor k = 1:Nk
    for j = 1:Nj
        % Solve the Reduced-Order Model:
        % [~,Y] = ode45(@(t,Y) myODE2(t,Y,C{j},X(k,:)),t,IC(:,j),options);
        [~,Y] = ode45(@(t,Y) myODE2(t,Y,C{j},X(k,:)),t,IC(:,j));
        % Evaluate MASE:
        epsilon(k,j) = sum(abs((dataA(:,j)-Y(:,1)))/den(j));
    end
end
% Return the mean of MASE values:
epsilon = sum(epsilon/Ni,2)/Nj;
end