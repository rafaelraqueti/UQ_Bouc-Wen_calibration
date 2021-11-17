%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rafael da Silva Raqueti
% Bruna Silveira Pavlack
% João Pedro Canisso Valese Norenberg
% Luccas Pereira Miguel
% Uncertainty Quantification
% Model calibration >> UQLab_Bayesian
%  - myMASE.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This spcrip returns the sum of the mean absolute scaled error (MASE) of
% internal variables outputs.

function epsilon = myMASE(X,Parameters)

% 1. Retrieve the static configuration parameters:
t = Parameters.t;           % Time vector;
dataA = Parameters.dataA;   % Data of viscoelastic internal variables;
dataC = Parameters.dataC;   % Data of right Cauchy-Green;
Ni = size(dataA,1);         % Number of time points;
Nj = size(dataA,2);         % Number of selected data;
Nk = size(X,1);             % Number of samples.

% 2. Calculate the model response on each sample:
epsilon = zeros(Nk,Nj);     % Pre-allocating for MASE values;
C = cell(1,Nj);             % Pre-allocating for input C(t);
IC = zeros(2,Nj);           % Pre-allocating for initial condition;
den = zeros(1,Nj);          % MASE denomiator;
for j = 1:Nj
    % Assign input C(t) and initial condition:
    C{1,j} = griddedInterpolant(t,dataC(:,j));
    IC(1,j) = dataA(1,j);
    % Compute MASE denominator:
    den(1,j) = sum(abs(dataA(2:end,j)-dataA(1:end-1,j)),1)/(Ni-1);
end
parfor k = 1:Nk
    for j = 1:Nj
        % Solve reduced-order model:
        [~,y] = ode45(@(t,Y) myode(t,Y,C{j},X(k,:)),t,IC(:,j));
        % Evaluate MASE:
        epsilon(k,j) = sum(abs((dataA(:,j)-y(:,1)))/den(j));
    end
end

% Return MASE
epsilon = sum(epsilon/Ni,2);

end