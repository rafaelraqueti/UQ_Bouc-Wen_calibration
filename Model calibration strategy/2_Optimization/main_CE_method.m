%%=======================================================================%%
%  Author: Rafael da Silva Raqueti
%  Advisor: Samuel da Silva
%  On the Calibration of Reduced-Order Models to Describe the 
%  Viscoelasticity in Steady-State Rolling Tires 
%  Methodology >> 2_Optimization
%    main_CE_method.m
%%=======================================================================%%

%  main_CE_method.m 
%    This is the main file to solve an optimization problem
%    through the Cross-Entropy Method.

% [1] Americo Cunha Jr. Enhancing the performance of a bistable energy
% harvesting device via the cross-entropy method. Nonlinear Dynamics,
% Springer Verlag, 2021, 103, pp.137-155. hal-01531845v4

%__________________________________________________________________________
%% 0 - INITIALIZATION

clearvars
rng(100,'twister')

%__________________________________________________________________________
%% 1 - MICHELIN DATASET PRE-PROCESSING


%__________________________________________________________________________
%% 2 - DEFINE BOUNDED REGION

% Parameter $c$
X1min = +0.0e+00;
X1max = +1.0e-02;
% X1max = +3.0e-02;

% Parameter $k$
X2min = +9.99e-01;
X2max = +1.001e+00;

% Parameter $\alpha$
X3min = +0.0e+00;
X3max = +1.0e+00;
% X3max = +2.0e+00;

% Parameter $\nu$
X4min = +1.0e+00;
X4max = +3.0e+00;

% Define the limits of the design variables:
Xmin = [X1min, X2min, X3min, X4min]; % Minimum;
Xmax = [X1max, X2max, X3max, X4max]; % Maximum.
Nl = length(Xmin);

%__________________________________________________________________________
%% 3 - CROSS-ENTROPY METHOD

% Step 1 - Set the Cross-Entropy method parameters:
Ns = 100;        % # of samples;
p = 0.04;        % 0 < p < 1, Ne < Ns;
Ne = ceil(p*Ns); % # of elite samples;
lmax = 5e+02;    % Maximum of iteration levels; 
tol = 1.0e-06;   % Stopping criterion;
a = 0.8;         % Smooth update schema.

% Define initial hyper-parameters values:
mu = Xmin+(Xmax-Xmin).*rand(1,Nl); % Mean;
sigma = 5*(Xmax-Xmin);               % Standard deviation.

% And set iteration level to 0:
l = 0;

% Define old hyper-parameters values:
mu0 = mu;       % Mean;
sigma0 = sigma; % Standard deviation;
Smax0 = Inf;    % Old function maximum value.

X = zeros(Ns,Nl); % Preallocating Ns random samples;
% Solve the optimization problem:
tic
while norm(sigma,Inf) > tol && l < lmax
    
    % 2. Update the level counter:
    l = l+1;
    
    % 3. Generate Ns iid samples:
    suppl = ((Xmin-mu)./sigma);
    supph = ((Xmax-mu)./sigma);
    for k = 1:Ns
        for par = 1:Nl
            X(k,par) = mu(par)+sigma(par)*trandn(suppl(par),supph(par));
        end
    end
    
    % 4. Evaluate the objective function:
    S = myMASE(X,Parameters);
    [Ssort, Isort] = sort(S,'ascend'); % Sort the results;
    
    % 5. Update estimators with the elite sample set:
    
    % I do not understand the importance of Smax in the code:
    Smax = min(Ssort(Ne));
    if Smax < Smax0
        Smax0 = Smax; % But I think it should be a very small value.
    end
    
    % New hyper-parameters values:
    mu = mean(X(Isort(1:Ne),:));    % Mean;
    sigma = std(X(Isort(1:Ne),:));  % Standard deviation.
    
    % Compute and print stopping criterion:
    cond = norm(sigma,Inf); disp(cond);
    if cond < tol
        break;
    else
        % Smooth updating scheme:
        mu = a*mu+(1-a)*mu0;            % Mean;
        sigma = a*sigma+(1-a)*sigma0;   % Standard deviation.
        % Store old hyper-parameters values:
        mu0 = mu;
        sigma0 = sigma;
    end
end
toc

%__________________________________________________________________________
%% 4 - POST-PROCESSING RESULTS

% Store optimal values:
XCE = mu;
epsilonCE = myMASE(XCE,Parameters);

% Results:
for j = 1:Nj
    % Compute responses:
    C = griddedInterpolant(t,tempC(:,j));
    IC = [tempA(1,j); 0];
    [~,Y] = ode45(@(t,Y) myODE2(t,Y,C,XCE),t,IC);
    % And plot:
    figure
    plot(t,tempA(:,j),'-b')
    hold on
    plot(t,Y(:,1),'--r')
end
