%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rafael da Silva Raqueti
% Bruna Silveira Pavlack
% João Pedro Canisso Valese Norenberg
% Luccas Pereira Miguel
% Uncertainty Quantification
% Model calibration >> myBayesian
%  - main_myBayesian.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This spcrip is the main file for Bayesian inference.

% Before beginning:
clearvars % Clear all variables from Workspace;
rng(100,'twister') % And control the random number generator.

%% 1 - DATASET

%% 2 - PRIOR DISTRIBUTION OF THE MODEL PARAMETERS

% These values are theta from cross-entropy method:
thetaCE = [...
    ...
    ];

% Parameter $c$
X1min = +0.0e+00;
X1max = +5.0e-03;

% Parameter $k$
X2min = +9.995e-01;
X2max = +1.0005e+00;

% Parameter $\alpha$
X3min = +0.0e+00;
X3max = +1.0e+00;

% Parameter $\nu$
X4min = +1.0e+00;
X4max = +3.0e+00;

% Discrepancy $\sigma^{2}$
X5min = +0.0e+00;
X5max = +1.0e+00;
X50   = +2.0e-02;

% Set bound vectors:
Xmin = [X1min, X2min, X3min, X4min, X5min]; % minimum values;
Xmax = [X1max, X2max, X3max, X4max, X5max]; % maximum values.
XCE = [thetaCE(1,:), X50];


%% 3 - SETTING RANDOM WALK METROPOLIS ALGORITHM PARAMETERS

% Parameters:
NMC = 1e+05;   % Number of MC samples;
sigma = 0.2;   % Random walk step size.

% Uniform distribution parameters and random walk step size are stored in
% the following structure array. Generated candidate must be in the setted
% bounds.
MHParameters = struct('Xmin',Xmin,'Xmax',Xmax,'sigma',sigma);
% Parameters is passed to getcandidate.m

% Pre-allocating for NMC samples and likelihood function values:
samples = zeros(NMC,length(Xmin));
Psamples = zeros(NMC,1);


%% 4 - RANDOM WALK METROPOLIS ALGORITHM

% Preallocating NMC x (Ns,Nout) computed responses. These values can be
% used to compute the convergence of the Markov chain.
Y = zeros(NMC,Ni,Nj);

% Get a first candidate and its likelihood function value:
X = XCE;
[YX,PX] = mylikelihoodfunc(X,Parameters);

% For now it is the maximum likelihood value;
Pmax = PX;
XMLE = X;

% And the number of accepted samples is null.
Naccept = 0;

% Monte Carlo simulation
tic
for k = 1:NMC
    
    % Get a candidate and its likelihood function value:
    Xcd = getcandidate(X,MHParameters);
    [Ycd,PXcd] = mylikelihoodfunc(Xcd,Parameters);
    
    % Compute the acceptance and generate a uniformly distributed
    % pseudorandom number:
    a = min(0, log10(PXcd)-log10(PX));
    u = log10(rand);
    % X,YX and PX get candidate values if the it is accepted:
    if (a >= u)
        X = Xcd;
        YX = Ycd;
        PX = PXcd;
        Naccept = Naccept+1;
        
        % For maximum likelihood estimation:
        if(PX > Pmax)
            XMLE = X;
        end
        
    end
    % Markov chain gets X, YX and PX. In the the case the current candidate
    % is not accepted, these values are already for the last accepted
    % candidate:
    samples(k,:) = X;
    Y(k,:,:) = YX;
    Psamples(k,1) = PX;
    
    % Just print its values for monitoring:
    fprintf('PX = %.5f\n', PX);
    fprintf('Acceptance rate: %f\n', Naccept/k); % 15-50 % --> Adjust random walk step size
    fprintf('%.5f completed\n', k/NMC);
    
end
toc

return
%% 5 - POST-PROCESSING SAMPLES

burnin = 0.5*NMC;               % 50% burnin samples;
x = samples(burnin+1:end,:);    % Get random variables samples.

%% 5.1 - CONVERGENCE STUDY: PLOT PARAMETERS MEAN

for par = 1:length(Xmin)
    xmean = zeros(1,size(x,1));
    for k = 1:size(x,1)
        xmean(k) = mean(x(1:k,par));
    end
    figure
    plot(x(:,par),'o')
    hold on
    plot(xmean,'--')
    plot(Xmin(par)*ones(1,size(x,1)),'--r')
    plot(Xmax(par)*ones(1,size(x,1)),'--r')
    
end

%% 5.2 - CONVERGENCE STUDY: PLOT PARAMETERS VARIANCE

for par = 1:length(Xmin)
    xvariance = zeros(1,size(x,1));
    for k = 1:size(x,1)
        xvariance(k) = mean(x(1:k,par).^2);
    end
    figure
    plot(xvariance)
    hold on
end

%% 5.3 - DISTRIBUTION OF THE SUM OF MASE MEASURES

epsilon = zeros(size(x,1),Nj);
den = zeros(1,Nj);
for j = 1:Nj
    den(1,j) = sum(abs(tempA(2:end,j)-tempA(1:end-1,j)),1)/(Ni-1);
end
parfor k = 1:size(x,1)
    for j = 1:Nj
        epsilon(k,j) = sum(abs((tempA(:,j)-Y(burnin+k,:,j)'))/den(j));
    end
end
% Return MASE
epsilon = sum(epsilon/Ni,2);
figure
histogram(epsilon)

%% 5.4 - PROPAGATION OF UNCERTAINTIES WITH 99% CONFIDENCE BOUNDS

tempY = zeros(Ni,NMC);
for j = 1:Nj
    figure
    for k = 1:size(x,1)
        for tau = 1:Ni
            tempY(tau,k) = Y(k,tau,j);
        end
        plot(t,tempY(:,k),'Color',[175 175 175]/255)
        hold on
    end
    grid on
    plot(t,tempA(:,j),'-b')
    xlabel('Time (s)','Interpreter','latex')
    ylabel('$A(t)$','Interpreter','latex')
    set(gca,'TickLabelInterpreter','latex','fontsize',12)
    set(gcf,'Position',[100 100 448 336])
end
