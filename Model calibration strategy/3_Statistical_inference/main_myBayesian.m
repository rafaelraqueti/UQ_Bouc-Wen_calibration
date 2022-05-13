%%=======================================================================%%
%  Author: Rafael da Silva Raqueti
%  Advisor: Samuel da Silva
%  On the Calibration of Reduced-Order Models to Describe the 
%  Viscoelasticity in Steady-State Rolling Tires 
%  Methodology >> 3_Statistical_inference
%    main_myBayesian.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  main_myBayesian.m
%    This is the main file to perform a statistical inference through
%    Bayesian inference.

%% 0 - INITIALIZATION

clearvars
rng('shuffle','twister')
% rng(100,'twister')

%__________________________________________________________________________
%% 1 - MICHELIN DATASET PRE-PROCESSING


%__________________________________________________________________________
%% 2 - PRIOR DISTRIBUTION OF MODEL PARAMETERS

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

% Discrepancy $\sigma^{2}$
X5min = +9.9e-04;
X5max = +1.0e-02;
X50   = +9.9e-04;

% Set bound vectors:
Xmin = [X1min, X2min, X3min, X4min, X5min]; % minimum values;
Xmax = [X1max, X2max, X3max, X4max, X5max]; % maximum values.
XCE = [XCE, X50];

%__________________________________________________________________________
%% 3 - SETTING RANDOM WALK METROPOLIS ALGORITHM PARAMETERS

% Metropolis-Hastings static configuration parameters:
NMC = 1e+03;    % Number of MC samples;
sigma = 0.04;   % Random walk step size.
% Uniform distribution parameters and random walk step size are stored in
% the following structure array. The generated candidate must be in the
% setted bounds:
MHParameters = struct('Xmin',Xmin,'Xmax',Xmax,'sigma',sigma);
% Parameters is passed to getcandidate.m

% Preallocate NMC samples and Likelihood function values:
samples = zeros(NMC,length(Xmin));
Psamples = zeros(NMC,1);

%__________________________________________________________________________
%% 4 - RANDOM WALK METROPOLIS ALGORITHM

% Preallocate NMCx(Ns,Nout) computed responses. These values can be
% used to compute the convergence of the Markov chain.
Y = zeros(NMC,Ni,Nj);

% Get a first candidate and its Likelihood function value:
X = XCE;
[YX,PX] = myLikelihood(X,Parameters);

% For now it is the maximum likelihood value:
Pmax = PX;
XMLE = X;

% At first, the number of accepted samples is null:
Naccept = 0;

% MONTE CARLO SIMULATION:
tic
for k = 1:NMC
    
    % Get a candidate and its Likelihood function value:
    Xcd = getcandidate(X,MHParameters);
    [Ycd,PXcd] = myLikelihood(Xcd,Parameters);
    
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
    
    % Print its values for monitoring:
    fprintf('PX = %.5f\n', PX);
    fprintf('Acceptance rate: %f\n', Naccept/k); % 15-50 % --> Adjust random walk step size.
    fprintf('%.5f completed\n', k/NMC);
end
toc

return
%% 5 - POST-PROCESSING

%% 5.1 - ROM PARAMETERS AND DISCREPANCY VARIANCE

parstr = {'$c$','$k$','$\alpha$','$\nu$','$\sigma_{\mbox{\boldmath $\epsilon$}}^{2}$'};

%% 5.1.1 - TRACE PLOTS

for par = 1:length(Xmin)
    figure
    plot(samples(:,par))
    grid on
    xlabel('Samples','Interpreter','latex')
    ylabel(parstr{par},'Interpreter','latex')
    set(gca,'TickLabelInterpreter','latex','fontsize',12)
    set(gcf,'Position',2/3*[100 100 448 268.8])
end

%% 5.1.2 - DEFINE BURN-IN AND GET SAMPLES

burnin = 0.0*NMC;            % 0% burnin samples;
x = samples(burnin+1:end,:); % Get random variables samples;
Ns = size(x,1);              % Get number of random variariable semples.

%% 5.1.3 - PLOT SAMPLES

for par = 1:length(Xmin)
    xmean = mean(x(:,par)); % Compute sample mean;
    xstd = std(x(:,par)); % Compute sample standard deviation;
    Pc = 95; % Define 95th percentile
    r_p = 0.5*(100+Pc); x_upp = prctile(x(:,par),r_p); % Get upper;
    r_m = 0.5*(100-Pc); x_low = prctile(x(:,par),r_m); % Get lower percentile.
    
    figure
    h(1) = plot(x(:,par),'oc');
    hold on
    h(2) = plot([1 Ns],[xmean xmean],'--k','linewidth',1.4);
    h(3) = plot([1 Ns],[xmean-xstd xmean-xstd],'-.b','linewidth',1.4);
    h(4) = plot([1 Ns],[x_low x_low],':m','linewidth',1.4);
    plot([1 Ns],[xmean+xstd xmean+xstd],'-.b','linewidth',1.4)
    plot([1 Ns],[x_upp x_upp],':m','linewidth',1.4)
    grid on
    xlabel('Samples','Interpreter','latex')
    ylabel(parstr{par},'Interpreter','latex')
    set(gca,'TickLabelInterpreter','latex','fontsize',12)
    legstr = {'Samples','Mean','Mean $\pm$ std','$95\%$ interval'};
    legend(h([1,2,3,4]),legstr,'Interpreter','latex')
    set(gcf,'Position',[100 100 448 336])
end

%% 5.1.4 - CONVERGENCE STUDY: 1ST MOMENT ESTIMATOR PLOT

for par = 1:length(Xmin)
    xcmean = cumsum(x(:,par))'./(1:Ns);
    figure
    h(1) = plot(x(:,par));
    hold on
    h(2) = plot(xcmean,'--k','linewidth',1.0);
    h(3) = plot(Xmin(par)*ones(1,Ns),'--r','linewidth',1.0);
    plot(Xmax(par)*ones(1,Ns),'--r','linewidth',1.0)
    grid on
    xlabel('Samples','Interpreter','latex')
    ylabel(parstr{par},'Interpreter','latex')
    set(gca,'TickLabelInterpreter','latex','fontsize',9)
    legstr = {'Samples','Mean','Support'};
    legend(h([1,2,3]),legstr,'Interpreter','latex')
    title('1st central moment','Interpreter','latex')
    set(gcf,'Position',2/3*[100 100 448 336])
end

%% 5.1.5 - CONVERGENCE STUDY: 2ND MOMENT ESTIMATOR PLOT

for par = 1:length(Xmin)
    % Second moment estimator
    xcmean = cumsum(x(:,par))'./(1:Ns);
    x2 = cumsum(x(:,par).^2)'./(1:Ns);
    xdev = sqrt(x2-xcmean.^2);
    figure
    h(1) = plot(x(:,par));
    hold on
    yyaxis left
    %h(2) = plot(xdev,'-k','linewidth',1.0);
    h(2) = plot(xcmean,'--k','linewidth',1.0);
    h(4) = plot(Xmin(par)*ones(1,Ns),'--r','linewidth',1.0);
    plot(Xmax(par)*ones(1,Ns),'--r','linewidth',1.0)
    ylabel(parstr{par},'Interpreter','latex')
    yyaxis right
    h(3) = plot(x2,'-k','linewidth',1.0);
    grid on
    xlabel('Samples','Interpreter','latex')
    set(gca,'TickLabelInterpreter','latex','fontsize',9)
    legstr = {'Chain 1','Mean','Variance','Support'};
    legend(h([1,2,3,4]),legstr,'Interpreter','latex')
    %title('2nd central moment','Interpreter','latex')
    %title('Standard deviation (SD)','Interpreter','latex')
    %set(gcf,'Position',2/3*[100 100 448 336])
    set(gcf,'Position',2/3*[100 100 1.5*448 336])
end

%% Stability

for par = 1:length(Xmin)
    figure
    h(1) = plot(chain_1(:,par),'-','linewidth',1.0);
    hold on
    h(2) = plot(chain_2(:,par),'-.','linewidth',1.0);
    h(3) = plot(chain_3(:,par),':','linewidth',1.0);
    h(4) = plot(Xmax(par)*ones(1,NMC),'--r','linewidth',1.0);
    %plot(Xmax(par)*ones(1,NMC),'--r','linewidth',1.0);    
    ylabel(parstr{par},'Interpreter','latex')
    grid on
    ylim([Xmin(par)-(Xmax(par)-Xmin(par))*(0.1),Xmax(par)+(Xmax(par)-Xmin(par))*(0.1)])
    xlabel('Samples','Interpreter','latex')
    set(gca,'TickLabelInterpreter','latex','fontsize',9)
    legstr = {'Chain 1','Chain 2','Chain 3','Support'};
    legend(h([1,2,3,4]),legstr,'Interpreter','latex')
    set(gcf,'Position',2/3*[100 100 1.5*448 336])
end

%% 5.1.6 - PLOT PARAMETERS DISTRIBUTIONS

Nbins = round(0.18*sqrt(Ns));
Nksd = round(0.18*Ns);
for par = 1:length(Xmin)
    a = Xmin(par); b = Xmax(par);
    [x_bins,x_freq,x_area] = randvar_pdf(x(:,par),Nbins);
    [x_ksd,x_supp1] = ksdensity(x(:,par),'Support',[a b]);
    [x_cdf,x_supp2] = ecdf(x(:,par));

    figure
    h(1) = bar(x_bins,x_freq,1.0);
    hold on
    yyaxis left
    h(2) = plot(x_supp1,x_ksd,'-k','LineWidth',1.0);
    h(3) = line([a b],[1/(b-a) 1/(b-a)],'Color','r','LineStyle','--','LineWidth',1.0);
    line([a a],[0 1/(b-a)],'Color','r','LineStyle','--','LineWidth',1.0)
    line([b b],[0 1/(b-a)],'Color','r','LineStyle','--','LineWidth',1.0)
    xlim([Xmin(par)-(Xmax(par)-Xmin(par))*(0.1),Xmax(par)+(Xmax(par)-Xmin(par))*(0.1)])
    xlabel(parstr{par},'Interpreter','latex')
    ylabel('Density','Interpreter','latex')
    yyaxis right
    h(4) = plot(x_supp2,x_cdf,':+k','LineWidth',1.0,'MarkerIndices',1:Ns/50:length(x_supp2));
    grid on
    ylabel('Cumulative density','Interpreter','latex')
    legend(h([3,1,2,4]),'Prior','Chain 1','EPDF','ECDF','Interpreter','latex')
    set(gca,'TickLabelInterpreter','latex','fontsize',9)
    set(gcf,'Position',3/4*[100 100 448 336])
end


%% 5.2 - ERROR MEASURE

epsilonstr = {'$\bar{E}(\mbox{\boldmath $\theta$})$'};

%% 5.2.1 - EVALUATE ERROR MEASURE

epsilon = zeros(size(x,1),Nj);
den = zeros(1,Nj);
parfor j = 1:Nj
    den(1,j) = sum(abs(tempA(2:end,j)-tempA(1:end-1,j)),1)/(Ni-1);
end
parfor k = 1:Ns
    for j = 1:Nj
        epsilon(k,j) = sum(abs((tempA(:,j)-Y(burnin+k,:,j)'))/den(j));
    end
end
epsilon = sum(epsilon/Ni,2)/Nj;

%% 5.2.2 - PLOT ERROR MEASURE SAMPLES

epsilonmean = mean(epsilon); % Compute sample mean;
epsilonstd = std(epsilon);   % Compute sample standard deviation;
Pc = 95; % Define 95th percentile
r_p = 0.5*(100+Pc); epsilon_upp = prctile(epsilon,r_p); % Get upper;
r_m = 0.5*(100-Pc); epsilon_low = prctile(epsilon,r_m); % Get lower percentile.

figure
h(1) = plot(epsilon,'oc');
hold on
h(2) = plot([1 Ns],[epsilonmean epsilonmean],'--k','linewidth',1.4);
h(3) = plot([1 Ns],[epsilonmean-epsilonstd epsilonmean-epsilonstd],'-.b','linewidth',1.4);
h(4) = plot([1 Ns],[epsilon_low epsilon_low],':m','linewidth',1.4);
plot([1 Ns],[epsilonmean+epsilonstd epsilonmean+epsilonstd],'-.b','linewidth',1.4)
plot([1 Ns],[epsilon_upp epsilon_upp],':m','linewidth',1.4)
grid on
xlabel('Samples','Interpreter','latex')
ylabel(epsilonstr,'Interpreter','latex')
set(gca,'TickLabelInterpreter','latex','fontsize',12)
legstr = {'Samples','Mean','Mean $\pm$ std','$95\%$ interval'};
legend(h([1,2,3,4]),legstr,'Interpreter','latex')
set(gcf,'Position',[100 100 448 336])

%% 5.2.3 - CONVERGENCE STUDY: 1ST MOMENT ESTIMATOR PLOT

epsiloncmean = cumsum(epsilon)'./(1:Ns);
figure
h(1) = plot(epsilon);
hold on
h(2) = plot(epsiloncmean,'--k','linewidth',1.0);
grid on
xlabel('Samples','Interpreter','latex')
ylabel(epsilonstr,'Interpreter','latex')
set(gca,'TickLabelInterpreter','latex','fontsize',9)
legstr = {'Samples','Mean'};
legend(h([1,2]),legstr,'Interpreter','latex')
set(gcf,'Position',2/3*[100 100 448 336])

%% 5.2.4 - CONVERGENCE STUDY: 2ND MOMENT ESTIMATOR PLOT

epsilon2 = cumsum(epsilon.^2)'./(1:Ns);
figure
h(1) = plot(epsilon);
hold on
h(2) = plot(epsilon2,'-k','linewidth',1.0);
grid on
xlabel('Samples','Interpreter','latex')
ylabel(epsilonstr,'Interpreter','latex')
set(gca,'TickLabelInterpreter','latex','fontsize',9)
legstr = {'Samples','Variance'};
legend(h([1,2]),legstr,'Interpreter','latex')
set(gcf,'Position',2/3*[100 100 448 336])

%%

% Second moment estimator
ecmean = cumsum(epsilon)'./(1:Ns);
e2 = cumsum(epsilon.^2)'./(1:Ns);
edev = sqrt(e2-ecmean.^2);
figure
h(1) = plot(epsilon);
hold on
yyaxis left
%h(2) = plot(xdev,'-k','linewidth',1.0);
h(2) = plot(ecmean,'--k','linewidth',1.0);
ylabel(epsilonstr,'Interpreter','latex')
yyaxis right
h(3) = plot(e2,'-k','linewidth',1.0);
grid on
xlabel('Samples','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex','fontsize',9)
legstr = {'Error','Mean','Variance'};
legend(h([1,2,3]),legstr,'Interpreter','latex')
%title('2nd central moment','Interpreter','latex')
%title('Standard deviation (SD)','Interpreter','latex')
%set(gcf,'Position',2/3*[100 100 448 336])
set(gcf,'Position',2/3*[100 100 1.5*448 336])

%% 5.2.5 - PLOT ERROR MEASURE DISTRIBUTION

Nbins = round(0.18*sqrt(Ns));
Nksd = round(0.18*Ns);

[epsilon_bins,epsilon_freq,epsilon_area] = randvar_pdf(epsilon,Nbins);
[epsilon_ksd,epsilon_supp1] = ksdensity(epsilon,'Support',[0 Inf]);
[epsilon_cdf,epsilon_supp2] = ecdf(epsilon);

figure
h(1) = bar(epsilon_bins,epsilon_freq,1.0);
hold on
yyaxis left
h(2) = plot(epsilon_supp1,epsilon_ksd,'-k','LineWidth',1.0);
axis([0 3 0 3])
xlabel(epsilonstr,'Interpreter','latex')
ylabel('Density','Interpreter','latex')
yyaxis right
h(3) = plot(epsilon_supp2,epsilon_cdf,':+k','LineWidth',1.0,'MarkerIndices',1:Ns/50:length(epsilon_supp2));
h(4) = xline(epsilonCE,'Color','m','LineStyle','--','LineWidth',1.0);
grid on
ylabel('Cumulative density','Interpreter','latex')
legend(h([1,2,3,4]),'Error','EPDF','ECDF','$\hat{E}$','Interpreter','latex')
title('Group~1','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex','fontsize',9)
set(gcf,'Position',2/3*[100 100 448 336])


%% 5.3 PROPAGATION OF UNCERTAINTIES

Groupv = [03,03,03,03; ... % Set
          01,02,03,06];    % Component
Njv = size(Groupv,2);
Yv = {};

% Pre-allocating data of:
Cv = zeros(Ni,Njv); % Right Cauchy-Green tensor;
Av = zeros(Ni,Njv); % Viscoelastic internal variables.
% And assigning selected data:
for j = 1:Njv
    Av(:,j) = dataA{16*(Branch-1)+Groupv(1,j)}(:,Groupv(2,j));
    Cv(:,j) = dataC{16*(Branch-1)+Groupv(1,j)}(:,Groupv(2,j));
end
options = odeset('RelTol',1e-06,'AbsTol',1e-08);
parfor j = 1:Njv
    for k = 1:Ns
        % Solve reduced-order model:
        C = griddedInterpolant(t,Cv(:,j)); % Input C(t);
        IC = [Av(1,j);0];                  % Initial conditions.
        [~,y] = ode45(@(t,Y) myODE2(t,Y,C,x(k,1:end-1)),t,IC,options);
        Yv{j}(:,k) = y(:,1);
    end
end

%% 5.3.1 - 99% CONFIDENCE BOUNDS

Pc = 99;
r_p = 0.5*(100 + Pc); Yv_upp = prctile(Yv{1}',r_p);
r_m = 0.5*(100 - Pc); Yv_low = prctile(Yv{1}',r_m);

% ROM response

figure
h(1) = fill([t' fliplr(t')],[Yv_low fliplr(Yv_upp)],[175 175 175]/255);
hold on
h(2) = plot(t,mean(Yv{1}'),'--k','LineWidth',1.4);
h(3) = plot(t,Av(:,1),'-b','LineWidth',1.4);
grid on
xlim([t(1) t(end)])
xlabel('Time [s]','Interpreter','latex')
ylabel('$A(t)$','Interpreter','latex')
legend(h([1,2,3]),'$95\%$ CI','Mean','Dataset','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex','fontsize',12)
set(gcf,'Position',[100 100 448 336])

% Hysteresis loop

figure
h(1) = fill([Cv(:,1)' fliplr(Cv(:,1)')],[Yv_low fliplr(Yv_upp)],[175 175 175]/255);
hold on
h(2) = plot(Cv(:,1),mean(Yv{1}'),'--k','LineWidth',1.4);
h(3) = plot(Cv(:,1),Av(:,1),'-b','LineWidth',1.4);
grid on
xlabel('$C(t)$','Interpreter','latex')
ylabel('$A(t)$','Interpreter','latex')
legend(h([1,2,3]),'$95\%$ CI','Mean','Dataset','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex','fontsize',12)
set(gcf,'Position',[100 100 448 336])

%% 5.3.2 - 95% CONFIDENCE BOUNDS

Pc = 95;
r_p = 0.5*(100 + Pc); Yv_upp = prctile(Yv{1}',r_p);
r_m = 0.5*(100 - Pc); Yv_low = prctile(Yv{1}',r_m);

% ROM response

figure
%yyaxis left
h(1) = fill([t' fliplr(t')],[Yv_low fliplr(Yv_upp)],[175 175 175]/255);
hold on
h(2) = plot(t,mean(Yv{1}'),'--k','LineWidth',1.0);
h(3) = plot(t,Av(:,1),'-b','LineWidth',1.0);
ylabel('$A(t)$','Interpreter','latex')
%yyaxis right
%h(4) = plot(t,Cv(:,2),':r','LineWidth',1.4);
%ylabel('$C(t)$','Interpreter','latex')
grid on
xlim([t(1) t(end)])
xlabel('Time [s]','Interpreter','latex')
legend(h([1,2,3]),'$95\%$ CI','Mean','$A^{\mathrm{DS}}$','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex','fontsize',9)
set(gcf,'Position',2/3*[100 100 448 336])
%%
% Hysteresis loop

figure
h(1) = fill([Cv(:,1)' fliplr(Cv(:,1)')],[Yv_low fliplr(Yv_upp)],[175 175 175]/255);
hold on
h(2) = plot(Cv(:,1),mean(Yv{1}'),'--k','LineWidth',1.0);
h(3) = plot(Cv(:,1),Av(:,1),'-b','LineWidth',1.0);
grid on
xlabel('$C(t)$','Interpreter','latex')
ylabel('$A(t)$','Interpreter','latex')
legend(h([1,2,3]),'$95\%$ CI','Mean','Dataset','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex','fontsize',9)
set(gcf,'Position',2/3*[100 100 448 336])

%% 5.3.3 - PLOT MLE RESPONSE

for j = 1:Nj
    figure
    C = griddedInterpolant(t,tempC(:,j));
    IC = [tempA(1,j);0];
    [~,y] = ode45(@(t,Y) myODE2(t,Y,C,XMLE(1,1:end-1)),t,IC);
    h(1) = plot(t,y(:,1),'--k','LineWidth',1.4);
    hold on
    h(2) = plot(t,tempA(:,j),'-b','LineWidth',1.4);
    grid on
    xlim([t(1) t(end)])
    xlabel('Time [s]','Interpreter','latex')
    ylabel('$A(t)$','Interpreter','latex')
    legend(h([1,2]),'MLE','Dataset','Interpreter','latex')
    set(gca,'TickLabelInterpreter','latex','fontsize',12)
    set(gcf,'Position',[100 100 448 336])
end
