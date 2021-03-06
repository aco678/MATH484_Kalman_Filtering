%% Prep

clc;
clear;
close all;

addpath(genpath(pwd))


%% Run Options

Video = 0; %1 for true

%% Create Dummy Data

%Parameter Setting
% alpha = infection rate: S --> I
alpha = 0.005;
% beta = recovery rate: I --> R
beta = 0.05;
% gamma = resusceptible rate: R --> S
gamma = 0.03;
params = [alpha, beta, gamma];


%Solve ODE
IC = [200; 1; 0];           %Actual initial conditions
mT_step = 5;                %Size of step between measurements
measT = [0:mT_step:100];    %Timing of all measurements
[t,truY] = ode45( @(t,St) SIRS(St,params,t), measT, IC );

%Isolate Data of Interest
intY = truY(:,2);
measCov = 0.1;
errY = normrnd(1,measCov^(1/2),length(intY),1);
dumY = intY.*errY;

%Visualize Initial Data
figure;
plot(measT,truY(:,2),'b-o')
hold on;
plot(measT,dumY,'k-o')
plot(measT,truY(:,1),'g-o')
plot(measT,truY(:,3),'m-o')
legend('True I','Obs I','True S','True R')
xlabel('Time (days)')
ylabel('Individuals')
title('Initial Data')
hold off

%Variance in Data
dataSmoo = smooth(dumY);
dataVar = std(dumY-dataSmoo).^2;

%Variance with SSA
TS = SSA(dumY);
TS = TS{1}+TS{2}; %Retain 2 eigenvalues for a fairly smooth figure
dataVarTS = std(dumY-TS).^2;

%Determine smoothing to be applied
figure
plot(measT,dumY,'k-o')
hold on;
plot(measT,dataSmoo,'r-o')
plot(measT,TS,'b-o')
legend('Raw Data','Smoothed Data','SSA')
xlabel('Time (days)')
ylabel('Individuals')
title({strcat('Smoothing -> Variance:',string(dataVar));...
    strcat('SSA -> Variance:',string(dataVarTS))})
hold off;


%% Prep Ensemble Kalman Filter

%Number of States
nS = 3;                     %number of states
nE = 20;                    %number of ensembles
nM = length(measT);         %number of measurements

xbar = zeros(nS, nE);       %
ybar = zeros(nM, nE);       %predicted measurements (forecasted meas)
ubar = zeros(nM, nE);
R = dataVar;                %measurement covariance matrix - matches the size of the measurement
Q = 10;                   %process covariance matrix - matches the size of the state vector - can be a single value if the variance is the same across the popn
y = dumY; %observation matrix
H = [0  1  0] ;            %observation model
% X; %num_ensemblesormalized anum_ensemblesomaly
% Y; %num_ensemblesormalized anum_ensemblesomaly
% KG; %kalman gain
% xpred; %predicted ensemble

timestep_vec = measT; %Times with measurements


% alpha = inum_ensemblesfectionum_ensembles rate
aGuess = 0.005;
% beta = recovery rate
bGuess = 0.05;
% gamma = resusceptible rate
gGuess = 0.03;
params = [aGuess, bGuess, gGuess];

IC = [200;1;0];

%% Initialization Step

store_xFor  = NaN(3*(nM-1),nE); %Don't have a forecast on the first run
store_xAsis = NaN(3*nM,nE);
store_xTraj = NaN(100*(nM-1),3*nE); %Each ensemble will have a single matrix of time points - needs two levels?

ens = ones(3,nE).*IC;

initVar = normrnd(0,R^(1/2),nS,nE);

ens = ens + initVar;

ens = max(ens,0); %floor function

store_xAsis(1:3,:) = ens;

trajTime = [];

for k = 1:nM-1
    
    time_start = timestep_vec(k);
    time_end = timestep_vec(k+1);
    
    ens = zeros(nS, nE);
    
    anVec = store_xAsis(3*(k-1)+1:3*k,:);
    
    timeGap = linspace(time_start, time_end, 100);
    
    trajTime = [trajTime timeGap];
    
    %Forecast Step
    for i = 1:nE
        [t,Y] = ode45( @(t,St) SIRS(St,params,t), timeGap, anVec(:,i) );
        ens(:,i) = Y(end,:);%%%NOTE - check to make sure this is correct
        store_xTraj(100*(k-1)+1:100*k,3*(i-1)+1:3*i) = Y;
    end
    
    %Process Noise
%     u = NaN(nS,nE); %Initializing
%     u = normrnd(0,Q^(1/2),nS,nE);
%     
%     x_for = ens + u;
    
    x_for = ens;
    
    store_xFor(3*(k-1)+1:3*(k),:) = x_for;
    
    %Calculate Variance of States
    
%     S_var = var(x_for(1,:));
%     I_var = var(x_for(2,:));
%     R_var = var(x_for(3,:));
%     
%     Cov_Upd = diag([S_var I_var R_var]); %We aren't doing this part right
    
    Cov_Upd = cov(x_for');
    
    %Isolate Estimated Measured State
    
    y_for = x_for(2,:); %All Infected Measurements - would be associated with H normally
    
    %Update
    
    KG = Cov_Upd * H' * (H * Cov_Upd * H' + R)^(-1);
    
    x_asis = x_for + KG * ( y(k+1) - y_for);
    
    store_xAsis(3*(k)+1:3*(k+1),:) = x_asis;
    
    
end

%% Kalman Filter Trajectory Plotting

%Analysis/Update PostProcessing
EnsAsisS = store_xAsis([1:3:nM*3-2],:);
EnsAsisI = store_xAsis([2:3:nM*3-1],:);
EnsAsisR = store_xAsis([3:3:nM*3],:);

SmeanAsis = mean(EnsAsisS,2);
ImeanAsis = mean(EnsAsisI,2);
RmeanAsis = mean(EnsAsisR,2);

%Forecast PostProcessing
EnsForS = store_xFor([1:3:(nM-1)*3-2],:);
EnsForI = store_xFor([2:3:(nM-1)*3-1],:);
EnsForR = store_xFor([3:3:(nM-1)*3],:);

SmeanFor = mean(EnsForS,2);
ImeanFor = mean(EnsForI,2);
RmeanFor = mean(EnsForR,2);

%Trajectory of Ensembles PostProcessing
EnsTrajS = store_xTraj(:,[1:3:nE*3-2]);
EnsTrajI = store_xTraj(:,[2:3:nE*3-1]);
EnsTrajR = store_xTraj(:,[3:3:nE*3]);

SmeanTraj = mean(EnsTrajS,2);
ImeanTraj = mean(EnsTrajI,2);
RmeanTraj = mean(EnsTrajR,2);

figure('Renderer', 'painters', 'Position', [700 550 800 450])
subplot(1,2,1)
plot(measT,SmeanAsis,'ro','linewidth',2);
hold on;
plot(measT(2:end),SmeanFor,'r*','MarkerSize',10,'linewidth',1.5);
plot(trajTime,SmeanTraj,'r-','linewidth',2);
plot(trajTime,EnsTrajS,'b-');
xlabel('Time (Days)');
ylabel('Individuals');
title('Susceptible');
legend('Mean Update','Mean Forecast','Mean Trajectory','Ensembles')
set(gca,'FontSiz',16)
hold off

subplot(1,2,2)
plot(measT,RmeanAsis,'ro','linewidth',2);
hold on;
plot(measT(2:end),RmeanFor,'r*','MarkerSize',10,'linewidth',1.5);
plot(trajTime,RmeanTraj,'r-','linewidth',2);
plot(trajTime,EnsTrajR,'b-');
xlabel('Time (Days)');
ylabel('Individuals');
title('Recovered');
% legend('Mean Update','Mean Trajectory','Ensembles')
set(gca,'FontSiz',16)
hold off

sgtitle("Unobserved States")

figure('Renderer', 'painters', 'Position', [200 550 500 450])
plot(measT,y,'k-o','linewidth',3);
hold on;
plot(measT,ImeanAsis,'ro','linewidth',2);
plot(measT(2:end),ImeanFor,'r*','MarkerSize',10,'linewidth',1.5);
plot(trajTime,ImeanTraj,'r-','linewidth',2);
plot(trajTime,EnsTrajI,'b-');
xlabel('Time (Days)');
ylabel('Individuals');
title('Infected');
legend('Data','Mean Update','Mean Forecast','Mean Trajectory','Ensembles')
set(gca,'FontSiz',16)
hold off

%% Residuals to True Values
% only applicable when we have a real solution

figure('Renderer', 'painters', 'Position', [200 50 1200 450])
subplot(1,3,1)
plot(measT, SmeanAsis - truY(:,1),'rx')
hold on
xL = xline(0);
xL.Color = "black";
xlabel('Time (min)');
ylabel('Residual (Participants)');
title('Suceptible')
set(gca,'FontSiz',16)


subplot(1,3,2)
plot(measT, ImeanAsis - truY(:,2),'rx')
hold on
xL = xline(0);
xL.Color = "black";
xlabel('Time (min)');
ylabel('Residual (Participants)');
title('Infected')
set(gca,'FontSiz',16)

subplot(1,3,3)
plot(measT, RmeanAsis - truY(:,3),'rx')
hold on
xL = xline(0);
xL.Color = "black";
xlabel('Time (min)');
ylabel('Residual (Participants)');
title('Recovered')

sgtitle('Residuals of State Compared to True State')
set(gca,'FontSiz',16)


%% Video Writing

if Video

    v = VideoWriter('KM_algorithm.avi');
    
    v.FrameRate = 5;
    
    open(v);

    figure('Renderer', 'painters', 'Position', [100 100 1600 900])

    stS = 3; %starting step of figure
    trajSS = length(trajTime)/(nM-1); %Step size for a trajectory

    trajSt = (stS-1)*trajSS + 1; %index of the start of the focused trajectory
    trajEn = (stS)*trajSS + 1;
    
    forSteps = 50;

    %Initial Starting Point - 3 seconds
    for k = 1:15

        subplot(1,2,1)
        plot(measT(1:stS),y(1:stS),'k-o','linewidth',3);
        hold on;
        plot(measT(1:stS),ImeanAsis(1:stS),'ro','linewidth',2);
        plot(measT(2:stS),ImeanFor(1:stS-1),'r*','MarkerSize',10,'linewidth',1.5);
        plot(trajTime(1:trajSt),ImeanTraj(1:trajSt,:),'r-','linewidth',2);
        plot(trajTime(1:trajSt),EnsTrajI(1:trajSt,:),'b-');
        xlabel('Time (Days)');
        ylabel('Individuals');
        title('Observed State: Infected');
        xlim([0,30]);
        legend('Data','Mean Update','Mean Forecast','Mean Trajectory','Ensembles')
        set(gca,'FontSiz',16)

        subplot(1,2,2)
        plot(measT(1:stS),RmeanAsis(1:stS),'ro','linewidth',2);
        hold on;
        plot(measT(2:stS),RmeanFor(1:stS-1),'r*','MarkerSize',10,'linewidth',1.5);
        plot(trajTime(1:trajSt),RmeanTraj(1:trajSt,:),'r-','linewidth',2);
        plot(trajTime(1:trajSt),EnsTrajR(1:trajSt,:),'b-');
        xlabel('Time (Days)');
        ylabel('Individuals');
        title('Unobserved State: Suceptible');
        xlim([0,30]);
        lgd = legend('Mean Update','Mean Forecast','Mean Trajectory','Ensembles');
        lgd.Location = "southeast";

        sgtitle({'Kalman Filter Algorithm';'Step: Preparing Next Step'},'FontSiz',20)

        set(gca,'FontSiz',16)
        hold off

        frame = getframe(gcf);
        writeVideo(v,frame);

    end

    %Propagate Ensembles - 10 seconds
    for k = 1:trajSS-1
        
        trajES = trajSt + k; %New Endpoint

        subplot(1,2,1)
        plot(measT(1:stS),y(1:stS),'k-o','linewidth',3);
        hold on;
        plot(measT(1:stS),ImeanAsis(1:stS),'ro','linewidth',2);
        plot(measT(2:stS),ImeanFor(1:stS-1),'r*','MarkerSize',10,'linewidth',1.5);
        plot(trajTime(1:trajES),ImeanTraj(1:trajES,:),'r-','linewidth',2);
        plot(trajTime(1:trajES),EnsTrajI(1:trajES,:),'b-');
        xlabel('Time (Days)');
        ylabel('Individuals');
        title('Observed State: Infected');
        xlim([0,30]);
        legend('Data','Mean Update','Mean Forecast','Mean Trajectory','Ensembles')
        set(gca,'FontSiz',16)

        subplot(1,2,2)
        plot(measT(1:stS),RmeanAsis(1:stS),'ro','linewidth',2);
        hold on;
        plot(measT(2:stS),RmeanFor(1:stS-1),'r*','MarkerSize',10,'linewidth',1.5);
        plot(trajTime(1:trajES),RmeanTraj(1:trajES,:),'r-','linewidth',2);
        plot(trajTime(1:trajES),EnsTrajR(1:trajES,:),'b-');
        xlabel('Time (Days)');
        ylabel('Individuals');
        title('Unobserved State: Suceptible');
        xlim([0,30]);
        lgd = legend('Mean Update','Mean Forecast','Mean Trajectory','Ensembles');
        lgd.Location = "southeast";

        sgtitle({'Kalman Filter Algorithm';'Step: Forecasting'},'FontSiz',20)

        set(gca,'FontSiz',16)
        hold off

        frame = getframe(gcf);
        writeVideo(v,frame);

    end

    %Show Forecast State - 3 seconds
    for k = 1:15

        subplot(1,2,1)
        plot(measT(1:stS+1),y(1:stS+1),'k-o','linewidth',3);
        hold on;
        plot(measT(1:stS),ImeanAsis(1:stS),'ro','linewidth',2);
        plot(measT(2:stS+1),ImeanFor(1:stS),'r*','MarkerSize',10,'linewidth',1.5);
        plot(trajTime(1:trajEn),ImeanTraj(1:trajEn,:),'r-','linewidth',2);
        plot(trajTime(1:trajEn),EnsTrajI(1:trajEn,:),'b-');
        xlabel('Time (Days)');
        ylabel('Individuals');
        title('Observed State: Infected');
        xlim([0,30]);
        legend('Data','Mean Update','Mean Forecast','Mean Trajectory','Ensembles')
        set(gca,'FontSiz',16)

        subplot(1,2,2)
        plot(measT(1:stS+1),RmeanAsis(1:stS+1),'ro','linewidth',2);
        hold on;
        plot(measT(2:stS+1),RmeanFor(1:stS),'r*','MarkerSize',10,'linewidth',1.5);
        plot(trajTime(1:trajEn),RmeanTraj(1:trajEn,:),'r-','linewidth',2);
        plot(trajTime(1:trajEn),EnsTrajR(1:trajEn,:),'b-');
        xlabel('Time (Days)');
        ylabel('Individuals');
        title('Unobserved State: Suceptible');
        xlim([0,30]);
        lgd = legend('Mean Update','Mean Forecast','Mean Trajectory','Ensembles');
        lgd.Location = "southeast";

        sgtitle({'Kalman Filter Algorithm';'Step: Forecasted Results'},'FontSiz',20)

        set(gca,'FontSiz',16)
        hold off

        frame = getframe(gcf);
        writeVideo(v,frame);

    end

    %Results of Update - 3 seconds
    for k = 1:15

        subplot(1,2,1)
        plot(measT(1:stS+1),y(1:stS+1),'k-o','linewidth',3);
        hold on;
        plot(measT(1:stS+1),ImeanAsis(1:stS+1),'ro','linewidth',2);
        plot(measT(2:stS+1),ImeanFor(1:stS),'r*','MarkerSize',10,'linewidth',1.5);
        plot(trajTime(1:trajEn),ImeanTraj(1:trajEn,:),'r-','linewidth',2);
        plot(trajTime(1:trajEn),EnsTrajI(1:trajEn,:),'b-');
        xlabel('Time (Days)');
        ylabel('Individuals');
        title('Observed State: Infected');
        xlim([0,30]);
        legend('Data','Mean Update','Mean Forecast','Mean Trajectory','Ensembles')
        set(gca,'FontSiz',16)

        subplot(1,2,2)
        plot(measT(1:stS+1),RmeanAsis(1:stS+1),'ro','linewidth',2);
        hold on;
        plot(measT(2:stS+1),RmeanFor(1:stS),'r*','MarkerSize',10,'linewidth',1.5);
        plot(trajTime(1:trajEn),RmeanTraj(1:trajEn,:),'r-','linewidth',2);
        plot(trajTime(1:trajEn),EnsTrajR(1:trajEn,:),'b-');
        xlabel('Time (Days)');
        ylabel('Individuals');
        title('Unobserved State: Suceptible');
        xlim([0,30]);
        lgd = legend('Mean Update','Mean Forecast','Mean Trajectory','Ensembles');
        lgd.Location = "southeast";

        sgtitle({'Kalman Filter Algorithm';'Step: Updating'},'FontSiz',20)

        set(gca,'FontSiz',16)
        hold off

        frame = getframe(gcf);
        writeVideo(v,frame);

    end

    %Show Restart of Propagation - 5 seconds
    for k = 1:floor(trajSS/2)

        trajES = trajEn + k; %New Endpoint

        subplot(1,2,1)
        plot(measT(1:stS+1),y(1:stS+1),'k-o','linewidth',3);
        hold on;
        plot(measT(1:stS+1),ImeanAsis(1:stS+1),'ro','linewidth',2);
        plot(measT(2:stS+1),ImeanFor(1:stS),'r*','MarkerSize',10,'linewidth',1.5);
        plot(trajTime(1:trajES),ImeanTraj(1:trajES,:),'r-','linewidth',2);
        plot(trajTime(1:trajES),EnsTrajI(1:trajES,:),'b-');
        xlabel('Time (Days)');
        ylabel('Individuals');
        title('Observed State: Infected');
        xlim([0,30]);
        legend('Data','Mean Update','Mean Forecast','Mean Trajectory','Ensembles')
        set(gca,'FontSiz',16)

        subplot(1,2,2)
        plot(measT(1:stS+1),RmeanAsis(1:stS+1),'ro','linewidth',2);
        hold on;
        plot(measT(2:stS+1),RmeanFor(1:stS),'r*','MarkerSize',10,'linewidth',1.5);
        plot(trajTime(1:trajES),RmeanTraj(1:trajES,:),'r-','linewidth',2);
        plot(trajTime(1:trajES),EnsTrajR(1:trajES,:),'b-');
        xlabel('Time (Days)');
        ylabel('Individuals');
        title('Unobserved State: Suceptible');
        xlim([0,30]);
        lgd = legend('Mean Update','Mean Forecast','Mean Trajectory','Ensembles');
        lgd.Location = "southeast";

        sgtitle({'Kalman Filter Algorithm';'Step: Repeat'},'FontSiz',20)

        set(gca,'FontSiz',16)
        hold off

        frame = getframe(gcf);
        writeVideo(v,frame);

    end


    close(v);
    
end


%% Function - SIRS

function dYdt = SIRS(States,params,t)

    S = States(1);
    I = States(2);
    R = States(3);
    
    alpha = params(1); %contagion rate
    gamma = params(3); %recovery rate
    beta = params(2); %desensitisation rate
    
    dSdt = -alpha*S*I + gamma*R;
    dIdt =  alpha*S*I - beta*I;
    dRdt = -gamma*R + beta*I;
    
    dYdt = [dSdt; dIdt; dRdt];
end