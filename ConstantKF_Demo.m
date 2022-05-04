%% Scalar Example (estimates a constant):

clc
clear
close all

%System defined as a constnat of -0.37727 volts:
clear s
s.x = -0.37727;
s.A = 1;

%process noise
s.Q = 1e-5; %variance

%voltmeter to measure the voltage itself:
s.H = 1;

%measurement error
s.R = 0.01^2; %variance
Rstd = sqrt(s.R); %random measurement noise stdev

%initial state
s.x = -0.37727;
s.P = 1;

%Generate random voltage and perform the filter operation 
tru = []; %tru voltage
for t = 1:100
    tru(end+1) = -0.37727;
    s(end).y = tru(end) + Rstd*randn; %creates a random measurement
    s(end+1) = kalmanf(s(end)); %perform a Kalman filter iteration
end

%plot measurment data
figure('Renderer', 'painters', 'Position', [200 200 1200 400]);; 
hy = plot([s(1:end-1).y],'r:','LineWidth',1); %measurments
hold on;
grid on;
hk = plot([s(2:end).x],'b--','LineWidth',1); %a-posteriori state estimates
ht = plot(tru,'g-','LineWidth',1); %true data
legend([hy hk ht], 'observations','Kalman output','true voltage')
title('Estimating a constant')
xlabel('Iter')
ylabel('State Value (Voltage)')
hold off

%% KalmanF
% updates a system state vector eestimate based upon an observation, using
% a discrete Kalman filter

% I pulled this from the textbook of Marc Asch

function s = kalmanf(s)

    if ~isfield(s,'y'); error('Observation vector missing'); end
    if ~isfield(s,'x'); s.x = nan*y; end
    if ~isfield(s,'P'); x.P = nan; end
    if ~isfield(s,'A'); x.A = eye(length(s.x)); end
    if ~isfield(s,'Q'); x.Q = zeros(length(s.x)); end
    if ~isfield(s,'R'); error('Observation covariance missing'); end
    if ~isfield(s,'H'); s.H = eye(length(x)); end
    
    if isnan(s.x)
        %initialize state estimate from first observation
        if diff(sixe(s.H))
            error('Observation matrix must be square and invertible for state auto-initialization')
        end
        s.x = inv(s.H)*x.y;
        x.P = inv(s.H)*s.R*inv(s.H');
        
    else
        
        %PRediction for state vector and covariance
        s.x = s.A*s.x;
        s.P = s.A * s.P * s.A' + s.Q;
        
        %Compute Kalman gain factor
        K = (s.P)*(s.H') * inv(s.H*s.P*s.H'+s.R);
        
        %Correction based on observation
        s.x = s.x + K*(s.y - s.H*s.x);
        s.P = s.P - K*s.H*s.P;
        
    end

end
