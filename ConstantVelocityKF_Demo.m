%% Constant Velocity Example (estimates a constant velocity):

clc
clear
close all

%%% The idea for this situation is we have an unknown, constant velocity
%%% and would like to predict both its position and its velocity from noisy
%%% observations of the position only

% x = [x x' x'']

%The position of the train is x(t) = x_0 + x_0' * t + 1/2 x'' t

clear s

%where the time step is 0.1
time_step = 0.1;
t0 = 0;
tend = 10;

s.M = [1    time_step  time_step^2;
       0       1       time_step; 
       0      0           1];

%process noise
s.Q = 1e-3*eye(size(s.M,3)); %variance

%voltmeter to measure the voltage itself:
s.H = [1 0 0];

%measurement error
s.R = 2; %variance
Rstd = sqrt(s.R); %random measurement noise stdev
%initial state

x_t = 2;  %m
v_t = 10; %m/s
a_t = -5; %m/s^2

s.x = [10; 0.5 *x_t; %%TBD]; %Initial conditions
s.P = eye(length(s.x));

%Generate random voltage and perform the filter operation 
tru = [0]; %tru voltage
for t = [t0:time_step:tend]
    tru(end+1) = tru(end)+x_t*time_step;
    s(end).y = tru(end) + Rstd*randn; %creates a random measurement
    s(end+1) = kalmanf(s(end)); %perform a Kalman filter iteration
end

%plot measurment data
figure('Renderer', 'painters', 'Position', [200 200 1200 400]);
hy = plot([t0:time_step:tend],[s(1:end-1).y],'r:','LineWidth',1); %measurments
hold on;
grid on;
kalman_s = [s(1:end).x];
hk = plot([t0:time_step:tend],kalman_s(1,2:end),'b--','LineWidth',1); %a-posteriori state estimates
ht = plot([t0:time_step:tend],tru(1:end-1),'g-','LineWidth',1); %true data
legend([hy hk ht], 'observations','Kalman output','true position')
title('Estimaton of Position with Constant Velocity')
xlabel('Time (s)')
ylabel('Position (m)')
hold off

figure('Renderer', 'painters', 'Position', [200 200 1200 400]);
h_vel = yline(x_t,'r-','LineWidth',1); %measurments
hold on;
grid on;
kalman_s = [s(1:end).x];
hk = plot([t0:time_step:tend],kalman_s(2,2:end),'b--','LineWidth',1); %a-posteriori state estimates
legend([h_vel hk], 'True Velocity','Kalman Estimated Velocity')
title('Estimation of Constant Velocity')
xlabel('Time (s)')
ylabel('Velocity (m/s)')
hold off

figure('Renderer', 'painters', 'Position', [200 200 1200 400]);
x_ax = yline(0,'k-','LineWidth',1); %xline
hold on;
grid on;
kalman_s = [s(1:end).x];
resid = kalman_s(1,2:end) - tru(1:end-1);
hr = plot([t0:time_step:tend],resid,'b--','LineWidth',1); %a-posteriori state estimates
legend([x_ax hr], 'X axis','Residuals')
title('Position Residuals')
xlabel('Time (s)')
ylabel('Position Difference (m)')
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
        if diff(size(s.H))
            error('Observation matrix must be square and invertible for state auto-initialization')
        end
        s.x = inv(s.H)*x.y;
        x.P = inv(s.H)*s.R*inv(s.H');
        
    else
        
        %PRediction for state vector and covariance
        s.x = s.M*s.x;
        s.P = s.M * s.P * s.M' + s.Q;
        
        %Compute Kalman gain factor
        K = (s.P)*(s.H') * inv(s.H*s.P*s.H'+s.R);
        
        %Correction based on observation
        s.x = s.x + K*(s.y - s.H*s.x);
        s.P = s.P - K*s.H*s.P;
        
    end

end
