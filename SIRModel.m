% from this source: http://ceadserv1.nku.edu/longa/classes/mat375/mat375prep/SIR%20Model.html

% alpha = infection rate
alpha = 0.005;
% beta = recovery rate
beta = 0.05;
% gamma = resusceptible rate
gamma = 0.1;

% S = y(1)
% I = y(2)
% R = y(3)
dydt = @(t,y) [-alpha*y(1)*y(2) + gamma*y(3);
               alpha*y(1)*y(2) - beta*y(2);
               -gamma*y(3) + beta*y(2)];
    
ic = [101, 1, 0];
tspan = [0 50]; 
[t,y] = ode45(dydt,tspan,ic);
plot(t,y)
legend('Susceptible', 'Infected', 'Recovered')