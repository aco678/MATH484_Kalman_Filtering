%%%Kalman Ensemble%%%
%take in data%
ens = zeros(n); %initialize ensemble
xbar = zeros(n);
ybar = zeros(n);
ubar = zeros(n);
R; %covariance matrix
y; %observation matrix
H; %observation model
X; %normalized anomaly
Y; %normalized anomaly
KG; %kalman gain
xpred; %predicted ensemble
M; %Forward Model
for k = 1:K
    u = sqrt(R(k)).*randn(m,1);
    for i = 1:m
        y(i,k) = y(k) + u(i);
    end
    xbar(k) = (1/m)*sum(ens(:,k));
    ubar(k) = (1/m)*sum(u(:,k));
    ybar(k) = (1/m)*sum(H(ens(:,k)));
    X(k) = (ens(:,k)-xbar(k))/(sqrt(m-1));
    Y(k) = (H(ens(:,k))-u(:,k)-ybar(k)+ubar(k))/(sqrt(m-1));
    KG(k) = X(k)*transpose(Y(k))*inv(Y(k)*transpose(Y(k)));
    for i = 1:m
        xpred(i,k) = ens(i,k)+KG(k)*(y(i,k)-H(ens(i:k)));
    end
    for i = 1:m
        ens(i,k+1) = M(xpred(i,k));
    end

    
end


   