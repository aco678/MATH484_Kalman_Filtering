function TS = ElementToTS(X)

%This function takes an element matrix and turns that matrix into a single
%vector/time series - ignoring the Hankel matrix intermediate

%Reworked from https://www.kaggle.com/jdarcy/introducing-ssa-for-time-series-decomposition

    X_rev = flip(X,2);
    
    [L,K] = size(X_rev);
    
    TS_rev = zeros(L+K-1,1);
    
    for i = -L+1:K-1
        dg = diag(X_rev,i);
        TS_rev(i+L) = mean(dg);
    end
    
    TS = flip(TS_rev);
    
end