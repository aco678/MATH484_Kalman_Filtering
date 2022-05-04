function FXP = wDot(TS1,TS2,L)

%This is a function to calculate weighted inner products

%Reworked from https://www.kaggle.com/jdarcy/introducing-ssa-for-time-series-decomposition

FXP = 0;
N = length(TS1);
K = N-L+1;
    for k = 1:N
        if 1<=k && k<=L
            w = k;
        elseif L+1<=k && k<=K
            w = L;
        elseif K+1<=k && k<=N
            w = N-k+1;
        end
        
        FXP = FXP + w*TS1(k)*TS2(k);
        
    end

end