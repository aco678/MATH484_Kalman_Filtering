function [TS,varargout] = SSA(F,varargin)

%This function has the capacity to take in time series data (F) and then 
%use the SSA algorithm to convert this data into its component time series 
%parts
    %INPUT:
        %F - time series data - this is required to be a column or row
        %vector
        
        %L - optional parameter to be able to choose the size of the number 
        %of lagged row vectors
        
    %OUTPUT:
        %TS - this is a cell array of all of the time series ordered by
        %influence
        
        %Optional:
        
            %Wcorr - this is a symmetric matrix with correlation values
            %between time series - if the correlation of two time series is
            %greater than 0.3, then they should be grouped togehter - 1 is
            %perfectly correlated and 0 is not correlated at all
            
            %Element - this is a cell array of all the elementary matrices
            %that composed the trajectory matrix created from the lagged
            %matrix created from the data set F

[l,w] = size(F);

if l~=1 && w~=1
    error('improper data size')
end

F = F(:);

N = length(F);    %length of data

if nargin==2
    L = varargin{1};
    if L==2
        warning('Small L - choose a different L')
    elseif L>=N/2
        warning('large L - choose a different L')
    end
elseif nargin == 1
    L = ceil(N/3);
else
    error('Not the correct number of inputs')
end

K = N - L +1;       %rows of K-lagged column vectors

X = zeros(L,K);     %Trajectory matrix

for i = 1:L
    X(i,:) = F(i:i+N-L);    
end

[U,S,V] = svd(X);
%U - orthonormal basis of time subseries - columns
%V - orthonormal basis of time subseries - rows
%S - singular values - relative importance

%% Elementary Matrices

d = rank(X);

X_elem = cell(1,d);

for i = 1:d
    X_hold = S(i,i)*U(:,i)*V(:,i)';
    X_elem{i} = X_hold;
end

%% Convert back to Timeseries

F_recon = cell(1,d);

for i= 1:L
    F_recon{i} = ElementToTS(X_elem{i});    
end

F_norms = zeros(1,d);
W_corr = eye(d);
for i= 1:d
    
    F_norms(i) = wDot(F_recon{i},F_recon{i},L).^(-1/2);
    
end

for i = 1:d
    for j = i+1:d
        
        W_corr(i,j) = abs(wDot(F_recon{i},F_recon{j},L)*F_norms(i)*F_norms(j));
        W_corr(j,i) = W_corr(i,j);
        
    end
end

%% Output Collection

TS = F_recon;
Element = X_elem;
Wcorr = W_corr;

s{1} = Wcorr;
s{2} = Element;

nout = max(nargout,1) - 1;
for k = 1:nout
    varargout{k} = s{k};
end
