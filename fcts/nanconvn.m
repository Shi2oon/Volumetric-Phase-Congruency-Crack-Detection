function c = nanconvn(a, k)
%apply a filter to an image containing NaNs
%2019.06.25

%locate the NaNs
nanID = isnan(a);

%NaNs --> 0s
a0 = a;
a0(nanID) = 0;

%nonNaNs --> 1s
a1 = zeros(size(a));
a1(~nanID) = 1;

%counting filter
k_count = zeros(size(k));
k_count(k>0) = 1;

%counting
a_nb = convn(a1,k_count,'same'); %TODO:to be checked !!!

%summation
a_k = convn(a0,k_count,'same');

%weighting
c = a_k ./ a_nb;

%
c(nanID) = NaN;

