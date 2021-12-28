function v = kappa(L, K)
%
% when a vector w is put into LxK Hankel matrix
% v count the repeat of each element of w apperaing in the Hankel matrix
% Example: w = (1:7); L=3; K =5;
%          v = [1 2 3 3 3 2 1]
N = L+K-1;
if L < K
    v = [(1:L)'; L*ones(K-L,1); ((N-K):-1:1)'];
elseif L == K
    v = [(1:L)'; ((N-K):-1:1)'];
else % L > K
    v = [(1:K)'; K*ones(L-K,1); ((N-L):-1:1)'];
end


end

