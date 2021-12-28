function [W] = wmatrix(w, K, L)
% 
% Give a vector weight w (length(w) = N)
% transform it to equivalent Hankel weight matrix of LxK (N = L+K-1)
% W is the equivalent Hankel matrix weight of the vector weight w
% Be more accurate:
% \sum w_i (x_i - y_i)^2 = \| W o (X-Y) \|^2
% X, Y are Hankel matrices from the vector x and y respectively
% W is also a Hankel matrix 

N = length(w);
if ( (L+K-1) ~= N)
    error('Size of Hankle matrix does not conform with N: Re-input L,K');
end

if L < K
    v = [(1:L)'; L*ones(K-L,1); ((N-K):-1:1)'];
elseif L == K
    v = [(1:L)'; ((N-K):-1:1)'];
else % L > K
    v = [(1:K)'; K*ones(L-K,1); ((N-L):-1:1)'];
end

v = sqrt(v);

w = sqrt(w) .* (1 ./v);
W = hankel(w(1:L), w(L:N));

end

