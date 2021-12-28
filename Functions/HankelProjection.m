function H = HankelProjection(W, A)
%
% This code is to solve the weighted projection problem onto the subspace
% of Hankel matrices
%
%               H = argmin \| W o (A-X) \| s.t. X is Hankel
% where
% W: Weight matrix >= 0
% A: Given matrix

H = W.*A;
[L, K] = size(W);
idx = hankel(1:L, L:(K-1)+L);

sum_H = accumarray(idx(:), H(:));
sum_W = accumarray(idx(:), W(:));

sum_W(sum_W==0) = 1;

sum_H = sum_H ./ sum_W;
H = hankel(sum_H(1:L), sum_H(L:end));

end

