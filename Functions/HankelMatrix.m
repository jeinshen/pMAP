function H = HankelMatrix(A)
% A is a non-Hankel matrix,
% Obtain Hankel matrix H by anti-digonally averaging A

[L, K] = size(A);
v = kappa(L, K);

idx = hankel(1:L, L:(K-1)+L);

sum_A = accumarray(idx(:), A(:));
sum_A = sum_A ./ v;

H = hankel(sum_A(1:L), sum_A(L:end));


end

