function M = LowRankProj(A, r)
%
% M = arg min \| (A-X) \| s.t. X \in \M_r (rank(X) \le r)

[L, K] = size(A);

if ( r > min (L, K) )
    error('rank r is too big: must be no bigger than L, K')
end

 if [L, K] ~= size(A)
     error('Dimesnions do agree');
 end

[U,S,V] = svds(A, r);
M = U*S*V';

end

