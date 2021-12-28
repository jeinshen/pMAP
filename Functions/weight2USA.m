function w = weight2USA(N, m, beta)
% N is the length of the time series
% m is the length of the last m data
% to generate weight vector accoridng to the second rule of GZ16 paper
% beta = 1.01; used in GZ16 paper

w = ones(N,1);

for i = 1:N-m
    w(i) = beta^(i-1);
end

alpha = w(N-m)/(m+1);

for i=1:m
    w(N-m+i) = -alpha*i + w(N-m);
end


end

