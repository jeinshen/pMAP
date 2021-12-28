function y=fhmvmultiply_1D(h,x,q)
% last modified: 20 Feb. 2017

n=length(h);

% reverse x
xrev=x(q:-1:1);

l=2^nextpow2(n);

fft_h=fft(h,l);
fft_xx=fft(xrev,l);

yy=ifft(fft_h.*fft_xx);

% extract the product
y=yy(q:n);