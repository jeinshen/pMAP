%%%%%%
% This demo generates a signal composed of k complex sinusoids, and randomly
% get m samples from n equispaced samples. Two algorithm, one based on
% semidefinite programming (ctscs_sdpt) and one based on discretization and
% Lass (ctscs_dast), are used to recover the missing entries and identify
% the frequencies of the individual sinusoids. 
%
% Please refer to 
% G. Tang, B. Bhaskar, P. Shah, and B. Recht, "Compressed sensing off the
% grid" for more details of the algorithms
%
% To run the codes, one needs sdpt3 (available at
% http://www.math.nus.edu.sg/~mattohkc/sdpt3.html)
% and cvx (available at http://cvxr.com/cvx/)
%%%%%%
%% 
clear;
clc;
javaaddpath java;
n=99; % # of equispaced samples
k=5;  % # of frequencies
m=25;  % # of random samples
upsample = 32; % discretization level is upsample x n

amps = exp(2*pi*1i*rand(k,1)).*(randn(k,1).^2+.5); %random amplitudes with random phase

%generate signal with random frequencies but separated by 2/n
[signal,true_c,true_poles] = moment_vector(n,k,'random',amps,2); 

I = randperm(n); 
I = sort(I(1:m));
xVals = signal(I); %get m random samples

fprintf('Running SDP...')
T0 = clock;
%recover missing samples and identify the frequencies from the random
%samples using semidefinite programming
[x_sdp,poles_sdp,coeffs_sdp] = ctscs_sdpt3(xVals,I,n); 
fprintf(' done. %.2f s\n',etime(clock,T0));

%%
fprintf('Running DAST...')
T0 = clock; 
%recover missing samples and identify the frequencies from the random
%samples using discretization and lasso
%[x_dast,poles_dast,coeffs_dast] = ctscs_dast(xVals,I,n,upsample); 
fprintf(' done. %.2f s\n',etime(clock,T0));


fprintf('\n');
fprintf('SDP Error: %.3e\n',norm(x_sdp-signal)/norm(signal))
%fprintf('DAST Error: %.3e\n',norm(x_dast-signal)/norm(signal))

%%
figure(1)
subplot(2,1,1)
stem(true_poles,abs(true_c),'b'); 
hold on
stem(poles_sdp,abs(coeffs_sdp),'r')
hold off
legend('True Poles','Est Poles');
title('SDP')

subplot(2,1,2)
stem(true_poles,abs(true_c),'b'); 
hold on
stem(poles_dast,abs(coeffs_dast),'r')
hold off
legend('True Poles','Est Poles');
title('DAST')


