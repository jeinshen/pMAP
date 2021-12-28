function [x,c,w] = moment_vector(n,k,varargin)
% x = moment_vector(n,k,frequency_spacing,amplitude)
% Generate n dimensional moment vector with
% x_{m+1} = sum{l=1}^{k} c_l exp(i 2 pi w_l m)
% REQUIRED PARAMETERS:
%   n, the number of samples
%   k, number of frequencies
% OPTIONAL PARAMETERS:
%   frequency_spacing is 'equispaced' (default) or 'random'
%   amplitude (default 1) is either k dimensional or a scalar
%   separationFactor (default 1). The spacing is at least separationFactor/n

opt = {'equispaced', ones(k,1),1};
opt(1:length(varargin)) = varargin;
[spacing,c,sepFactor]=opt{:};
if strcmp(spacing,'equispaced')
  w  = linspace(0,2*pi,k+1).';
  w = mod(w + 2*pi*rand,2*pi); % Random Offset, off grid now
  w=w(1:k);
elseif strcmp(spacing,'random')
  w = 2*pi*spaced_frequencies(sepFactor/n,k);
  w=sort(w);
else
  error('InvalidArgument',...
  'The options for frequency spacing are [equispaced|random]')
end

x = sum( diag(c)*exp(1i * w * (0:n-1)),1).';
end
