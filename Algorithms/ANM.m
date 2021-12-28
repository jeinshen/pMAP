function [ x_sdp,poles_sdp,coeffs_sdp ] = ANM( obs, I, n, ox1 )
%ANM This ANM fnction is designed to solve the spectral sparse signal
%       problem with sdpt3 package. 
% 

I = sort(I)';

obs = ox1(I);

[x_sdp,poles_sdp,coeffs_sdp] = ctscs_sdpt3(obs,I,n); 


end

